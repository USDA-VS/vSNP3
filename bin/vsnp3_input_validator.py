#!/usr/bin/env python

from vsnp3_version import __version__

import os
import sys
import re
import gzip
import time
from collections import defaultdict, Counter
from concurrent import futures
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set


@dataclass
class ValidationError:
    """Container for validation error details"""
    file_path: str
    error_type: str
    error_message: str
    severity: str  # 'critical', 'warning', 'info'


@dataclass
class VCFHeaderScan:
    '''
    What one pass over the START of a VCF found.

    The scan stops after the column header and a few data records, so it costs
    one buffered read of a few kilobytes per file however large the file is.
    Everything step 2 needs to decide whether a file can be analysed, and
    against which reference it was called, is in the header.
    '''
    path: str
    size: int = 0
    fileformat: bool = False
    column_header: bool = False
    columns: int = 0
    contigs: Tuple[str, ...] = ()        # ##contig IDs, in header order
    data_chroms: Tuple[str, ...] = ()    # CHROM of the data records scanned
    data_lines: int = 0
    short_lines: int = 0                 # scanned data records with < 8 columns
    data_before_header: bool = False
    non_ascii: bool = False              # any byte >= 0x80 in the scanned part
    scanned_to_end: bool = False
    truncated: bool = False              # the file ends part-way through a record
    error: Optional[str] = None          # OSError text, when the file could not be read

    @property
    def reference_label(self) -> Optional[str]:
        if self.contigs:
            return self.contigs[0]
        if self.data_chroms:
            return self.data_chroms[0]
        return None


@dataclass
class VCFValidationResults:
    """Container for VCF validation results"""
    total_files: int = 0
    valid_files: List[str] = field(default_factory=list)
    corrupted_files: List[str] = field(default_factory=list)
    reference_mismatches: List[Dict[str, str]] = field(default_factory=list)
    empty_files: List[str] = field(default_factory=list)
    format_errors: List[str] = field(default_factory=list)
    permission_errors: List[str] = field(default_factory=list)
    unreadable_files: List[str] = field(default_factory=list)
    errors: List[ValidationError] = field(default_factory=list)
    # Kept for callers that read it.  Nothing is rewritten in place any more, so
    # this stays empty: a non-ASCII byte is tolerated when the file is read, and
    # is listed in non_ascii_files instead of being "fixed".
    fixed_files: List[Dict[str, str]] = field(default_factory=list)
    # Valid files that carry no data records at all: a sample with no variant
    # and no zero-coverage position, i.e. identical to the reference everywhere
    # it was covered.  Analysed, not excluded.
    header_only_files: List[str] = field(default_factory=list)
    non_ascii_files: List[str] = field(default_factory=list)
    multi_sample_files: Dict[str, int] = field(default_factory=dict)
    reference_contigs: Tuple[str, ...] = ()
    reference_file_count: int = 0
    reference_basis: str = ''
    # True when fewer than half of the readable files share the dominant
    # reference.  A minority on another reference is excluded and reported; a
    # majority lost means the run itself is misconfigured.
    reference_majority_lost: bool = False
    # file -> one sentence saying why it is not in valid_files
    reasons: Dict[str, str] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)
    seconds: float = 0.0

    @property
    def total_valid(self) -> int:
        return len(self.valid_files)

    @property
    def total_invalid(self) -> int:
        return self.total_files - self.total_valid

    @property
    def has_critical_errors(self) -> bool:
        return len(self.reference_mismatches) > 0 or len(self.corrupted_files) > 0

    @property
    def excluded_files(self) -> List[str]:
        '''Every file not in valid_files, in input order, each with a reason.'''
        return list(self.reasons)

    def reason(self, file_path: str) -> str:
        return self.reasons.get(file_path, '')

    def add_error(self, file_path: str, error_type: str, message: str, severity: str = 'warning'):
        """Add a validation error"""
        error = ValidationError(file_path, error_type, message, severity)
        self.errors.append(error)
        if file_path:
            self.reasons.setdefault(file_path, message)

        # Categorize the error
        if error_type == 'reference_mismatch':
            if file_path not in [item['file'] for item in self.reference_mismatches]:
                self.reference_mismatches.append({
                    'file': file_path,
                    'message': message
                })
        elif error_type == 'corrupted':
            if file_path not in self.corrupted_files:
                self.corrupted_files.append(file_path)
        elif error_type == 'empty':
            if file_path not in self.empty_files:
                self.empty_files.append(file_path)
        elif error_type == 'format_error':
            if file_path not in self.format_errors:
                self.format_errors.append(file_path)
        elif error_type == 'permission':
            if file_path not in self.permission_errors:
                self.permission_errors.append(file_path)
        elif error_type == 'unreadable':
            if file_path not in self.unreadable_files:
                self.unreadable_files.append(file_path)

    def add_fixed_file(self, file_path: str, original_error: str, fix_applied: str):
        """Track a file that had an error but was successfully auto-fixed"""
        self.fixed_files.append({
            'file': file_path,
            'original_error': original_error,
            'fix_applied': fix_applied,
        })


# ---------------------------------------------------------------- header scan

_CONTIG_ID = re.compile(rb'^##contig=<[^>]*?\bID=([^,>]+)')
# Stop after this many data records: enough to see which contigs the records
# are on, and the header has already been read in full by then.
_WANT_DATA_LINES = 20
# A file with no column header inside this many bytes is not a VCF.
_MAX_SCAN_BYTES = 8 * 1024 * 1024


def scan_vcf_header(path: str) -> VCFHeaderScan:
    '''
    Read the header and the first few records of one VCF.

    Binary, so the result does not depend on the process locale.  Every vSNP3
    module forces LC_ALL=C, under which Python's default text encoding is ASCII,
    and a text-mode open() then raised UnicodeDecodeError on the first byte of
    an accented letter or an en dash in a ##reference path.  That is how a
    perfectly good VCF came to be reported as "corrupted" and then rewritten
    in place to remove the offending characters.  Here a non-ASCII byte is
    noted and nothing else.
    '''
    scan = VCFHeaderScan(path=path)
    try:
        st = os.stat(path)
    except OSError as e:
        scan.error = f'{type(e).__name__}: {e}'
        return scan
    scan.size = st.st_size
    if scan.size == 0:
        scan.scanned_to_end = True
        return scan

    opener = gzip.open if path.endswith('.gz') else open
    contigs: List[str] = []
    chroms: List[str] = []
    bytes_read = 0
    try:
        with opener(path, 'rb') as f:
            for raw in f:
                bytes_read += len(raw)
                if not raw.isascii():
                    scan.non_ascii = True
                line = raw.rstrip(b'\r\n')
                if not line:
                    continue
                if line.startswith(b'##'):
                    if line.startswith(b'##fileformat=VCF'):
                        scan.fileformat = True
                    elif line.startswith(b'##contig='):
                        m = _CONTIG_ID.match(line)
                        if m:
                            contigs.append(m.group(1).decode('ascii', 'replace'))
                    if bytes_read > _MAX_SCAN_BYTES:
                        break
                    continue
                if line.startswith(b'#'):
                    if line.startswith(b'#CHROM'):
                        scan.column_header = True
                        scan.columns = line.count(b'\t') + 1
                    continue
                if not scan.column_header:
                    scan.data_before_header = True
                    break
                scan.data_lines += 1
                parts = line.split(b'\t')
                if len(parts) < 8:
                    scan.short_lines += 1
                else:
                    chrom = parts[0].decode('ascii', 'replace')
                    if chrom not in chroms:
                        chroms.append(chrom)
                if scan.data_lines >= _WANT_DATA_LINES or bytes_read > _MAX_SCAN_BYTES:
                    break
            else:
                scan.scanned_to_end = True
    except (OSError, EOFError) as e:
        scan.error = f'{type(e).__name__}: {e}'
    # The scan above reads only the start of the file, so a file cut off part-way
    # (an interrupted copy, a full disk) looked complete and was analysed with its
    # remaining records missing -- every SNP after the cut read as the reference.
    # Every VCF writer ends the last record with a newline; a file that does not
    # end with one stopped in the middle of a record.  One read of the last byte.
    # A cut that falls exactly on a line boundary cannot be seen this way.
    if not scan.error and not path.endswith('.gz') and scan.column_header:
        try:
            with open(path, 'rb') as f:
                f.seek(-1, os.SEEK_END)
                scan.truncated = f.read(1) != b'\n'
        except OSError as e:
            scan.error = f'{type(e).__name__}: {e}'
    scan.contigs = tuple(contigs)
    scan.data_chroms = tuple(chroms)
    return scan


def classify_scan(scan: VCFHeaderScan) -> Tuple[str, str]:
    '''
    (category, message) for one scan.  Category is one of
    'valid', 'empty', 'unreadable', 'permission', 'corrupted'.
    '''
    if scan.error:
        if 'Permission' in scan.error:
            return 'permission', f'VCF is not readable: {scan.error}'
        return 'unreadable', f'VCF could not be read: {scan.error}'
    if scan.size == 0:
        return 'empty', 'file is 0 bytes'
    if scan.data_before_header:
        return 'corrupted', 'not a VCF: a non-header line appears before the #CHROM column header'
    if not scan.column_header:
        if scan.scanned_to_end:
            return 'corrupted', 'not a VCF: no #CHROM column header line in the file'
        return 'corrupted', (f'not a VCF: no #CHROM column header line in the first '
                             f'{_MAX_SCAN_BYTES // (1024 * 1024)} MB')
    if scan.columns < 8:
        return 'corrupted', (f'#CHROM line has {scan.columns} column(s); a VCF needs '
                             f'at least the 8 fixed columns')
    if scan.short_lines:
        return 'corrupted', (f'{scan.short_lines} of the first {scan.data_lines} data '
                             f'record(s) have fewer than 8 tab-separated columns')
    if scan.truncated:
        return 'corrupted', ('truncated: the file ends part-way through a record (no final '
                             'newline), so records after the cut are missing; copy the file '
                             'again from its source')
    return 'valid', 'VCF format validation passed'


class InputValidator:
    """Comprehensive input validation for vSNP3"""

    def __init__(self, debug: bool = False):
        self.debug = debug
        self.validation_log = []

    def log_validation(self, message: str, level: str = 'INFO'):
        """Log validation messages"""
        log_entry = f"[{level}] {message}"
        self.validation_log.append(log_entry)
        if self.debug:
            print(log_entry)

    def validate_file_exists(self, file_path: str, file_type: str = "file") -> Tuple[bool, str]:
        """Validate file exists and is readable"""
        if not file_path:
            return False, f"No {file_type} path provided"

        if not os.path.exists(file_path):
            return False, f"{file_type} does not exist: {file_path}"

        if not os.path.isfile(file_path):
            return False, f"Path is not a file: {file_path}"

        if not os.access(file_path, os.R_OK):
            return False, f"{file_type} is not readable: {file_path}"

        # Check if file is empty
        if os.path.getsize(file_path) == 0:
            return False, f"{file_type} is empty: {file_path}"

        return True, "File validation passed"

    def validate_directory_exists(self, dir_path: str) -> Tuple[bool, str]:
        """Validate directory exists and is accessible"""
        if not dir_path:
            return False, "No directory path provided"

        expanded_path = os.path.expanduser(dir_path)
        abs_path = os.path.abspath(expanded_path)

        if not os.path.exists(abs_path):
            return False, f"Directory does not exist: {abs_path}"

        if not os.path.isdir(abs_path):
            return False, f"Path is not a directory: {abs_path}"

        if not os.access(abs_path, os.R_OK):
            return False, f"Directory is not readable: {abs_path}"

        return True, "Directory validation passed"

    def sanitize_path(self, path: str) -> Tuple[bool, str, str]:
        """Sanitize and validate file paths for security"""
        if not path:
            return False, "", "Empty path provided"

        try:
            # Expand user directory and resolve path
            expanded = os.path.expanduser(path)
            resolved = os.path.abspath(expanded)

            # Check for path traversal attempts
            if ".." in path or path.startswith("/"):
                # Allow absolute paths but warn about relative traversal
                if ".." in path:
                    return False, "", f"Path traversal detected in: {path}"

            # Check for null bytes (common injection technique)
            if "\x00" in path:
                return False, "", f"Null byte detected in path: {path}"

            # Resolve symbolic links for security check
            if os.path.islink(resolved):
                real_path = os.path.realpath(resolved)
                self.log_validation(f"Symbolic link detected: {resolved} -> {real_path}", "WARNING")
                resolved = real_path

            return True, resolved, "Path validation passed"

        except Exception as e:
            return False, "", f"Path validation failed: {str(e)}"

    def validate_fastq_format(self, file_path: str) -> Tuple[bool, str]:
        """Validate FASTQ file format"""
        is_valid, msg = self.validate_file_exists(file_path, "FASTQ")
        if not is_valid:
            return False, msg

        try:
            # Check if file is gzipped
            opener = gzip.open if file_path.endswith('.gz') else open
            mode = 'rt' if file_path.endswith('.gz') else 'r'

            # Explicit encoding: under LC_ALL=C the default is ASCII on Linux,
            # and a non-ASCII byte in a read name raised UnicodeDecodeError.
            with opener(file_path, mode, encoding='utf-8', errors='replace') as f:
                # Read first 4 lines to validate FASTQ format
                lines = []
                for i, line in enumerate(f):
                    if i >= 4:
                        break
                    lines.append(line.strip())

                if len(lines) < 4:
                    return False, f"FASTQ file too short (less than 4 lines): {file_path}"

                # Validate FASTQ structure
                if not lines[0].startswith('@'):
                    return False, f"Invalid FASTQ header (line 1): {file_path}"

                if not lines[2].startswith('+'):
                    return False, f"Invalid FASTQ separator (line 3): {file_path}"

                # Check sequence and quality lengths match
                if len(lines[1]) != len(lines[3]):
                    return False, f"Sequence and quality length mismatch: {file_path}"

                # Basic sequence validation (DNA characters)
                valid_bases = set('ACGTN')
                if not all(base.upper() in valid_bases for base in lines[1]):
                    return False, f"Invalid DNA sequence characters: {file_path}"

        except Exception as e:
            return False, f"FASTQ format validation failed: {str(e)}"

        return True, "FASTQ format validation passed"

    def validate_fasta_format(self, file_path: str) -> Tuple[bool, str]:
        """Validate FASTA file format"""
        is_valid, msg = self.validate_file_exists(file_path, "FASTA")
        if not is_valid:
            return False, msg

        try:
            with open(file_path, 'r', encoding='utf-8', errors='replace') as f:
                first_line = f.readline().strip()

                if not first_line.startswith('>'):
                    return False, f"Invalid FASTA header (must start with '>'): {file_path}"

                # Read a few more lines to validate structure
                sequence_lines = 0
                for line_num, line in enumerate(f, 2):
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        # Another header - valid multi-sequence FASTA
                        break
                    else:
                        # Sequence line
                        sequence_lines += 1
                        valid_bases = set('ACGTRYSWKMBDHVN.-')  # Include IUPAC codes
                        if not all(base.upper() in valid_bases for base in line):
                            return False, f"Invalid sequence characters on line {line_num}: {file_path}"

                        if sequence_lines > 10:  # Don't check entire file
                            break

                if sequence_lines == 0:
                    return False, f"FASTA file contains no sequence data: {file_path}"

        except Exception as e:
            return False, f"FASTA format validation failed: {str(e)}"

        return True, "FASTA format validation passed"

    # ------------------------------------------------------------------- VCF

    def validate_vcf_format(self, vcf_path: str) -> Tuple[bool, str, Optional[str]]:
        """
        (is_valid, message, reference label) for one VCF.

        A header-only file is valid: it is a sample with nothing to report,
        not a broken one.
        """
        scan = scan_vcf_header(vcf_path)
        category, message = classify_scan(scan)
        if category != 'valid':
            return False, f"VCF validation error in {os.path.basename(vcf_path)}: {message}", scan.reference_label
        return True, message, scan.reference_label

    def _is_same_reference_genome(self, ref1: str, ref2: str, min_similarity: float = 0.7) -> bool:
        """
        Check if two reference chromosome names belong to the same reference genome.

        Only used when no file carries ##contig header lines, so the reference
        has to be inferred from record CHROM names.  Handles multi-segment
        references (influenza segments, two-chromosome bacteria) whose contig
        names share a long common prefix.

        Examples:
            A/owl/CA/25-003495-001/2024_PB2 vs A/owl/CA/25-003495-001/2024_PB1 -> True (94%)
            A/owl/CA/25-003495-001/2024_PB2 vs A/mallard/CA/25-003495-001/2024_PB2 -> False (6%)
            NC_003317.1 vs NC_003318.1 -> True (73%)
        """
        if ref1 == ref2:
            return True

        common = os.path.commonprefix([ref1, ref2])
        min_len = min(len(ref1), len(ref2))

        if min_len == 0:
            return False

        return (len(common) / min_len) >= min_similarity

    @staticmethod
    def _dominant_reference(scans: List[VCFHeaderScan]) -> Tuple[frozenset, int, str]:
        '''
        The reference most of the files were called against.

        Decided from every readable file, not a random sample of 50: the old
        sample made the set of "expected" references, and therefore which files
        were flagged, differ from one run to the next on the same input.

        Identified by the ##contig ID set, which names every sequence of the
        reference FASTA whether or not the sample has a record on it.  The
        record CHROM names are the fallback for VCFs written without ##contig
        lines, where a sample with calls on one segment only looks as if it
        used a one-segment reference.
        '''
        by_contigs = Counter(frozenset(s.contigs) for s in scans if s.contigs)
        if by_contigs:
            dominant = max(by_contigs, key=lambda k: (by_contigs[k], tuple(sorted(k))))
            return dominant, by_contigs[dominant], 'contig'
        by_first = Counter(s.data_chroms[0] for s in scans if s.data_chroms)
        if by_first:
            dominant = max(by_first, key=lambda k: (by_first[k], k))
            return frozenset([dominant]), by_first[dominant], 'record'
        return frozenset(), 0, ''

    def _compatible(self, scan: VCFHeaderScan, dominant: frozenset, basis: str) -> bool:
        if not dominant:
            return True
        if scan.contigs:
            # Shares at least one sequence ID: positions on shared contigs align.
            return bool(set(scan.contigs) & dominant)
        if scan.data_chroms:
            if basis == 'contig':
                return set(scan.data_chroms) <= dominant
            anchor = next(iter(dominant))
            return all(self._is_same_reference_genome(c, anchor) for c in scan.data_chroms)
        return True                     # header only: nothing to contradict

    def validate_vcf_list(self, vcf_list: List[str]) -> VCFValidationResults:
        """
        Validate a list of VCF files and check that they share a reference.

        One short header read per file, run in a thread pool because the work
        is entirely I/O.  Results are assembled in input order so the outcome
        and the log are identical from run to run.
        """
        started = time.monotonic()
        results = VCFValidationResults()
        results.total_files = len(vcf_list)

        if not vcf_list:
            results.add_error("", "empty_list", "No VCF files provided", "critical")
            return results

        self.log_validation(f"Validating {len(vcf_list)} VCF files...")

        workers = max(1, min(32, len(vcf_list)))
        with futures.ThreadPoolExecutor(max_workers=workers) as pool:
            scans = dict(zip(vcf_list, pool.map(scan_vcf_header, vcf_list)))

        valid_scans: List[VCFHeaderScan] = []
        for vcf_file in vcf_list:
            scan = scans[vcf_file]
            category, message = classify_scan(scan)
            base = os.path.basename(vcf_file)
            if category == 'valid':
                results.valid_files.append(vcf_file)
                valid_scans.append(scan)
                if scan.data_lines == 0:
                    results.header_only_files.append(vcf_file)
                    self.log_validation(f"Header only (no records; reference-identical where covered): {base}", "INFO")
                if scan.non_ascii:
                    results.non_ascii_files.append(vcf_file)
                    self.log_validation(f"Non-ASCII bytes in header, ignored on read: {base}", "INFO")
                if scan.columns > 10:
                    results.multi_sample_files[vcf_file] = scan.columns - 9
                if not scan.fileformat:
                    self.log_validation(f"No ##fileformat line (accepted): {base}", "INFO")
                continue
            severity = 'warning' if category == 'empty' else 'critical'
            results.add_error(vcf_file, category, message, severity)
            self.log_validation(f"{category.capitalize()}: {base} - {message}",
                                "WARNING" if category == 'empty' else "ERROR")

        # Reference consistency, decided from all readable files.
        dominant, dominant_count, basis = self._dominant_reference(valid_scans)
        results.reference_contigs = tuple(sorted(dominant))
        results.reference_file_count = dominant_count
        results.reference_basis = basis
        if dominant:
            shown = ', '.join(results.reference_contigs[:4])
            if len(results.reference_contigs) > 4:
                shown += f', ... ({len(results.reference_contigs)} contigs)'
            self.log_validation(
                f"Reference: {shown} -- shared by {dominant_count} of {len(valid_scans)} "
                f"readable files (from {'##contig header lines' if basis == 'contig' else 'record CHROM names'})", "INFO")
            mismatched = [s for s in valid_scans if not self._compatible(s, dominant, basis)]
            readable = len(valid_scans)
            for scan in mismatched:
                found = ', '.join((scan.contigs or scan.data_chroms)[:4])
                message = (f"called against a different reference: file has "
                           f"{found}; {dominant_count} of {readable} files have {shown}")
                results.add_error(scan.path, "reference_mismatch", message, "critical")
                self.log_validation(f"Reference mismatch: {os.path.basename(scan.path)} uses {found}", "ERROR")
                results.valid_files.remove(scan.path)
            # The files kept must be a strict majority.  `dominant_count * 2 <
            # readable` let an even split through: three B. suis and three
            # M. bovis files kept whichever half's contig names sort higher and
            # excluded the other half, rather than aborting as a mixed set.
            kept = readable - len(mismatched)
            if mismatched and kept * 2 <= readable:
                results.reference_majority_lost = True
                results.notes.append(
                    f'Only {kept} of {readable} readable files are on the most common '
                    f'reference ({shown}), which is not a majority.')

        results.seconds = time.monotonic() - started
        self.log_validation(f"Validation complete: {results.total_valid}/{results.total_files} "
                            f"files valid in {results.seconds:.1f} s")
        return results

    def print_vcf_validation_summary(self, results: VCFValidationResults):
        """Print VCF validation results to terminal"""
        print("\n" + "="*80)
        print("📊 VCF FILE VALIDATION SUMMARY")
        print("="*80)
        print(f"Total files processed: {results.total_files}")
        print(f"Valid files: {results.total_valid}")
        if results.header_only_files:
            print(f"  - with no records (reference-identical where covered): {len(results.header_only_files)}")
        if results.non_ascii_files:
            print(f"  - with non-ASCII bytes in the header, ignored on read: {len(results.non_ascii_files)}")
        print(f"Invalid files: {results.total_invalid}")
        if results.reference_contigs:
            shown = ', '.join(results.reference_contigs[:4])
            if len(results.reference_contigs) > 4:
                shown += f', ... ({len(results.reference_contigs)} contigs)'
            print(f"Reference: {shown}")
        print(f"Validation time: {results.seconds:.1f} s")

        if results.total_invalid == 0:
            print("✅ All VCF files passed validation!")
            print("="*80)
            return

        if results.reference_mismatches:
            print(f"\n🔥 EXCLUDED: Reference Mismatch ({len(results.reference_mismatches)} files)")
            print("-" * 50)
            for mismatch in results.reference_mismatches:
                print(f"  ❌ {os.path.basename(mismatch['file'])}")
                print(f"     {mismatch['message']}")

        if results.corrupted_files:
            print(f"\n💥 EXCLUDED: Corrupted Files ({len(results.corrupted_files)} files)")
            print("-" * 50)
            for file_path in results.corrupted_files:
                print(f"  ❌ {os.path.basename(file_path)}: {results.reason(file_path)}")

        if results.empty_files:
            print(f"\n⚠️  EXCLUDED: Empty Files ({len(results.empty_files)} files)")
            print("-" * 50)
            for file_path in results.empty_files:
                print(f"  ⚠️  {os.path.basename(file_path)}: {results.reason(file_path)}")

        if results.permission_errors:
            print(f"\n🔒 EXCLUDED: Permission Issues ({len(results.permission_errors)} files)")
            print("-" * 50)
            for file_path in results.permission_errors:
                print(f"  🔒 {os.path.basename(file_path)}: {results.reason(file_path)}")

        if results.unreadable_files:
            print(f"\n❓ EXCLUDED: Unreadable Files ({len(results.unreadable_files)} files)")
            print("-" * 50)
            for file_path in results.unreadable_files:
                print(f"  ❓ {os.path.basename(file_path)}: {results.reason(file_path)}")

        for note in results.notes:
            print(f"\n### {note}")
        print("="*80)

    def write_validation_log(self, log_file: str, results: VCFValidationResults):
        """Write detailed validation log to file"""
        try:
            with open(log_file, 'w', encoding='utf-8') as f:
                f.write("VCF VALIDATION LOG\n")
                f.write("==================\n\n")
                f.write(f"Total files processed: {results.total_files}\n")
                f.write(f"Valid files: {results.total_valid}\n")
                if results.header_only_files:
                    f.write(f"  - with no records: {len(results.header_only_files)}\n")
                if results.non_ascii_files:
                    f.write(f"  - with non-ASCII bytes in the header: {len(results.non_ascii_files)}\n")
                f.write(f"Invalid files: {results.total_invalid}\n")
                if results.reference_contigs:
                    f.write(f"Reference contigs: {', '.join(results.reference_contigs)}\n")
                    f.write(f"  shared by {results.reference_file_count} readable file(s)\n")
                f.write(f"Validation time: {results.seconds:.1f} s\n\n")

                # Write validation messages
                f.write("VALIDATION MESSAGES:\n")
                f.write("-" * 40 + "\n")
                for log_msg in self.validation_log:
                    f.write(f"{log_msg}\n")
                f.write("\n")

                # Write detailed error information
                if results.errors:
                    f.write("EXCLUDED FILES:\n")
                    f.write("-" * 40 + "\n")
                    for error in results.errors:
                        f.write(f"File: {error.file_path}\n")
                        f.write(f"Type: {error.error_type}\n")
                        f.write(f"Severity: {error.severity}\n")
                        f.write(f"Message: {error.error_message}\n")
                        f.write("-" * 40 + "\n")

                if results.non_ascii_files:
                    f.write(f"\nNON-ASCII BYTES IN HEADER ({len(results.non_ascii_files)}), ignored on read, file unchanged:\n")
                    f.write("-" * 40 + "\n")
                    for vcf_file in results.non_ascii_files:
                        f.write(f"{vcf_file}\n")

                if results.header_only_files:
                    f.write(f"\nHEADER-ONLY FILES ({len(results.header_only_files)}), analysed as reference-identical:\n")
                    f.write("-" * 40 + "\n")
                    for vcf_file in results.header_only_files:
                        f.write(f"{vcf_file}\n")

                # Write valid files list
                if results.valid_files:
                    f.write(f"\nVALID FILES ({len(results.valid_files)}):\n")
                    f.write("-" * 40 + "\n")
                    for vcf_file in results.valid_files:
                        f.write(f"{vcf_file}\n")

        except Exception as e:
            print(f"Warning: Could not write validation log: {str(e)}")


def validate_excel_file(file_path: str) -> Tuple[bool, str]:
    """Validate Excel file can be opened"""
    try:
        import openpyxl
        workbook = openpyxl.load_workbook(file_path, read_only=True)
        workbook.close()
        return True, "Excel file validation passed"
    except ImportError:
        return False, "openpyxl library not available for Excel validation"
    except Exception as e:
        return False, f"Excel file validation failed: {str(e)}"


# Main validation functions for backward compatibility
def validate_file_inputs(fastq_r1=None, fastq_r2=None, fasta=None, gbk_list=None, debug=False):
    """Validate all file inputs for vSNP3 Step1"""
    validator = InputValidator(debug=debug)
    errors = []

    # Validate FASTQ files
    if fastq_r1:
        is_valid, msg = validator.validate_fastq_format(fastq_r1)
        if not is_valid:
            errors.append(f"FASTQ R1: {msg}")

    if fastq_r2:
        is_valid, msg = validator.validate_fastq_format(fastq_r2)
        if not is_valid:
            errors.append(f"FASTQ R2: {msg}")

    # Validate FASTA files
    if fasta:
        fasta_list = fasta if isinstance(fasta, list) else [fasta]
        for fasta_file in fasta_list:
            is_valid, msg = validator.validate_fasta_format(fasta_file)
            if not is_valid:
                errors.append(f"FASTA: {msg}")

    # Validate GBK files
    if gbk_list:
        for gbk_file in gbk_list:
            is_valid, msg = validator.validate_file_exists(gbk_file, "GBK")
            if not is_valid:
                errors.append(f"GBK: {msg}")

    return len(errors) == 0, errors


if __name__ == "__main__":
    # Test the validator
    import sys
    if len(sys.argv) < 2:
        print("Usage: python vsnp3_input_validator.py <vcf_file_or_directory>")
        sys.exit(1)

    validator = InputValidator(debug=True)
    test_path = sys.argv[1]

    if os.path.isfile(test_path):
        # Test single VCF file
        is_valid, msg, ref = validator.validate_vcf_format(test_path)
        print(f"File: {test_path}")
        print(f"Valid: {is_valid}")
        print(f"Message: {msg}")
        print(f"Reference: {ref}")

    elif os.path.isdir(test_path):
        # Test directory of VCF files
        import glob
        vcf_files = sorted(glob.glob(os.path.join(test_path, "*.vcf")))
        results = validator.validate_vcf_list(vcf_files)
        validator.print_vcf_validation_summary(results)
