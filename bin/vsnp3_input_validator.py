#!/usr/bin/env python

__version__ = "3.34"

import os
import sys
import re
import gzip
import pathlib
from collections import defaultdict
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, Set
import pandas as pd


@dataclass
class ValidationError:
    """Container for validation error details"""
    file_path: str
    error_type: str
    error_message: str
    severity: str  # 'critical', 'warning', 'info'


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
    fixed_files: List[Dict[str, str]] = field(default_factory=list)

    @property
    def total_valid(self) -> int:
        return len(self.valid_files)

    @property
    def total_invalid(self) -> int:
        return self.total_files - self.total_valid

    @property
    def has_critical_errors(self) -> bool:
        return len(self.reference_mismatches) > 0 or len(self.corrupted_files) > 0

    def add_error(self, file_path: str, error_type: str, message: str, severity: str = 'warning'):
        """Add a validation error"""
        error = ValidationError(file_path, error_type, message, severity)
        self.errors.append(error)

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

            with opener(file_path, mode) as f:
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
            with open(file_path, 'r') as f:
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

    def extract_vcf_reference(self, vcf_path: str) -> Tuple[Optional[str], List[str]]:
        """Extract reference information from VCF file"""
        references = set()
        errors = []

        try:
            opener = gzip.open if vcf_path.endswith('.gz') else open
            mode = 'rt' if vcf_path.endswith('.gz') else 'r'

            with opener(vcf_path, mode) as f:
                header_found = False
                data_lines = 0

                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    if not line:
                        continue

                    # Check VCF header
                    if line.startswith('##fileformat=VCF'):
                        header_found = True
                        continue

                    # Skip other header lines
                    if line.startswith('##'):
                        continue

                    # Column headers
                    if line.startswith('#CHROM'):
                        continue

                    # Data lines - extract chromosome/reference
                    if not line.startswith('#'):
                        data_lines += 1
                        parts = line.split('\t')

                        if len(parts) < 8:
                            errors.append(f"Invalid VCF format on line {line_num}: insufficient columns")
                            continue

                        chrom = parts[0]
                        references.add(chrom)

                        # Only check first few data lines for efficiency
                        if data_lines > 5:
                            break

                if not header_found:
                    errors.append("Missing VCF header (##fileformat=VCF)")

                if data_lines == 0:
                    errors.append("VCF file contains no data lines")

        except Exception as e:
            errors.append(f"Failed to read VCF file: {str(e)}")

        # Return primary reference (most common)
        if references:
            primary_ref = sorted(references)[0]  # Use first alphabetically for consistency
            return primary_ref, errors

        return None, errors

    def validate_vcf_format(self, vcf_path: str) -> Tuple[bool, str, Optional[str]]:
        """Comprehensive VCF file validation"""
        # Basic file validation
        is_valid, msg = self.validate_file_exists(vcf_path, "VCF")
        if not is_valid:
            return False, msg, None

        # Extract reference and validate format
        reference, errors = self.extract_vcf_reference(vcf_path)

        if errors:
            error_msg = f"VCF validation errors in {os.path.basename(vcf_path)}: {'; '.join(errors)}"
            return False, error_msg, reference

        if reference is None:
            return False, f"No reference found in VCF file: {vcf_path}", None

        return True, "VCF format validation passed", reference

    def _is_same_reference_genome(self, ref1: str, ref2: str, min_similarity: float = 0.7) -> bool:
        """
        Check if two reference chromosome names belong to the same reference genome.

        Handles multi-segment/multi-chromosome references (e.g., influenza segments,
        two-chromosome bacteria) where different chromosomes share a common base name
        but have chromosome-specific suffixes.

        Uses longest common prefix similarity: if two reference names share at least
        min_similarity of the shorter name as a common prefix, they are considered
        the same reference genome.

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

    def _analyze_reference_consensus(self, vcf_list: List[str], sample_size: int = 50) -> Tuple[Set[str], Dict[str, int]]:
        """
        Analyze a sample of VCF files to determine consensus reference pattern

        Returns:
        - Set of expected/acceptable references
        - Dictionary with reference counts for analysis
        """
        if not vcf_list:
            return set(), {}

        # Sample a subset for analysis to avoid processing huge datasets
        import random
        sample_files = random.sample(vcf_list, min(sample_size, len(vcf_list)))

        reference_counts = defaultdict(int)

        self.log_validation(f"Analyzing reference consensus from {len(sample_files)} sample files...", "INFO")

        for vcf_file in sample_files:
            # Quick reference extraction
            reference, errors = self.extract_vcf_reference(vcf_file)
            if reference and not errors:
                reference_counts[reference] += 1

        if not reference_counts:
            self.log_validation("No valid references found in sample", "WARNING")
            return set(), {}

        # Calculate statistics
        total_samples = len(sample_files)
        dominant_ref = max(reference_counts.keys(), key=lambda x: reference_counts[x])
        dominant_percentage = (reference_counts[dominant_ref] / total_samples) * 100

        self.log_validation(f"Reference analysis: {len(reference_counts)} unique references found", "INFO")
        self.log_validation(f"Most common: {dominant_ref} ({reference_counts[dominant_ref]}/{total_samples} = {dominant_percentage:.1f}%)", "INFO")

        # Determine acceptable references based on consensus
        expected_references = set()

        # If one reference dominates (>70%), expect primarily that reference
        if dominant_percentage > 70:
            expected_references.add(dominant_ref)
            self.log_validation(f"Single dominant reference detected: {dominant_ref}", "INFO")

            # Also accept references that appear in at least 5% of samples (could be legitimate segments/variants)
            for ref, count in reference_counts.items():
                percentage = (count / total_samples) * 100
                if percentage >= 5.0:
                    expected_references.add(ref)
                    if ref != dominant_ref:
                        self.log_validation(f"Accepting secondary reference: {ref} ({percentage:.1f}%)", "INFO")

        # If references are more evenly distributed, accept the top ones
        else:
            # Accept references that appear in at least 10% of samples
            for ref, count in reference_counts.items():
                percentage = (count / total_samples) * 100
                if percentage >= 10.0:
                    expected_references.add(ref)
                    self.log_validation(f"Accepting reference: {ref} ({percentage:.1f}%)", "INFO")

        return expected_references, dict(reference_counts)

    def _attempt_encoding_fix(self, vcf_path: str) -> Tuple[bool, str]:
        """
        Attempt to fix encoding issues in a VCF file by removing non-ASCII characters.

        Reads the file with latin-1 encoding (which tolerates all byte values 0x00-0xFF),
        strips non-ASCII characters from every line, and writes the cleaned content back
        to the original file in-place.

        Returns:
            Tuple of (success: bool, message: str)
        """
        try:
            with open(vcf_path, 'r', encoding='latin-1') as f:
                content = f.read()

            # Count non-ASCII characters before cleaning
            non_ascii_count = sum(1 for c in content if ord(c) >= 128)
            if non_ascii_count == 0:
                return False, "No non-ASCII characters found to remove"

            # Strip non-ASCII chars from every line (preserving line endings)
            cleaned_lines = []
            for line in content.splitlines(keepends=True):
                cleaned_line = line.encode('ascii', errors='ignore').decode('ascii')
                cleaned_lines.append(cleaned_line)

            cleaned_content = ''.join(cleaned_lines)

            with open(vcf_path, 'w', encoding='ascii') as f:
                f.write(cleaned_content)

            return True, f"Removed {non_ascii_count} non-ASCII character(s) from VCF file"

        except Exception as e:
            return False, f"Encoding fix failed: {str(e)}"

    def validate_vcf_list(self, vcf_list: List[str]) -> VCFValidationResults:
        """Validate list of VCF files comprehensively with consensus-based reference checking"""
        results = VCFValidationResults()
        results.total_files = len(vcf_list)

        if not vcf_list:
            results.add_error("", "empty_list", "No VCF files provided", "critical")
            return results

        self.log_validation(f"Validating {len(vcf_list)} VCF files...")

        # Step 1: Analyze reference consensus from a sample of files
        expected_references, consensus_analysis = self._analyze_reference_consensus(vcf_list)

        # Track references found during full validation
        reference_counts = defaultdict(list)

        for vcf_file in vcf_list:
            # Validate file access (no individual logging)
            is_accessible, access_msg = self.validate_file_exists(vcf_file, "VCF")
            if not is_accessible:
                if "does not exist" in access_msg:
                    results.add_error(vcf_file, "unreadable", access_msg, "critical")
                    self.log_validation(f"Unreadable file: {os.path.basename(vcf_file)}", "WARNING")
                elif "not readable" in access_msg:
                    results.add_error(vcf_file, "permission", access_msg, "critical")
                    self.log_validation(f"Permission error: {os.path.basename(vcf_file)}", "WARNING")
                elif "empty" in access_msg:
                    results.add_error(vcf_file, "empty", access_msg, "warning")
                    self.log_validation(f"Empty file: {os.path.basename(vcf_file)}", "WARNING")
                continue

            # Validate VCF format and extract reference (no individual logging)
            is_valid_vcf, vcf_msg, reference = self.validate_vcf_format(vcf_file)
            if not is_valid_vcf:
                if "no data lines" in vcf_msg.lower() or "empty" in vcf_msg.lower():
                    results.add_error(vcf_file, "empty", vcf_msg, "warning")
                    self.log_validation(f"Empty VCF: {os.path.basename(vcf_file)}", "WARNING")
                    continue
                elif "codec can't decode" in vcf_msg or "ordinal not in range" in vcf_msg:
                    # Encoding error — attempt in-place fix
                    original_error_msg = vcf_msg
                    fix_success, fix_msg = self._attempt_encoding_fix(vcf_file)
                    if fix_success:
                        # Re-validate the now-cleaned file
                        is_valid_vcf, vcf_msg, reference = self.validate_vcf_format(vcf_file)
                        if is_valid_vcf:
                            results.add_fixed_file(vcf_file, original_error_msg, fix_msg)
                            self.log_validation(
                                f"Auto-fixed encoding: {os.path.basename(vcf_file)} - {fix_msg}", "INFO"
                            )
                            # Fall through to reference tracking and valid_files append below
                        else:
                            results.add_error(vcf_file, "corrupted", vcf_msg, "critical")
                            self.log_validation(
                                f"Corrupted VCF (fix attempted, still invalid): {os.path.basename(vcf_file)} - {vcf_msg}", "ERROR"
                            )
                            continue
                    else:
                        results.add_error(vcf_file, "corrupted", original_error_msg, "critical")
                        self.log_validation(
                            f"Corrupted VCF: {os.path.basename(vcf_file)} - {original_error_msg}", "ERROR"
                        )
                        continue
                else:
                    results.add_error(vcf_file, "corrupted", vcf_msg, "critical")
                    self.log_validation(f"Corrupted VCF: {os.path.basename(vcf_file)} - {vcf_msg}", "ERROR")
                    continue

            # Track reference usage
            if reference:
                reference_counts[reference].append(vcf_file)

            # File passed all validations (including auto-fixed files)
            results.valid_files.append(vcf_file)

        # Step 2: Check reference consistency using consensus analysis
        if expected_references and len(reference_counts) > 0:
            # Determine the dominant reference from ALL validated files (more reliable than sample)
            dominant_ref_all = max(reference_counts.keys(), key=lambda x: len(reference_counts[x]))

            # Find files whose reference is incompatible with the dominant reference.
            # A reference is accepted if it is either:
            #   1. An exact match to a sample-derived expected reference, OR
            #   2. A segment/chromosome of the same reference genome as the dominant
            #      (detected via long common prefix - handles multi-FASTA references
            #      like influenza where samples may only have SNPs on one segment)
            mismatched_files = []

            for ref, files in reference_counts.items():
                if ref not in expected_references and not self._is_same_reference_genome(ref, dominant_ref_all):
                    mismatched_files.extend(files)

            # Only report mismatches if they represent a small minority
            # This prevents legitimate multi-reference datasets from being flagged
            total_valid_files = len(results.valid_files)
            mismatch_percentage = (len(mismatched_files) / total_valid_files * 100) if total_valid_files > 0 else 0

            if mismatched_files and mismatch_percentage < 30:  # Only flag if <30% of files are "mismatched"
                self.log_validation(f"Reference mismatch detected - {len(mismatched_files)} files use unexpected references", "ERROR")

                # Determine the most common expected reference for error reporting
                if consensus_analysis:
                    most_common_expected = max(consensus_analysis.keys(), key=lambda x: consensus_analysis[x])
                else:
                    most_common_expected = list(expected_references)[0] if expected_references else "unknown"

                for vcf_file in mismatched_files:
                    # Find what reference this file actually uses
                    file_reference = None
                    for ref, files in reference_counts.items():
                        if vcf_file in files:
                            file_reference = ref
                            break

                    if file_reference:
                        mismatch_msg = f"Expected reference: {most_common_expected}, Found: {file_reference}"
                        results.add_error(vcf_file, "reference_mismatch", mismatch_msg, "critical")
                        self.log_validation(f"Reference mismatch: {os.path.basename(vcf_file)} uses {file_reference}", "ERROR")
                        # Remove from valid files
                        if vcf_file in results.valid_files:
                            results.valid_files.remove(vcf_file)

            elif mismatch_percentage >= 30:
                self.log_validation(f"Multiple references detected ({len(reference_counts)} unique), but distribution suggests legitimate multi-reference dataset", "INFO")
                for ref, files in reference_counts.items():
                    percentage = (len(files) / total_valid_files * 100)
                    self.log_validation(f"Reference {ref}: {len(files)} files ({percentage:.1f}%)", "INFO")

        self.log_validation(f"Validation complete: {results.total_valid}/{results.total_files} files valid")
        return results

    def print_vcf_validation_summary(self, results: VCFValidationResults):
        """Print comprehensive VCF validation results to terminal"""
        print("\n" + "="*80)
        print("📊 VCF FILE VALIDATION SUMMARY")
        print("="*80)
        clean_count = results.total_valid - len(results.fixed_files)
        print(f"Total files processed: {results.total_files}")
        print(f"Valid files: {results.total_valid}")
        if results.fixed_files:
            print(f"  - Clean (no issues): {clean_count}")
            print(f"  - Auto-fixed (encoding corrected): {len(results.fixed_files)}")
        print(f"Invalid files: {results.total_invalid}")

        if results.total_invalid == 0 and not results.fixed_files:
            print("✅ All VCF files passed validation!")
            return

        # Show auto-fixed files
        if results.fixed_files:
            print(f"\n🔧 AUTO-FIXED: Encoding Issues ({len(results.fixed_files)} files) - included in analysis")
            print("-" * 50)
            for entry in results.fixed_files:
                print(f"  🔧 {os.path.basename(entry['file'])}")
                print(f"     Fix: {entry['fix_applied']}")

        if results.total_invalid == 0:
            print("="*80)
            return

        # Show remaining validation issues
        if results.reference_mismatches:
            print(f"\n🔥 CRITICAL: Reference Mismatch ({len(results.reference_mismatches)} files)")
            print("-" * 50)
            for mismatch in results.reference_mismatches:
                print(f"  ❌ {os.path.basename(mismatch['file'])}")
                print(f"     {mismatch['message']}")

        if results.corrupted_files:
            print(f"\n💥 CRITICAL: Corrupted Files ({len(results.corrupted_files)} files)")
            print("-" * 50)
            for file_path in results.corrupted_files:
                print(f"  ❌ {os.path.basename(file_path)}")

        if results.empty_files:
            print(f"\n⚠️  WARNING: Empty Files ({len(results.empty_files)} files)")
            print("-" * 50)
            for file_path in results.empty_files:
                print(f"  ⚠️  {os.path.basename(file_path)}")

        if results.permission_errors:
            print(f"\n🔒 ERROR: Permission Issues ({len(results.permission_errors)} files)")
            print("-" * 50)
            for file_path in results.permission_errors:
                print(f"  🔒 {os.path.basename(file_path)}")

        if results.unreadable_files:
            print(f"\n❓ ERROR: Unreadable Files ({len(results.unreadable_files)} files)")
            print("-" * 50)
            for file_path in results.unreadable_files:
                print(f"  ❓ {os.path.basename(file_path)}")

        print("="*80)

    def write_validation_log(self, log_file: str, results: VCFValidationResults):
        """Write detailed validation log to file"""
        try:
            with open(log_file, 'w') as f:
                f.write("VCF VALIDATION LOG\n")
                f.write("==================\n\n")
                clean_count = results.total_valid - len(results.fixed_files)
                f.write(f"Total files processed: {results.total_files}\n")
                f.write(f"Valid files: {results.total_valid}\n")
                if results.fixed_files:
                    f.write(f"  - Clean (no issues): {clean_count}\n")
                    f.write(f"  - Auto-fixed (encoding corrected): {len(results.fixed_files)}\n")
                f.write(f"Invalid files: {results.total_invalid}\n\n")

                # Write validation messages
                f.write("VALIDATION MESSAGES:\n")
                f.write("-" * 40 + "\n")
                for log_msg in self.validation_log:
                    f.write(f"{log_msg}\n")
                f.write("\n")

                # Write auto-fixed files section
                if results.fixed_files:
                    f.write(f"AUTO-FIXED FILES ({len(results.fixed_files)}) - included in analysis:\n")
                    f.write("-" * 40 + "\n")
                    for entry in results.fixed_files:
                        f.write(f"File: {entry['file']}\n")
                        f.write(f"Fix applied: {entry['fix_applied']}\n")
                        f.write(f"Original error: {entry['original_error']}\n")
                        f.write("-" * 40 + "\n")
                    f.write("\n")

                # Write detailed error information
                if results.errors:
                    f.write("DETAILED ERRORS:\n")
                    f.write("-" * 40 + "\n")
                    for error in results.errors:
                        f.write(f"File: {error.file_path}\n")
                        f.write(f"Type: {error.error_type}\n")
                        f.write(f"Severity: {error.severity}\n")
                        f.write(f"Message: {error.error_message}\n")
                        f.write("-" * 40 + "\n")

                # Write valid files list
                if results.valid_files:
                    fixed_paths = {e['file'] for e in results.fixed_files}
                    f.write(f"\nVALID FILES ({len(results.valid_files)}):\n")
                    f.write("-" * 40 + "\n")
                    for vcf_file in results.valid_files:
                        suffix = " [auto-fixed]" if vcf_file in fixed_paths else ""
                        f.write(f"{vcf_file}{suffix}\n")

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
        vcf_files = glob.glob(os.path.join(test_path, "*.vcf"))
        results = validator.validate_vcf_list(vcf_files)
        validator.print_vcf_validation_summary(results)