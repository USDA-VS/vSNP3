#!/usr/bin/env python

from vsnp3_version import __version__

import os
import sys

# Answer --version before the heavy imports below (pandas, Biopython, dask),
# which take seconds on a cold network file system: the vSNP GUI probes the
# version with a 10 s timeout and recorded an empty version when they were
# slow.  Same text as the argparse version action.
if __name__ == '__main__' and any(a in ('-v', '--version') for a in sys.argv[1:]):
    print(f'{os.path.basename(__file__)}: version {__version__}')
    raise SystemExit(0)
import io
import html
import shutil
import subprocess
import re
import json
import pickle
import locale
import argparse
import textwrap
import importlib.metadata
import pandas as pd
import numpy as np
import zipfile
import glob
from datetime import datetime

from collections import defaultdict, Counter
from concurrent import futures
import multiprocessing

# Move set_start_method inside if __name__ == "__main__" to avoid issues with Python 3.12
# multiprocessing.set_start_method('spawn', True)

import warnings
# Targeted, not blanket.  The previous filterwarnings('ignore') hid real defects:
# a SettingWithCopyWarning on a write-to-a-slice in make_groupings, and pandas
# deprecations for positional Series access that is removed in pandas 3.0.  Only
# the known-noisy messages are suppressed so anything new is visible.
for _msg in (r'.*invalid value encountered.*',
             r'.*divide by zero encountered.*',
             r'.*DataFrame is highly fragmented.*',
             r'.*Passing a BlockManager.*'):
    warnings.filterwarnings('ignore', message=_msg)
warnings.filterwarnings('ignore', category=DeprecationWarning, module='openpyxl')

from vsnp3_version import SYNTHESIZED_QUAL
from vsnp3_version import QUAL_THRESHOLD, N_THRESHOLD, MQ_THRESHOLD, MIXED_AC1_RATIO
from vsnp3_file_setup import Setup, PhaseTimer, format_elapsed, worker_count
from vsnp3_file_setup import redact_paths
from vsnp3_group_on_defining_snps import Group, DefiningSnpsError
from vsnp3_fasta_to_snps_table import bootstrap_replicates
from vsnp3_reference_options import Ref_Options
from vsnp3_remove_from_analysis import Remove_From_Analysis
from vsnp3_input_validator import InputValidator, VCFValidationResults

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

global_date_stamp=None
global_working_dir='.'


# ================================================================ VCF reading
#
# One read per file.  The 3.36 path read every VCF three or four times: a
# header check, a line-by-line normalizing rewrite to vcf_normalized/ (a dozen
# regex substitutions on every line, most of them no-ops), a pandas read of
# that copy, a second rewrite when a header note was needed, plus a copy to
# vcf_starting_files/ and a copy into the zip.  Then the INFO column was split
# into a dict three times per record to pull out AC, DP and MQ.
#
# Now the file is read once, in binary so the process locale cannot make it
# undecodable, normalized in memory with the regex work done only on lines that
# need it, parsed by pandas directly from memory, and AC/DP/MQ are extracted
# with one vectorized regex each.  The normalized copy is written only for
# -f/--fix_vcfs and --debug.

_NON_ASCII = bytes(range(128, 256))
_FIXED_COLUMNS = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
# An empty line inside the body: the point between two consecutive newlines.
_BLANK_LINE = re.compile(r'(?<=\n)(?=\n)')


def _normalize_header_line(line):
    '''The 3.36 vcf_fix treatment of one ## header line.'''
    line = line.replace(';MQM=', ';MQ=').replace('ID=MQM,', 'ID=MQ,')
    if '"' in line:
        line = line.replace('"AC=', 'AC=')
        line = line.replace('""', '"').replace('""', '"').replace('""', '"')
        line = re.sub('"$', '', line)
        line = line.replace('GQ:PL\t"', 'GQ:PL\t')
        line = re.sub('^"', '', line)
        if line.endswith('"'):
            line = line[:-1]
    return line.split('\t')[0]


def normalize_vcf_text(text, basename, assume_gt_only_quality=False, first_sample_only=False):
    '''
    Split a VCF into (header_lines, column_names, body_text, notes).

    The same repairs vcf_fix made in 3.36 -- freebayes MQM renamed to the
    standard MQ, stray quoting stripped, runs of spaces turned into tabs --
    but applied in memory and only where a line actually contains the thing
    being repaired.  notes carries what was done so it can be recorded.
    '''
    notes = {'synthesized_qual': False, 'multi_sample': 0, 'platform': None,
             'longread_model': None}
    header_lines = []
    column_line = None
    pos = 0
    n = len(text)
    while pos < n:
        nl = text.find('\n', pos)
        end = n if nl < 0 else nl
        line = text[pos:end].rstrip()
        pos = end + 1
        if not line:
            continue
        if line.startswith('##'):
            # A nanopore VCF names its platform and the parameters of the
            # evidence QUAL it carries (vsnp3_longread.py); kept so the summary
            # can say which samples are long-read and how their QUAL was made.
            if line.startswith('##vsnp3_platform='):
                notes['platform'] = line.split('=', 1)[1].strip()
            elif line.startswith('##vsnp3_longread_model='):
                notes['longread_model'] = line.split('=', 1)[1].strip()
            header_lines.append(_normalize_header_line(line))
            continue
        if line.startswith('#'):
            column_line = line
            break
        raise ValueError('a data record appears before the #CHROM column header line')
    if column_line is None:
        raise ValueError('no #CHROM column header line')
    body_text = text[pos:]

    column_line = re.sub(' +', '\t', column_line.replace('"', ''))
    columns = column_line.split('\t')
    if len(columns) < 8:
        raise ValueError(f'#CHROM line has {len(columns)} column(s); a VCF needs at least 8')
    if len(columns) > 10:
        notes['multi_sample'] = len(columns) - 9
        if not first_sample_only:
            raise ValueError(
                f'multi-sample VCF ({len(columns) - 9} samples: '
                f'{", ".join(columns[9:])}). vSNP3 expects one sample per '
                f'file; split with `bcftools view -s <sample>`, or pass '
                f'--first_sample_only to analyse {columns[9]} alone.')

    needs_line_pass = ('"' in body_text or ' ' in body_text or '\r' in body_text
                       or assume_gt_only_quality)
    if needs_line_pass:
        out = []
        for line in body_text.splitlines():
            line = line.rstrip()
            if not line:
                continue
            if '"' in line:
                line = line.replace('"', '')
            if ' ' in line:
                line = re.sub(' +', '\t', line)
            if assume_gt_only_quality:
                # Records carrying only a GT have no quality data, so every
                # threshold downstream rejects them.  On request substitute a
                # QUAL and leave everything else alone.
                fields = line.split('\t')
                if len(fields) > 9 and fields[8] == 'GT':
                    fields[5] = str(SYNTHESIZED_QUAL)
                    notes['synthesized_qual'] = True
                    line = '\t'.join(fields)
            out.append(line)
        body_text = '\n'.join(out)
    body_text = body_text.replace(';MQM=', ';MQ=')
    return header_lines, columns, body_text, notes


def write_normalized_copy(path, header_lines, columns, body_text, notes):
    '''
    The vcf_normalized/ file 3.36 wrote for every input, now only on request.
    Decisions taken during normalization are recorded in the header so they
    are not invisible.
    '''
    with open(path, 'w') as out:
        for line in header_lines:
            print(line, file=out)
        if notes.get('synthesized_qual'):
            print(f'##vsnp3_synthesized_quality=QUAL={SYNTHESIZED_QUAL} substituted for '
                  f'records whose FORMAT was GT only (--assume_gt_only_quality)', file=out)
        if notes.get('multi_sample'):
            print('##vsnp3_first_sample_only=other samples in this VCF were dropped', file=out)
        print('\t'.join(columns[:10]), file=out)
        for line in body_text.splitlines():
            if line:
                print('\t'.join(line.split('\t')[:10]), file=out)


def dataframe_from_vcf_body(columns, body_text):
    '''
    Parse the data records into the frame the rest of step 2 works on.

    Returns (df, stats).  Columns: CHROM, POS, REF, ALT, QUAL, AC, DP, MQ,
    CIGAR (multi-base records only), abs_pos.  The sample column is not read at all: nothing downstream uses
    it, and in a zero-coverage VCF it was a per-record string carried through
    every later stage for no purpose.
    '''
    if body_text.strip():
        # Blank lines are not records.  pandas skips them, so they are left
        # out of the count too: counted, a VCF ending in a blank line reported
        # one record that "could not be read and was skipped".
        n_lines = (body_text.count('\n') + (0 if body_text.endswith('\n') else 1)
                   - len(_BLANK_LINE.findall(body_text)) - body_text.startswith('\n'))
        # Only the eight fixed columns are named and read, so a record with
        # more fields than the column header is read like any other and the
        # extra fields are ignored, as FORMAT and the sample column are.  The
        # header's FORMAT and sample columns used to be named as well; pandas
        # then compared the first data line with that count and rejected the
        # whole file when that one line had an extra field, while the same
        # line further down was read normally.  3.36 read every such record
        # too, truncated to ten fields.  A garbage line still has no numeric
        # POS and is dropped and counted as malformed below.
        df = pd.read_csv(io.StringIO(body_text), sep='\t', header=None,
                         names=_FIXED_COLUMNS, usecols=range(8), dtype=str,
                         na_filter=False, on_bad_lines='skip', skip_blank_lines=True,
                         engine='c')
    else:
        n_lines = 0
        df = pd.DataFrame({c: pd.Series(dtype=object) for c in _FIXED_COLUMNS})
    skipped = max(0, n_lines - len(df))

    df['POS'] = pd.to_numeric(df['POS'], errors='coerce')
    malformed = int(df['POS'].isna().sum())
    if malformed:
        df = df[df['POS'].notna()]
    df['POS'] = df['POS'].astype(np.int64)
    # QUAL and MQ stay nullable.  These were .fillna(0), which turned an ABSENT
    # quality value into the assertion "quality is zero".  For QUAL that means
    # the position falls below n_threshold and is rewritten to the REFERENCE
    # base - so a missing measurement became a positive claim that the sample
    # matches the reference.  A missing value is a no-call; Int64 carries that
    # distinction through to make_groupings, which tests for it.  np.trunc keeps
    # the truncation .astype(int) performed on fractional QUAL such as 4523.5.
    df['QUAL'] = np.trunc(pd.to_numeric(df['QUAL'], errors='coerce')).astype('Int64')

    info = df['INFO']
    df['AC'] = pd.to_numeric(info.str.extract(r'(?:^|;)AC=([^;]*)', expand=False),
                             errors='coerce').fillna(0).astype(np.int64)
    df['DP'] = pd.to_numeric(info.str.extract(r'(?:^|;)DP=([^;]*)', expand=False),
                             errors='coerce').fillna(0).astype(np.int64)
    # An absent MQ is unknown, not zero.  As zero it silently failed the
    # MQ >= 56 filter, so those positions were dropped with no indication.
    df['MQ'] = np.trunc(pd.to_numeric(info.str.extract(r'(?:^|;)MQ=([^;]*)', expand=False),
                                      errors='coerce')).astype('Int64')
    # For a record whose REF spans more than one base (MNP, complex or indel),
    # the CIGAR says which reference base each ALT base aligns to.  Kept for those
    # records only: make_groupings reads a SNP freebayes reported inside one of
    # them from it.  NaN for every other record.
    multi = df['REF'].str.len() > 1
    df['CIGAR'] = info.where(multi).str.extract(r'(?:^|;)CIGAR=([^;]*)', expand=False)
    # Both alleles one base.  Group._passing and _confident_mixed test this on
    # every sample; computed once here instead of two str.len() passes there.
    df['single_base'] = (df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)

    called = df['QUAL'].notna()
    if called.any() and df.loc[called, 'MQ'].isna().all():
        # A sample whose records carry QUAL but never MQ cannot have the mapping
        # quality filter applied.  A file with no called records at all (only
        # zero-coverage records, or nothing) is a different thing: that sample
        # is identical to the reference wherever it has coverage, and it is
        # analysed as such rather than excluded.
        raise ValueError(
            'no MQ field in INFO on any called record, so the mapping quality '
            'filter cannot be applied. Sample excluded rather than contributing '
            'zero positions silently. Re-call with a caller that emits MQ (or MQM).')

    df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
    stats = {
        'records': int(len(df)),
        'zero_coverage': int((~called).sum()),
        'called': int(called.sum()),
        'ac1': int((df['AC'] == 1).sum()),
        'ac2': int((df['AC'] == 2).sum()),
        'malformed': malformed,
        'skipped_lines': skipped,
    }
    df = df.drop(columns=['INFO', 'ID', 'FILTER'])
    df = df.drop_duplicates(subset=['abs_pos'])
    return df, stats


def parse_vcf(vcf, assume_gt_only_quality=False, first_sample_only=False, normalized_dir=None):
    '''
    Read, normalize and parse one VCF.  Returns (vcf, df, stats, error).

    Never raises: a file that cannot be parsed comes back with df None and
    the reason in error, so the caller can account for it.  The input file is
    opened read-only and is never modified.
    '''
    base = os.path.basename(vcf)
    stats = {'non_ascii_bytes': 0, 'multi_sample': 0, 'synthesized_qual': False,
             'platform': None, 'longread_model': None}
    try:
        with open(vcf, 'rb') as f:
            data = f.read()
        if not data.isascii():
            # Usually an accented letter or an en dash in a ##reference path.
            # Under the C locale every module forces, a text-mode open() decoded
            # as ASCII and raised on the first such byte, which is how a good
            # file came to be called corrupt and rewritten in place.  The bytes
            # are dropped from the in-memory copy, counted, and the file is left
            # exactly as it was.
            kept = data.translate(None, _NON_ASCII)
            stats['non_ascii_bytes'] = len(data) - len(kept)
            data = kept
        text = data.decode('ascii')
        del data
        header_lines, columns, body_text, notes = normalize_vcf_text(
            text, base, assume_gt_only_quality, first_sample_only)
        del text
        stats['multi_sample'] = notes['multi_sample']
        stats['synthesized_qual'] = notes['synthesized_qual']
        stats['platform'] = notes.get('platform')
        stats['longread_model'] = notes.get('longread_model')
        if normalized_dir:
            write_normalized_copy(os.path.join(normalized_dir, base),
                                  header_lines, columns, body_text, notes)
        df, parse_stats = dataframe_from_vcf_body(columns, body_text)
        stats.update(parse_stats)
        return vcf, df, stats, None
    except Exception as e:                                      # noqa: BLE001
        return vcf, None, stats, f'{type(e).__name__}: {e}'


def _parse_task(task):
    '''Pool entry point: a plain tuple in, nothing pickled but the result.'''
    return parse_vcf(*task)


class VCF_to_DF():
    '''
    Validate, read and parse the VCF files, keeping an account of every one.
    '''

    def __init__(self, vcf_list=None, debug=False, assume_gt_only_quality=False,
                 first_sample_only=False, write_normalized=False, timer=None):
        self.startTime = datetime.now()
        self.timer = timer or PhaseTimer()
        self.vcf_original_count = len(vcf_list)
        self.assume_gt_only_quality = assume_gt_only_quality
        self.first_sample_only = first_sample_only
        self.chrom = None
        # basename -> (bucket, reason) for every file that will not be analysed
        self.excluded = {}
        # basename -> parse statistics for every file that was parsed
        self.sample_stats = {}
        dataframes = {}

        self.normalized_dir = os.path.join(os.getcwd(), 'vcf_normalized')
        normalized_dir = None
        if write_normalized or debug:
            os.makedirs(self.normalized_dir, exist_ok=True)
            normalized_dir = self.normalized_dir

        validator = InputValidator(debug=debug)
        print(f"\n🔍 VALIDATING {self.vcf_original_count} VCF FILES...")
        print("="*60)
        validation_results = validator.validate_vcf_list(vcf_list)
        self.validation_results = validation_results  # Store for HTML reporting
        validator.print_vcf_validation_summary(validation_results)

        log_file = f'vcf_validation_log-{global_date_stamp}.txt'
        validator.write_validation_log(log_file, validation_results)
        print(f"\n📋 Detailed validation log written to: {log_file}")

        for vcf_file in validation_results.excluded_files:
            self.excluded[os.path.basename(vcf_file)] = (
                'failed validation', validation_results.reason(vcf_file))

        if validation_results.reference_mismatches:
            print("\n" + "="*80)
            print("🔥 VCF REFERENCE MISMATCH")
            print("="*80)
            print("\nThese VCF files were called against a different reference and are "
                  "excluded from the analysis:")
            print("-" * 60)
            for mismatch in validation_results.reference_mismatches:
                print(f"❌ {os.path.basename(mismatch['file'])}")
                print(f"   {mismatch['message']}")
            print("\n💡 All VCF files in one step 2 run must use the same reference "
                  "genome.  Re-run step 1 for these samples with the reference the "
                  "rest of the set uses, or move them to their own directory.")
            print("📋 Full details in validation log: " + log_file)
            print("="*80 + "\n")
            if validation_results.reference_majority_lost:
                # Not a few strays: the set itself is mixed.  Continuing would
                # analyse whichever reference happened to be most common and
                # call the rest excluded, which is not what anyone asked for.
                print("### No reference is shared by a majority of the files, so this "
                      "directory holds a mixed set.  Nothing was analysed.")
                for note in validation_results.notes:
                    print(f"### {note}")
                sys.exit(1)
        self.timer.mark('Validate VCF headers')

        valid_vcf_list = validation_results.valid_files
        print(f"\n✅ PROCEEDING WITH {len(valid_vcf_list)} VALID VCF FILES")
        if self.excluded:
            print(f"⚠️  EXCLUDED {len(self.excluded)} FILE(S), see the summary report")
        print("="*60 + "\n")

        tasks = [(vcf, assume_gt_only_quality, first_sample_only, normalized_dir)
                 for vcf in valid_vcf_list]
        # Capped by VSNP3_MAX_CPUS (see worker_count).  The vSNP GUI's install
        # patch looks for that name here and leaves this file alone when found.
        workers = worker_count(len(tasks))
        every = max(1, len(tasks) // 10)

        def collect(index, result):
            vcf, df, stats, error = result
            base = os.path.basename(vcf)
            if stats.get('multi_sample'):
                print(f'  {base}: multi-sample VCF, keeping only the first sample', flush=True)
            if error is not None:
                # Previously the input file was deleted here.  A parse failure is a
                # reason to exclude a sample and say why, not to destroy the evidence.
                self.excluded[base] = ('could not be parsed', error)
            else:
                dataframes[base] = df
                self.sample_stats[base] = stats
                if self.chrom is None and len(df):
                    self.chrom = df['CHROM'].iloc[0]
            if (index + 1) % every == 0 or index + 1 == len(tasks):
                print(f'  parsed {index + 1}/{len(tasks)} VCF files '
                      f'({(index + 1) * 100 // len(tasks)}%)', flush=True)

        if debug or len(tasks) < 4:
            print(f'Processing {len(tasks)} VCF files in this process')
            for index, task in enumerate(tasks):
                collect(index, _parse_task(task))
        else:
            print(f'Processing: Pool processing with {workers} cpus...')
            chunksize = max(1, min(8, len(tasks) // (workers * 4)))
            with futures.ProcessPoolExecutor(max_workers=workers) as pool:
                for index, result in enumerate(pool.map(_parse_task, tasks, chunksize=chunksize)):
                    collect(index, result)

        if self.chrom is None and validation_results.reference_contigs:
            self.chrom = validation_results.reference_contigs[0]

        self.dataframes = dataframes
        non_ascii = sum(1 for s in self.sample_stats.values() if s['non_ascii_bytes'])
        if non_ascii:
            print(f'  {non_ascii} file(s) contained non-ASCII bytes, ignored on read; '
                  f'the files were not changed.')
        longread = sorted(b for b, s in self.sample_stats.items()
                          if s.get('platform') == 'nanopore')
        if longread:
            print(f'  {len(longread)} long-read (nanopore) sample(s), whose QUAL is the vSNP3 '
                  f'evidence QUAL and is read by the same thresholds: '
                  f'{", ".join(longread[:8])}{", ..." if len(longread) > 8 else ""}')
        self.timer.mark('Read and parse VCFs')
        print(f'\n\nDictionary of dataframes to memory runtime: {datetime.now() - self.startTime}\n')

        if debug:
            # For vsnp3_group_on_defining_snps.py -p.  Written only on request:
            # 3.36 pickled every dataframe to disk on every run and deleted the
            # file a moment later.
            with open('dictionary_of_dataframes.pickle', 'wb') as handle:
                pickle.dump(dataframes, handle, protocol=pickle.HIGHEST_PROTOCOL)

    @property
    def vcf_bad_list(self):
        '''Excluded files as "file  reason" strings, for callers that print them.'''
        return [f'{base}  {reason}' for base, (bucket, reason) in self.excluded.items()]


# Density filtering defaults (--density_threshold / --density_window), written
# once: the argparse constants and the summary text used to carry them separately.
DENSITY_THRESHOLD_DEFAULT = 3
DENSITY_WINDOW_DEFAULT = 20


# ============================================================ sample accounting

class SampleAccounting:
    '''
    Where every input file went.

    Each basename ends in exactly one of: analysed, or one of the not-analysed
    buckets below, each with a reason.  A sample that is analysed but has no
    SNP passing the thresholds is still analysed; it is listed separately so
    the reader can see why its row carries no calls of its own.
    '''
    BUCKETS = ('failed validation', 'could not be parsed',
               'removed by remove_from_analysis', 'duplicate sample name')

    def __init__(self, vcf_files):
        self.files = [os.path.basename(f) for f in vcf_files]
        self.excluded = {}       # basename -> (bucket, reason)
        self.analysed = {}       # sample name -> basename
        self.no_snps = {}        # sample name -> detail dict from Group
        self.remove_unmatched = []   # remove-list entries that named no file

    def exclude(self, basename, bucket, reason):
        self.excluded.setdefault(basename, (bucket, reason))

    def by_bucket(self):
        out = {bucket: [] for bucket in self.BUCKETS}
        for base, (bucket, reason) in self.excluded.items():
            out.setdefault(bucket, []).append((base, reason))
        return out

    @property
    def unaccounted(self):
        seen = set(self.excluded) | set(self.analysed.values())
        return [f for f in self.files if f not in seen]


# ================================================================ HTML summary

def _conda_versions():
    '''
    {package: version} for the active conda environment, from ONE call.

    3.36 ran `conda list <program>` twice for each of 13 programs -- 26
    launches of the conda interpreter -- after the run time had already been
    computed, which is one reason the reported run time disagreed with the
    file's modification time.
    '''
    try:
        out = subprocess.check_output(['conda', 'list', '--json'],
                                      stderr=subprocess.DEVNULL,
                                      universal_newlines=True, timeout=180)
        return {e['name'].lower(): e.get('version', 'nd')
                for e in json.loads(out) if e.get('name')}
    except (subprocess.SubprocessError, OSError, ValueError):
        return {}


class HTML_Summary():

    CSS = ('table { font-family: arial, sans-serif; border-collapse: collapse; } '
           'td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; '
           'vertical-align: top; } '
           '.box { font-size:11px; border:1px solid #8a5a00; background:#fdf3e0; color:#8a5a00; '
           'padding:6px; margin:6px 0; } '
           '.bad { color: #a52121; } .ok { color: #1a7f4b; } .warn { color: #8a5a00; } '
           '.small { font-size: 10px; margin-left: 20px; } '
           '.mixed { color: #a52121; background: #fdf3e0; font-weight: bold; }')

    def __init__(self, timer=None, vcf_to_df=None, reference=None, groupings_dict=None,
                 raxml_version=None, all_vcf_boolen=None, args=None,
                 validation_results=None, name_collisions=None, metadata_ambiguous=None,
                 accounting=None, metadata_problem=None, group=None, reference_note=None):
        self.out = []
        self.timer = timer
        self.vcf_to_df = vcf_to_df
        self.reference = reference
        self.groupings_dict = groupings_dict
        self.raxml_version = raxml_version
        self.all_vcf = all_vcf_boolen
        self.args = args
        self.validation_results = validation_results
        self.name_collisions = name_collisions or {}
        self.metadata_ambiguous = metadata_ambiguous or {}
        self.accounting = accounting
        self.metadata_problem = metadata_problem
        self.reference_note = reference_note
        self.group = group
        self.path = f'{global_working_dir}/vSNP_step2_summary-{global_date_stamp}.html'

        self.header()
        self.inputs()
        self.warnings()
        self.thresholds()
        self.long_read_samples()
        self.tree_failures()
        self.group_failures()
        self.sample_accounting()
        self.groupings()
        self.not_analysed()
        self.system_information()
        self.program_versions()
        # Timing is written last, after everything above has been gathered,
        # so the elapsed time covers the whole run.  3.36 computed it before
        # 26 conda calls and the version checks, so the report could say 5.5
        # hours while its own modification time said 13.
        self.timing()
        self.p("</body>\n</html>")
        with open(self.path, 'w', encoding='utf-8') as htmlfile:
            htmlfile.write('\n'.join(self.out))
        # The same accounting, machine-readable, for the vSNP GUI and other
        # callers that otherwise scrape this HTML.
        self.json_path = f'{global_working_dir}/vSNP_step2_sample_accounting-{global_date_stamp}.json'
        try:
            self.write_json()
        except Exception as e:                                  # noqa: BLE001
            print(f'### The accounting JSON could not be written: {type(e).__name__}: {e}')

    # ------------------------------------------------------------ helpers

    def p(self, text):
        self.out.append(text)

    @staticmethod
    def esc(value):
        return html.escape(redact_paths(str(value)))

    # ------------------------------------------------------------ sections

    def header(self):
        self.p('<html>\n<head><meta charset="UTF-8"><style>' + self.CSS + '</style></head>\n'
               '<body style="font-size:12px;">')
        self.p(f"<h2>Script ran using <u>{self.esc(self.reference)} </u> variables:</h2>")

    def inputs(self):
        args = self.args
        self.p('<div style="font-size:11px; font-weight:normal;">')
        if args.metadata:
            self.p(f"Metadata:  {self.esc(args.metadata)}<br>")
        else:
            self.p("No metadata for describing samples in trees and tables<br>")
        if args.defining_snps:
            self.p(f"Defining SNPs:  {self.esc(args.defining_snps)}<br>")
        else:
            self.p("No defining SNPs files for grouping and filtering<br>")
        if args.gbk:
            for each in args.gbk:
                self.p(f"gbk:  {self.esc(each)}<br>")
        else:
            self.p("No gbk for annotation<br>")
        if args.remove_by_name:
            self.p(f"Remove from analysis:  {self.esc(args.remove_by_name)}<br>")
        self.p('</div><br>')

    def warnings(self):
        grp = self.group
        if self.metadata_problem:
            self.p('<div class="box">')
            self.p('<b>Metadata not used</b><br>')
            self.p(f'{self.esc(self.metadata_problem)}<br>')
            self.p('The run continued with the VCF file names as sample names.')
            self.p('</div>')
        if self.reference_note:
            self.p('<div class="box">')
            self.p('<b>Reference chosen from the FASTA headers</b><br>')
            self.p(f'{self.esc(self.reference_note)}')
            self.p('</div>')
        blank_rows = getattr(grp, 'metadata_blank_rows', 0)
        if blank_rows:
            self.p('<div class="box">')
            self.p(f'<b>Metadata rows ignored: {blank_rows}</b><br>')
            self.p('These rows of the metadata worksheet have a blank cell in column A or B. '
                   'A blank cell is not a name, so the row was skipped and any VCF it was '
                   'meant to rename keeps its file name.')
            self.p('</div>')
        renames = getattr(grp, 'fallback_renames', None) or {}
        if renames:
            self.p('<div class="box">')
            self.p(f'<b>Sample names changed: {len(renames)}</b><br>')
            self.p('These VCF file names contain characters a Newick tree cannot carry '
                   '(space, quote, bracket, colon, comma or semicolon), which failed RAxML for '
                   'the whole group. The sample is named with "_" in their place:<br>')
            for base, (before, after) in sorted(renames.items()):
                self.p(f'&nbsp;&nbsp;<b>{self.esc(base)}</b>: {self.esc(before)} &rarr; '
                       f'{self.esc(after)}<br>')
            self.p('</div>')
        notes = getattr(grp, 'worksheet_notes', None) or []
        if notes:
            self.p('<div class="box">')
            self.p('<b>Defining SNP worksheet</b><br>')
            for note in notes:
                self.p(f'{self.esc(note)}<br>')
            self.p('</div>')
        # Remove-list entries that matched no VCF are not shown here: the list is
        # usually the reference type's, naming samples of other runs too, so the
        # box was mostly noise.  The terminal and the accounting JSON keep them.
        rewrite = getattr(grp, 'rewrite_errors', None) or {}
        if rewrite:
            self.p('<div class="box">')
            self.p(f'<b>Call rewrite failed for {len(rewrite)} group/sample pair(s)</b><br>')
            self.p('In these samples the positions below the N threshold may keep their ALT '
                   'instead of the reference base. Review them before relying on the table:<br>')
            for (group_name, sample), reason in sorted(rewrite.items()):
                self.p(f'&nbsp;&nbsp;<b>{self.esc(group_name)}</b> / {self.esc(sample)}: '
                       f'{self.esc(reason)}<br>')
            self.p('</div>')

        # Duplicate sample names.  The run completed on one file per name, so the
        # reader has to be told which file that was: the tables and tree carry the
        # name, and nothing else in the report would show that a second VCF claimed
        # it.  Placed high in the summary rather than in a footnote, because it
        # means a sample the user submitted is absent from the analysis.
        if self.name_collisions:
            self.p('<div class="box">')
            self.p(f'<b>Duplicate sample names: {len(self.name_collisions)}</b><br>')
            self.p('More than one VCF resolved to the same sample name, so only one '
                   'was analysed under each. The reason is given per name below; '
                   'correct it and rerun to include the files left out.<br>')
            for name, info in sorted(self.name_collisions.items()):
                dropped = ', '.join(self.esc(d) for d in info['dropped'])
                self.p(f'&nbsp;&nbsp;<b>{self.esc(name)}</b>: analysed '
                       f'{self.esc(info["used"])} &mdash; not analysed: {dropped}<br>')
                self.p(f'&nbsp;&nbsp;&nbsp;&nbsp;<i>because '
                       f'{self.esc(info.get("cause", ""))}</i><br>')
            self.p('</div>')

        # A file matching several worksheet rows keeps one name, so no count
        # changes and nothing else in the report would show it.  Separate box from
        # the one above: a different mistake, needing a different correction.
        if self.metadata_ambiguous:
            self.p('<div class="box">')
            self.p(f'<b>Ambiguous metadata rows: {len(self.metadata_ambiguous)}</b><br>')
            self.p('These VCF files match more than one row of the metadata '
                   'worksheet, and those rows give different names. The first '
                   'matching row was used; remove the duplicate rows to make the '
                   'choice explicit.<br>')
            for base, info in sorted(self.metadata_ambiguous.items()):
                targets = ', '.join(self.esc(t) for t in info['targets'])
                self.p(f'&nbsp;&nbsp;<b>{self.esc(base)}</b>: '
                       f'"{self.esc(info["key"])}" appears more than once, giving {targets} '
                       f'&mdash; used {self.esc(info["targets"][0])}<br>')
            self.p('</div>')

    def thresholds(self):
        args = self.args
        self.p(f'<span style="font-size:11px; font-weight:bold;">SNP calling thresholds:  '
               f'REF: QUAL <u>&lt;{args.n_threshold}</u>, N: QUAL <u>{args.n_threshold}-{args.qual_threshold}</u>, '
               f'ALT: QUAL <u>&gt;{args.qual_threshold}</u>, Ambigious: <u>AC=1</u>, '
               f'MQ: <u>&ge;{args.mq_threshold}</u></span>')
        if args.density_threshold is not None or args.density_window is not None:
            threshold = args.density_threshold if args.density_threshold is not None else DENSITY_THRESHOLD_DEFAULT
            window = args.density_window if args.density_window is not None else DENSITY_WINDOW_DEFAULT
            self.p(f"<h4>Density filtering enabled: SNPs removed when &gt;={threshold} SNPs "
                   f"found within {window} bp window</h4>")
        if self.all_vcf:
            grp = self.group
            all_name = getattr(grp, 'all_vcf_name', None) or 'all_vcf'
            built = (all_name in (getattr(grp, 'group_dataframe_dict', None) or {})
                     and all_name not in (getattr(grp, 'group_failures', None) or {}))
            if built:
                self.p("<h4>All_VCFs is available</h4>")
            else:
                # Said whenever -a was given, even when the all-samples group
                # produced no table (fewer than four samples, no informative
                # position, or a failure listed below).
                self.p(f"<h4>-a was given, but the all-samples group ({self.esc(all_name)}) "
                       f"produced no table or tree</h4>")
        self.tree_mode()

    def long_read_samples(self):
        # A nanopore VCF says so in its header (##vsnp3_platform), and its QUAL
        # is the evidence QUAL of vsnp3_longread.py rather than freebayes'.  The
        # thresholds read it the same way, but whoever reads the tables should
        # know which samples were sequenced on which platform, and how their
        # QUAL was made.
        stats = getattr(self.vcf_to_df, 'sample_stats', None) or {}
        longread = sorted((base, s) for base, s in stats.items()
                          if s.get('platform') == 'nanopore')
        if not longread:
            return
        self.p('<div class="box">')
        self.p(f'<b>Long-read (nanopore) samples: {len(longread)}</b><br>')
        self.p('These VCFs were called by bcftools from minimap2 alignments. Their QUAL '
               'is the vSNP3 evidence QUAL (vsnp3_longread.py): the Phred-scaled '
               'probability that the site is reference-only, from the read counts, the '
               'error rate measured in that sample and the variant prior. It is uncapped, '
               'on the freebayes scale, and the thresholds above read it exactly as they '
               'read freebayes QUAL. Calls supported on one strand only, or in or beside '
               'a homopolymer run, were limited to N. The model line of each sample:<br>')
        for base, s in longread:
            model = s.get('longread_model') or 'model line not recorded in the VCF header'
            self.p(f'&nbsp;&nbsp;<b>{self.esc(base)}</b>: '
                   f'<span class="small">{self.esc(model)}</span><br>')
        self.p('</div>')

    def tree_mode(self):
        # VSNP3_BOOTSTRAP changes which RAxML analysis builds every tree, and it
        # is set outside the command line (vsnp_gui sets it from its form), so
        # the report says which analysis ran and what the variable held.
        replicates, raw = bootstrap_replicates()
        if replicates:
            self.p(f'<div class="box"><b>RAxML trees: rapid bootstrap, {replicates} '
                   f'replicates</b> (VSNP3_BOOTSTRAP={self.esc(raw)}). Internal nodes carry '
                   f'bootstrap support. This is a different analysis from the default '
                   f'best-tree-only run, and its tree topology and table row order can '
                   f'differ from it.</div>')
        elif raw is None or not str(raw).strip():
            self.p('<div style="font-size:11px;">RAxML trees: best tree only '
                   '(VSNP3_BOOTSTRAP not set)</div>')
        else:
            self.p(f'<div style="font-size:11px;">RAxML trees: best tree only '
                   f'(VSNP3_BOOTSTRAP={self.esc(raw)}, read as 0)</div>')

    def tree_failures(self):
        # Tables and RAxML trees for these groups were written; only the HTML
        # view is missing.  Printed during the run as well.
        failures = getattr(self.group, 'html_tree_failures', None) or {}
        if failures:
            self.p('<div class="box">')
            self.p(f'<b>HTML tree not written: {len(failures)} group(s)</b><br>')
            self.p('The tables and RAxML tree of each group below were written; only '
                   'the HTML tree view is missing.<br>')
            for group_name, reason in sorted(failures.items()):
                self.p(f'&nbsp;&nbsp;<b>{self.esc(group_name)}</b>: {self.esc(reason)}<br>')
            self.p('</div>')

    def group_failures(self):
        # A group whose RAxML or tables failed was printed in the terminal
        # only, and the Groupings table listed its members as analysed.
        failures = getattr(self.group, 'group_failures', None) or {}
        if failures:
            self.p('<div class="box">')
            self.p(f'<b>No table or tree: {len(failures)} group(s)</b><br>')
            self.p('These groups were selected but produced no SNP table and no tree. '
                   'Their samples are analysed in every other group they belong to, and '
                   'the group is marked in the Groupings table below. The reason per group:<br>')
            for group_name, reason in sorted(failures.items()):
                self.p(f'&nbsp;&nbsp;<b>{self.esc(group_name)}</b>: {self.esc(reason)}<br>')
            self.p('</div>')

    def sample_accounting(self):
        acc = self.accounting
        vr = self.validation_results
        buckets = acc.by_bucket()
        analysed = len(acc.analysed)
        not_analysed = len(acc.excluded)
        self.p("<h2>Sample accounting</h2>")
        self.p('<div style="font-size:11px;">Every VCF file found is listed in exactly one '
               'row below. A sample that is analysed but has no SNP passing the thresholds '
               'is still in the tables and trees (All group, and "Group Not Found"), with '
               'the reference call at every position.</div><br>')
        self.p("<table style='width: 60%;'>")
        self.p("<tr style='background-color: #f2f2f2;'><th>Category</th><th>Count</th></tr>")
        self.p(f"<tr><td><strong>VCF files found</strong></td><td>{len(acc.files)}</td></tr>")
        self.p(f"<tr><td><strong>Analysed</strong></td><td class='ok'>{analysed}</td></tr>")
        self.p(f"<tr><td>&nbsp;&nbsp;of which, no SNPs passing thresholds</td>"
               f"<td class='{'warn' if acc.no_snps else ''}'>{len(acc.no_snps)}</td></tr>")
        where = ' <span style="font-weight:normal;">(each file is listed at the end of this report)</span>' if not_analysed else ''
        self.p(f"<tr><td><strong>Not analysed</strong>{where}</td>"
               f"<td class='{'bad' if not_analysed else ''}'>{not_analysed}</td></tr>")
        for bucket in SampleAccounting.BUCKETS:
            if buckets.get(bucket):
                self.p(f"<tr><td>&nbsp;&nbsp;{bucket}</td><td>{len(buckets[bucket])}</td></tr>")
        self.p("</table><br>")
        if acc.unaccounted:
            self.p(f"<div class='box'><b>{len(acc.unaccounted)} file(s) are not accounted for "
                   f"above. This is a defect in vSNP3; please report it.</b><br>"
                   + '<br>'.join(self.esc(f) for f in acc.unaccounted) + '</div>')

        if acc.no_snps:
            self.p(f"<h3 class='warn'>Analysed with no SNPs passing thresholds ({len(acc.no_snps)})</h3>")
            self.p("<div style='font-size:11px;'>These samples are in the tables and trees but "
                   "contribute no positions of their own. A sample identical to the reference "
                   "looks like this. So does a sample whose calls are all mixed (AC=1), which "
                   "usually means a mixed or contaminated sample: those are flagged.</div>")
            # Which group tables and trees each one is in, so the reader does
            # not have to infer it from the Groupings table further down.
            member_of = defaultdict(list)
            for group_name, members in (self.groupings_dict or {}).items():
                if group_name != 'Group Not Found':
                    for sample in members:
                        member_of[sample].append(group_name)
            self.p("<div class='small'>")
            for name, detail in sorted(acc.no_snps.items()):
                flag = ''
                if detail.get('mixed_only'):
                    flag = (" <b class='bad'>&#9888; every call is AC=1 (mixed): possible mixed "
                            "sample or contamination</b>")
                groups = member_of.get(name)
                if groups:
                    where = ('in the tables and tree of: '
                             + ', '.join(self.esc(g) for g in sorted(groups)))
                else:
                    where = ('in no table: no defining SNP group matched and -a was not '
                             'used, so there is no All table to carry it')
                self.p(f"&#9679; <b>{self.esc(name)}</b> ({self.esc(detail.get('file', ''))}) "
                       f"&mdash; {self.esc(detail.get('reason', ''))}{flag}<br>"
                       f"&nbsp;&nbsp;&nbsp;<i>{where}</i><br>")
            self.p("</div>")

        grp = self.group
        left_out = {g: smp for g, smp in (getattr(grp, 'undetermined_in_group', None) or {}).items() if smp}
        if left_out:
            total = sum(len(smp) for smp in left_out.values())
            self.p(f"<h3 class='warn'>Left out of a group's tree and table ({total})</h3>")
            self.p("<div style='font-size:11px;'>These samples belong to the group by its "
                   "defining SNP but have no determined base (every call is '-' or N) at the "
                   "group's positions, so RAxML cannot place them; they are left out of that "
                   "group's alignment and table. They are analysed in every other group they "
                   "belong to.</div>")
            self.p("<div class='small'>")
            for group_name, samples in sorted(left_out.items()):
                self.p(f"&#9679; <b>{self.esc(group_name)}</b>: "
                       + ', '.join(self.esc(x) for x in sorted(samples)) + "<br>")
            self.p("</div>")
        # Possible mixed samples are marked in the Groupings table, not listed
        # here.  The calls read from inside multi-base records are not shown in
        # the summary; the terminal and the tables carry them.

        self.p("<h3>VCF file validation</h3>")
        if vr:
            self.p("<table style='width: 60%;'>")
            self.p("<tr style='background-color: #f2f2f2;'><th>Validation</th><th>Count</th></tr>")
            self.p(f"<tr><td>Files checked</td><td>{vr.total_files}</td></tr>")
            self.p(f"<tr><td>Valid</td><td class='ok'>{vr.total_valid}</td></tr>")
            self.p(f"<tr><td>Not valid (listed under Not analysed)</td>"
                   f"<td class='{'bad' if vr.total_invalid else ''}'>{vr.total_invalid}</td></tr>")
            if vr.header_only_files:
                self.p(f"<tr><td>Valid, no records at all (reference-identical where covered)</td>"
                       f"<td>{len(vr.header_only_files)}</td></tr>")
            if vr.reference_contigs:
                shown = ', '.join(vr.reference_contigs[:6])
                if len(vr.reference_contigs) > 6:
                    shown += f', ... ({len(vr.reference_contigs)} contigs)'
                self.p(f"<tr><td>Reference contigs (shared by {vr.reference_file_count} files)</td>"
                       f"<td>{self.esc(shown)}</td></tr>")
            self.p(f"<tr><td>Validation time</td><td>{vr.seconds:.1f} s</td></tr>")
            self.p("</table>")
            non_ascii = [b for b, s in self.vcf_to_df.sample_stats.items() if s.get('non_ascii_bytes')]
            if non_ascii:
                total = sum(self.vcf_to_df.sample_stats[b]['non_ascii_bytes'] for b in non_ascii)
                self.p(f"<div style='font-size:11px;'>{len(non_ascii)} file(s) contained "
                       f"{total} non-ASCII byte(s), usually in a header path; ignored when read. "
                       f"The files were not changed. Listed in the validation log.</div>")
            odd = [(b, s) for b, s in self.vcf_to_df.sample_stats.items()
                   if s.get('malformed') or s.get('skipped_lines')]
            if odd:
                self.p("<div style='font-size:11px;'>Records that could not be read and were "
                       "skipped (the rest of the file was analysed):</div><div class='small'>")
                for base, s in odd:
                    self.p(f"{self.esc(base)}: {s.get('malformed', 0) + s.get('skipped_lines', 0)} "
                           f"record(s)<br>")
                self.p("</div>")
            if vr.total_invalid == 0 and not vr.header_only_files:
                self.p("<div class='ok' style='font-size:11px;'>&#9989; All VCF files passed validation.</div>")
        self.p("<br>")

    def groupings(self):
        group_vcfs_dict = defaultdict(list) #invert the key, values
        for group, dataframes in self.groupings_dict.items():
            for vcf in dataframes.keys():
                group_vcfs_dict[vcf].append(group)
        group_vcfs_dict = dict(sorted(group_vcfs_dict.items())) #sorts on key(vcf name)

        failed = getattr(self.group, 'group_failures', None) or {}
        # Possible mixed samples (vsnp3_version.mixed_signature): the name is
        # highlighted and the counts follow the sample's groups.  Analysed as usual.
        mixed = getattr(self.group, 'mixed_samples', None) or {}
        self.p(f'<h2>Groupings with {len(group_vcfs_dict):,} listed:</h2>')
        if mixed:
            q = getattr(self.group, 'qual_threshold', QUAL_THRESHOLD)
            mq = getattr(self.group, 'mq_threshold', MQ_THRESHOLD)
            self.p(f"<div style='font-size:11px;'><span class='mixed'>Highlighted</span>: "
                   f"possible mixed sample ({len(mixed)}), with more than {MIXED_AC1_RATIO:g} "
                   f"confident mixed calls (AC=1, QUAL &ge; {q}, MQ &ge; {mq}) per passing SNP. "
                   f"Analysed as usual; review before relying on its place in the tables and "
                   f"trees.</div>")
        self.p("<table>")
        self.p("<tr align=\"left\"><th>Sample Name</th><tr>")
        for key, value in group_vcfs_dict.items():
            flag = mixed.get(key)
            if flag:
                row = [f"<td><span class='mixed'>{self.esc(key)}</span></td>"]
            else:
                row = [f"<td>{self.esc(key)}</td>"]
            for group in value:
                if group == "Group Not Found":
                    row.append(f'<td><span style="color: red;">{self.esc(group)}</span></td>')
                elif group in failed:
                    row.append(f'<td><span style="color: red;">{self.esc(group)} (no table or tree)</span></td>')
                else:
                    row.append(f"<td>{self.esc(group)}</td>")
            if flag:
                row.append(f"<td class='bad'>{flag.get('confident_ac1', 0):,} confident mixed calls "
                           f"against {flag.get('passing', 0):,} passing SNPs</td>")
            self.p("<tr>" + "\t".join(row) + "</tr>")
        self.p("</table>")

    def not_analysed(self):
        '''
        Every VCF that was not analysed, by reason, each marked with a red X.

        Placed after the Groupings table, where 3.36 listed the files removed
        by remove_from_analysis.xlsx.  That list was a duplicate of this one
        (and silent about the other reasons a file is left out), so it is
        gone.  The counts stay in the Sample accounting table at the top.
        '''
        acc = self.accounting
        buckets = acc.by_bucket()
        not_analysed = len(acc.excluded)
        if not not_analysed:
            return
        self.p(f"<h2 class='bad'>Not analysed ({not_analysed})</h2>")
        self.p('<div style="font-size:11px;">These VCF files were found but are in no table '
               'or tree. Counts by reason are in the Sample accounting table at the top of '
               'this report.</div>')
        for bucket in SampleAccounting.BUCKETS:
            items = buckets.get(bucket)
            if not items:
                continue
            self.p(f"<h4>{bucket} ({len(items)})</h4>")
            self.p("<div class='small'>")
            for base, reason in items:
                self.p(f"❌ <b>{self.esc(base)}</b> &mdash; {self.esc(reason)}<br>")
            self.p("</div>")
        self.p("<br>")

    def write_json(self):
        '''
        The sample accounting as JSON beside the summary, so a caller (the
        vSNP GUI reconciles what it staged against what was analysed) can read
        it instead of scraping the HTML.
        '''
        acc = self.accounting
        grp = self.group
        groupings = defaultdict(list)
        for group_name, members in (self.groupings_dict or {}).items():
            for sample in members:
                groupings[sample].append(group_name)
        data = {
            'vsnp3_version': __version__,
            'reference_type': str(self.reference),
            'summary_html': os.path.basename(self.path),
            'started': self.timer.started_at.strftime('%Y-%m-%d %H:%M:%S'),
            'vcf_files_found': list(acc.files),
            'analysed': {name: base for name, base in sorted(acc.analysed.items())},
            'analysed_no_snps_passing': {name: detail for name, detail in sorted(acc.no_snps.items())},
            'not_analysed': {base: {'reason_class': bucket, 'reason': reason}
                             for base, (bucket, reason) in sorted(acc.excluded.items())},
            'unaccounted': list(acc.unaccounted),
            'remove_list_entries_unmatched': list(getattr(acc, 'remove_unmatched', None) or []),
            'duplicate_sample_names': self.name_collisions,
            'sample_names_changed': getattr(grp, 'fallback_renames', None) or {},
            'left_out_of_group': {g: sorted(smp) for g, smp in
                                  (getattr(grp, 'undetermined_in_group', None) or {}).items() if smp},
            'possible_mixed_samples': getattr(grp, 'mixed_samples', None) or {},
            'groups_without_table_or_tree': getattr(grp, 'group_failures', None) or {},
            'html_tree_not_written': getattr(grp, 'html_tree_failures', None) or {},
            'groupings': {sample: sorted(groups) for sample, groups in sorted(groupings.items())},
        }
        with open(self.json_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=1, default=str, sort_keys=False)

    def timing(self):
        timer = self.timer
        finished = datetime.now()
        self.p("<hr>")
        self.p("<h2>Run time</h2>")
        self.p("<table style='width: 60%;'>")
        self.p(f"<tr><td><b>Started</b></td><td>{timer.started_at.strftime('%Y-%m-%d %H:%M:%S')}</td></tr>")
        self.p(f"<tr><td><b>Finished</b></td><td>{finished.strftime('%Y-%m-%d %H:%M:%S')}</td></tr>")
        self.p(f"<tr><td><b>Total run time</b></td><td>{format_elapsed(timer.elapsed)}</td></tr>")
        self.p("</table><br>")
        total = max(timer.elapsed, 1e-9)
        phases = list(timer.phases)
        # Whatever is not in a named phase: this report, the version lookup,
        # file deletion.  Listed so the rows add up to the total.
        other = total - sum(seconds for _, seconds in phases)
        if other > 0.5:
            phases.append(('Summary report and cleanup', other))
        self.p("<table style='width: 60%;'>")
        self.p("<tr style='background-color: #f2f2f2;'><th>Phase</th><th>Time</th><th>Share</th></tr>")
        for name, seconds in phases:
            self.p(f"<tr><td>{self.esc(name)}</td><td>{format_elapsed(seconds)}</td>"
                   f"<td>{seconds * 100 / total:.0f}%</td></tr>")
        self.p("</table>")

    def system_information(self):
        try:
            import platform
            self.p("\n<h2>System Information:</h2>")
            os_name = platform.system()
            os_version = platform.version()
            os_release = platform.release()
            arch = platform.machine()
            processor = platform.processor()
            self.p(f"<b>Operating System:</b> {self.esc(os_name)} {self.esc(os_release)} {self.esc(os_version)}<br>")

            if os_name == 'Darwin':  # macOS
                cpu_type = "ARM (Apple Silicon)" if arch == 'arm64' else "Intel"
                mac_ver = platform.mac_ver()
                self.p(f"<b>macOS Details:</b> macOS {mac_ver[0]}, {cpu_type}<br>")
            elif os_name == 'Linux':
                try:
                    import distro
                    linux_distro = distro.name(pretty=True)
                except ImportError:
                    linux_distro = "Unknown Linux Distribution"
                    try:
                        with open('/etc/os-release') as f:
                            for line in f:
                                if line.startswith('PRETTY_NAME='):
                                    linux_distro = line.split('=', 1)[1].strip().strip('"')
                                    break
                    except OSError:
                        pass
                hpc_info = None
                for env_var, hpc_type in (('SLURM_CLUSTER_NAME', 'SLURM'), ('PBS_HOME', 'PBS'),
                                          ('SGE_ROOT', 'SGE'), ('LSB_JOBID', 'LSF')):
                    if env_var in os.environ:
                        hpc_info = f"{hpc_type} HPC environment"
                        break
                if hpc_info is None:
                    for path, hpc_type in (('/opt/slurm', 'SLURM'), ('/opt/pbs', 'PBS'),
                                           ('/opt/sge', 'SGE'), ('/opt/lsf', 'LSF')):
                        if os.path.exists(path):
                            hpc_info = f"{hpc_type} HPC environment"
                            break
                if hpc_info:
                    self.p(f"<b>Linux Details:</b> {self.esc(linux_distro)}, {hpc_info}<br>")
                else:
                    self.p(f"<b>Linux Details:</b> {self.esc(linux_distro)}<br>")
            elif os_name == 'Windows':
                win_ver = platform.win32_ver()
                self.p(f"<b>Windows Details:</b> Windows {win_ver[0]} (Build {win_ver[1]})<br>")

            self.p(f"<b>CPU Architecture:</b> {self.esc(arch)}<br>")
            self.p(f"<b>Processor:</b> {self.esc(processor)}<br>")
            self.p(f"<b>CPU count:</b> {multiprocessing.cpu_count()}<br>")
            try:
                import cpuinfo
                cpu_info = cpuinfo.get_cpu_info()
                self.p(f"<b>CPU Model:</b> {self.esc(cpu_info['brand_raw'])}<br>")
            except Exception:                                   # noqa: BLE001
                pass
        except Exception as e:                                  # noqa: BLE001
            self.p(f"Error collecting system information: {self.esc(e)} <br>")

    def program_versions(self):
        self.p("<hr>")
        self.p("\n<h2>Program versions:</h2>")
        self.p(f'vSNP3: {__version__} <br>')
        self.p(f'Python: {self.esc(sys.version)} <br>')
        try:
            conda_versions = _conda_versions()
            conda_env = os.environ.get('CONDA_DEFAULT_ENV')
            if conda_env:
                self.p(f'<b>Versions from conda environment: {self.esc(conda_env)}</b> <br>')
            else:
                self.p('<b>Versions from the Python environment</b> <br>')
            program_list = ['biopython', 'dask', 'humanize', 'numpy', 'pandas', 'openpyxl',
                            'xlsxwriter', 'parallel', 'pigz', 'regex', 'py-cpuinfo', 'plotly']
            for program in program_list:
                if program.lower() in conda_versions:
                    self.p(f'{program}: {self.esc(conda_versions[program.lower()])} (conda) <br>')
                    continue
                try:
                    version = importlib.metadata.version(program)
                    self.p(f'{program}: {self.esc(version)} (python package) <br>')
                except importlib.metadata.PackageNotFoundError:
                    self.p(f'{program}: nd <br>')
            if self.raxml_version:
                self.p(f'raxml: {self.esc(self.raxml_version)} <br>')
            elif 'raxml' in conda_versions:
                self.p(f'raxml: {self.esc(conda_versions["raxml"])} (conda) <br>')
            else:
                self.p('raxml: nd <br>')
        except Exception as e:                                  # noqa: BLE001
            self.p(f"Error checking versions: {self.esc(e)} <br>")


if __name__ == "__main__": # execute if directly access by the interpreter
    # Set multiprocessing start method here to be compatible with Python 3.12
    multiprocessing.set_start_method('spawn', True)

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Store VCF files from vSNP step1 to step 2 directory.  VCF files must be stored by reference type.  Make a VCF file directory database that will build over time as samples are ran in step 1

    For example...
    <path/to/files>
        referenceA_dir
            step1_dir
                sample1_dir
                    <alignemnt_files>
                sample2_dir
                    <alignemnt_files>
            step2_dir
                vcf_source_dir
                    sample1.vcf
                    sample2.vcf
                comparison1_dir
                comparison2_dir

    <path/to/dependencies> (added path using vsnp3_path_adder.py)
        referenceA_dir
            defining_snps_for_referenceA.xlsx
            metadata_for_referenceA.xlsx
            FASTA/s for referenceA
            GBK/s for referenceA

    When running samples through step 1 and 2 of vSNP, or when running a routine analysis, set up dependencies using vsnp3_path_adder.py.  See vsnp3_path_adder.py -h for adding a reference type and for more information

    Usage:

    vsnp3_step2.py -a -d -t ASFV_Georgia_2007

    vsnp3_step2.py -wd <path/to/vcf_directory> -abs_pos chrom1:123456 -group test_group -m <path/to/metadata.xlsx>

    vsnp3_step2.py -wd ../vcf_source

    vsnp3_step2.py -a --remove

    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-wd', '--wd', action='store', dest='wd', required=False, default='.', help='Optional: path to VCF files. By default .vcf in current working directory are used.')
    parser.add_argument('-o', '--output', action='store', dest='output_dir', required=False, default=None, help="Optional: Provide a name.  This name will be a directory output files are writen to.  Name can be a directory path, but doesn't have to be. By default VCF files are worked on in your current working directory")
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', default=None, required=False, help='Optional: A valid reference_type name will be automatically found, but a valid reference_type name can be supplied.  See vsnp3_path_adder.py -s')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple gbk files can be specified with wildcard')
    parser.add_argument('-s', '--defining_snps', action='store', dest='defining_snps', default=None, required=False, help='Optional: Defining SNPs with positions to filter.  See template_define_filter.xlsx in vsnp dependency folder.  Recommended having this file in reference type folder')
    parser.add_argument('-m', '--metadata', action='store', dest='metadata', default=None, required=False, help='Optional: Two column Excel file, Column One: full VCF file name, Column Two: Updated name.  Recommended having this file in reference type folder')
    parser.add_argument('-remove_by_name', '--remove_by_name', action='store', dest='remove_by_name', required=False, help='Optional: Excel file containing samples to remove from analysis Column 1: to match sample name minus extension. No header allowed.   Recommended having this file in reference type folder')
    parser.add_argument('-n', '--no_filters', action='store_true', dest='no_filters', default=False, help='Optional: turn off filters')
    parser.add_argument('-w', '--qual_threshold', action='store', dest='qual_threshold', default=QUAL_THRESHOLD, required=False, help='Optional: Minimum QUAL threshold for calling a SNP')
    parser.add_argument('-x', '--n_threshold', action='store', dest='n_threshold', default=N_THRESHOLD, required=False, help='Optional: Minimum N threshold.  SNPs between this and qual_threshold are reported as N')
    parser.add_argument('-y', '--mq_threshold', action='store', dest='mq_threshold', default=MQ_THRESHOLD, required=False, help='Optional: At least one position per group must have this minimum MQ threshold to be called.')
    parser.add_argument('-f', '--fix_vcfs', action='store_true', dest='fix_vcfs', help='Optional: Write normalized copies of the VCF files to "vcf_normalized" and exit.  The input files are not modified.')
    parser.add_argument('--assume_gt_only_quality', action='store_true', dest='assume_gt_only_quality', default=False, help=f'Optional: for records whose FORMAT is GT only, and which therefore carry no quality data, substitute QUAL={SYNTHESIZED_QUAL} so they are not rejected by every threshold.  Genotypes and INFO are left untouched and the substitution is recorded in the VCF header.  For assembly or consensus derived VCFs.')
    parser.add_argument('--first_sample_only', action='store_true', dest='first_sample_only', default=False, help='Optional: accept multi-sample VCF files by analysing only their first sample.  Without this a multi-sample VCF is refused rather than silently reduced to sample one.')
    parser.add_argument('-k', '--keep_ind_vcfs', action='store_true', dest='keep_ind_vcfs', default=False, help='Optional: keep the VCF files that were read.  By default they are deleted once archived to "vcf_starting_files.zip", since leaving them clutters the working directory.  Files outside the working directory that this run did not copy are never deleted regardless, so pointing --wd at a shared VCF database cannot empty it.')
    parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', required=False, help='Optional: create table with all isolates')
    parser.add_argument('-i', '--find_new_filters', action='store_true', dest='find_new_filters', help='Optional: find new positions to apply to the filter file.  Positions must be manually added to filter file.  They are not added by running this command.  Only text files are output showing position detail. Curant before adding filters')
    parser.add_argument('-abs_pos', '--abs_pos', action='store', dest='abs_pos', required=False, help='Optional: Make a group on defining SNP.  Must be supplied with --group option.  Format as chrom in VCF, chrom:10000.')
    parser.add_argument('-group', '--group', action='store', dest='group', required=False, help='Optional: Name a group on defining SNP.  Must be supplied with --abs_pos option')
    parser.add_argument('-hash', '--hash_groups', action='store_true', dest='hash_groups', required=False, help='Optional: The option will run defining snps marked with a # in the defining snps file.  The # is removed and the defining snps are run.')
    parser.add_argument('--show_groups', action='store_true', dest='show_groups', help='Show group names in SNP table')
    parser.add_argument('-html_tree', '--html_tree', action='store_true', dest='html_tree', help='Optional: Generate HTML tree visualization (automatically enables -dp)')
    parser.add_argument('-dp', '--dp', action='store_true', dest='dp', help='Optional: Include average depth of coverage in tables')
    parser.add_argument('--density_threshold', nargs='?', const=DENSITY_THRESHOLD_DEFAULT, type=int, dest='density_threshold', help='Optional: Minimum number of SNPs required to trigger density filtering (default: 3)')
    parser.add_argument('--density_window', nargs='?', const=DENSITY_WINDOW_DEFAULT, type=int, dest='density_window', help='Optional: Window size in base pairs for density filtering (default: 20)')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling.  A pickle file will be kept for troubleshooting to be used directly in vsnp3_group_on_defining_snps.py.  This saves processing time')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    # Handle density filtering parameter defaults
    # If only one parameter is provided, set the other to its default
    if args.density_threshold is not None and args.density_window is None:
        args.density_window = 20
        print(f"Density filtering enabled with threshold={args.density_threshold}, using default window={args.density_window}")
    elif args.density_window is not None and args.density_threshold is None:
        args.density_threshold = 3
        print(f"Density filtering enabled with window={args.density_window}, using default threshold={args.density_threshold}")
    elif args.density_threshold is not None and args.density_window is not None:
        print(f"Density filtering enabled with custom threshold={args.density_threshold}, window={args.density_window}")

    # Determine if density filtering should be enabled
    filter_density = args.density_threshold is not None or args.density_window is not None

    setup = Setup(debug=args.debug)
    timer = PhaseTimer()
    global_date_stamp = setup.date_stamp
    global_working_dir = setup.cwd
    # WORKING DIRECTORY AND VCF LIST VALIDATION
    print("\n🔍 VALIDATING WORKING DIRECTORY AND VCF FILES...")
    print("="*60)

    validator = InputValidator(debug=args.debug)
    cwd_test = False

    if args.wd == '.':
        cwd_test = True
        working_dir = os.getcwd()
        print(f"Using current directory: {working_dir}")
        vcf_list = sorted(glob.glob('*vcf'))
        wd_vcf_list = vcf_list
    else:
        # Validate and process specified working directory
        wd = os.path.expanduser(args.wd)
        wd = os.path.abspath(wd)

        # Validate directory exists and is accessible
        is_dir_valid, dir_msg = validator.validate_directory_exists(wd)
        if not is_dir_valid:
            print(f"\n❌ WORKING DIRECTORY VALIDATION FAILED:")
            print(f"   {dir_msg}")
            print(f"\n💡 SOLUTION:")
            print(f"   • Check that directory exists: {wd}")
            print(f"   • Verify you have read permissions")
            print(f"   • Use absolute path if needed")
            print("="*60)
            sys.exit(1)

        print(f"✅ Working directory validated: {wd}")
        vcf_pattern = os.path.join(wd, '*vcf')
        vcf_list = sorted(glob.glob(vcf_pattern))
        wd_vcf_list = vcf_list

    # Validate VCF file list is not empty
    if not vcf_list:
        print(f"\n❌ NO VCF FILES FOUND!")
        if args.wd == '.':
            search_location = "current directory"
        else:
            search_location = f"directory: {wd}"

        print(f"   Searched in: {search_location}")
        print(f"   Pattern: *vcf (any file name ending in vcf)")
        print(f"\n💡 SOLUTIONS:")
        print(f"   • Verify VCF files exist in the directory")
        print(f"   • Check file extensions (.vcf)")
        print(f"   • Ensure vsnp3_step1.py has been run first")
        print(f"   • Check directory permissions")
        print("="*60)
        sys.exit(1)

    print(f"✅ Found {len(vcf_list)} VCF files for processing")
    if args.debug:
        print("VCF files found:")
        for vcf_file in vcf_list[:10]:  # Show first 10 files
            print(f"   {os.path.basename(vcf_file)}")
        if len(vcf_list) > 10:
            print(f"   ... and {len(vcf_list) - 10} more files")

    print("="*60)

    if args.output_dir:
        wd_vcf_list = []
        output_dir = args.output_dir
        output_dir = os.path.expanduser(output_dir)
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        for each_vcf in vcf_list:
            shutil.copy(each_vcf, output_dir)
            wd_vcf_list.append(os.path.join(output_dir, os.path.basename(each_vcf)))
        os.chdir(output_dir)
        setup.cwd = os.getcwd()
        global_working_dir = setup.cwd
        timer.mark('Copy VCFs to output directory')

    accounting = SampleAccounting(wd_vcf_list)

    vcf_to_df = VCF_to_DF(vcf_list=wd_vcf_list, debug=args.debug,
                          assume_gt_only_quality=args.assume_gt_only_quality,
                          first_sample_only=args.first_sample_only,
                          write_normalized=args.fix_vcfs, timer=timer)
    for base, (bucket, reason) in vcf_to_df.excluded.items():
        accounting.exclude(base, bucket, reason)
    if args.fix_vcfs:
        print(f'\nNormalized copies written to {vcf_to_df.normalized_dir}')
        print('The input VCF files were not modified.')
        if vcf_to_df.excluded:
            print(f'\n{len(vcf_to_df.excluded)} file(s) could not be normalized:')
            for each in vcf_to_df.vcf_bad_list:
                print(f'  {each}')
            sys.exit(1)
        sys.exit(0)
    if not args.debug:
        shutil.rmtree(vcf_to_df.normalized_dir, ignore_errors=True)

    # Create the file to indicate the script is running
    notification_file = "step2_is_running__individual_folders_may_be_complete"
    with open(notification_file, 'w') as f:
        f.write("Script is still running.")
    print(f"Created file: {notification_file}")

    print('\nFiles excluded before analysis:' if vcf_to_df.excluded else '\nNo files excluded before analysis.')
    for each in vcf_to_df.vcf_bad_list:
        print(f'\t {each}')

    if args.reference_type:
        ro = Ref_Options(args.reference_type)
    else:
        ro = Ref_Options(vcf_to_df.chrom)

    if args.abs_pos and not args.group:
        print('\n### -abs_pos must be used with -group option\n')
        sys.exit(1)
    if args.group and not args.abs_pos:
        print('\n### -group must be used with -abs_pos option\n')
        sys.exit(1)

    # Prioritize explicitly provided files over reference type defaults
    # Only use reference type defaults if the specific arguments were not provided
    # This makes sure explicitly provided files via -b, -s, -m, and -remove_by_name take precedence
    if not args.gbk and ro.gbk:
        args.gbk = ro.gbk
        print(f"Using reference type GBK: {args.gbk}")
    elif args.gbk:
        print(f"Using explicitly provided GBK: {args.gbk}")

    if not args.defining_snps and ro.defining_snps:
        args.defining_snps = ro.defining_snps
        print(f"Using reference type defining SNPs: {args.defining_snps}")
    elif args.defining_snps:
        print(f"Using explicitly provided defining SNPs: {args.defining_snps}")

    if not args.metadata and ro.metadata:
        args.metadata = ro.metadata
        print(f"Using reference type metadata: {args.metadata}")
    elif args.metadata:
        print(f"Using explicitly provided metadata: {args.metadata}")
    else:
        print("No metadata file for this reference type; VCF file names are used as sample names")

    if not args.remove_by_name and ro.remove:
        args.remove_by_name = ro.remove
        print(f"Using reference type remove list: {args.remove_by_name}")
    elif args.remove_by_name:
        print(f"Using explicitly provided remove list: {args.remove_by_name}")

    print(f'Before sample filter: {len(vcf_to_df.dataframes)}')
    if args.remove_by_name:
        remove_from_analysis = Remove_From_Analysis(working_directory=global_working_dir, excel_remove=args.remove_by_name, extension="vcf")
        for each in remove_from_analysis.remove_list:
            sample_basename = os.path.basename(each)
            # Only samples that are actually present; the remove list names
            # several spellings of each entry (with/without .vcf, with _zc).
            if sample_basename in vcf_to_df.dataframes:
                vcf_to_df.dataframes.pop(sample_basename, None)
                vcf_to_df.sample_stats.pop(sample_basename, None)
                accounting.exclude(sample_basename, 'removed by remove_from_analysis',
                                   f'listed in {os.path.basename(args.remove_by_name)}')
        # An entry that named no file is said, in the terminal and the accounting
        # JSON (not the summary HTML).  It used to be silent, so a misspelt
        # remove-list entry left the sample in every tree with nothing to show
        # for the intent.
        found = set(accounting.files)
        accounting.remove_unmatched = [
            entry for entry in getattr(remove_from_analysis, 'entries', [])
            if not any(spelling in found
                       for spelling in (entry, f'{entry}.vcf', f'{entry}_zc.vcf'))]
        if accounting.remove_unmatched:
            shown = ', '.join(accounting.remove_unmatched[:20])
            more = ' ...' if len(accounting.remove_unmatched) > 20 else ''
            print(f'\n### {len(accounting.remove_unmatched)} entr'
                  f'{"y" if len(accounting.remove_unmatched) == 1 else "ies"} in '
                  f'{os.path.basename(args.remove_by_name)} matched no VCF file: {shown}{more}')
    print(f'After sample filter: {len(vcf_to_df.dataframes)}')

    # Archive the input exactly as read, straight from the source files.  3.36
    # copied every VCF to vcf_starting_files/, zipped the directory, then
    # deleted it: a full extra write and delete of the whole input set.
    starting_files = os.path.join(setup.cwd, 'vcf_starting_files')
    with zipfile.ZipFile(f'{starting_files}.zip', 'w', zipfile.ZIP_DEFLATED) as zf:
        for each_vcf in wd_vcf_list:
            if os.path.exists(each_vcf):
                zf.write(each_vcf, os.path.basename(each_vcf))
        if args.defining_snps:
            zf.write(args.defining_snps, os.path.basename(args.defining_snps)) #package with starting files for the record
    timer.mark('Archive starting files')


    # Before creating the Group instance, check if html_tree is True and set dp accordingly
    if args.html_tree:
        args.dp = True

    try:
        group = Group(cwd=global_working_dir, metadata=args.metadata, defining_snps=args.defining_snps, excel_remove=args.remove_by_name, gbk_list=args.gbk, dataframes=vcf_to_df.dataframes, all_vcf=args.all_vcf, find_new_filters=args.find_new_filters, no_filters=args.no_filters, qual_threshold=int(args.qual_threshold), n_threshold=int(args.n_threshold), mq_threshold=int(args.mq_threshold), abs_pos=args.abs_pos, group=args.group, show_groups=args.show_groups, hash_groups=args.hash_groups, html_tree=args.html_tree, dp=args.dp, filter_density=filter_density, density_threshold=args.density_threshold, density_window=args.density_window, debug=args.debug, sample_stats=vcf_to_df.sample_stats)
    except DefiningSnpsError as e:
        print(f'\n### {e}\n')
        sys.exit(1)
    timer.absorb(group.phase_times)

    # Account for what the grouping stage decided.
    for name, info in group.name_collisions.items():
        for dropped in info['dropped']:
            accounting.exclude(dropped, 'duplicate sample name',
                               f'resolves to "{name}", analysed from {info["used"]}: {info["cause"]}')
    accounting.analysed = dict(group.sample_sources)
    accounting.no_snps = dict(group.no_snp_samples)
    if accounting.unaccounted:
        print(f'\n### {len(accounting.unaccounted)} VCF file(s) are not accounted for: '
              f'{", ".join(accounting.unaccounted)}. This is a defect in vSNP3; please report it.')
    if accounting.no_snps:
        print(f'\n{len(accounting.no_snps)} sample(s) had no SNP passing the thresholds and are '
              f'included with the reference call at every position:')
        for name, detail in sorted(accounting.no_snps.items()):
            flag = '  <-- every call is AC=1: possible mixed sample' if detail.get('mixed_only') else ''
            print(f'\t {name}: {detail["reason"]}{flag}')
    if group.mixed_samples:
        print(f'\n{len(group.mixed_samples)} sample(s) look like a mixture of strains (more than '
              f'{MIXED_AC1_RATIO:g} confident mixed calls per passing SNP); analysed as usual and '
              f'highlighted in the Groupings table of the summary:')
        for name, detail in sorted(group.mixed_samples.items()):
            print(f'\t {name}: {detail["confident_ac1"]:,} confident mixed calls, '
                  f'{detail["passing"]:,} passing SNPs')
    if any(group.embedded_cells.values()) or any(group.embedded_now_uninformative.values()):
        print(f'\n{sum(group.embedded_cells.values()):,} table cell(s) were read from inside a '
              f'multi-base record; {sum(len(p) for p in group.embedded_now_uninformative.values()):,} '
              f'group position(s) left a table because that made every sample agree.')

    setup.print_time()
    HTML_Summary(timer=timer, vcf_to_df=vcf_to_df, reference=ro.select_ref, groupings_dict=group.groupings_dict, raxml_version=group.raxml_version, all_vcf_boolen=args.all_vcf, args=args, validation_results=vcf_to_df.validation_results, name_collisions=group.name_collisions, metadata_ambiguous=group.metadata_ambiguous, accounting=accounting, metadata_problem=group.metadata_problem, group=group, reference_note=ro.detection_note)
    print(f'Total run time: {format_elapsed(timer.elapsed)}')

    # Done last, once the summary is written.  Until 3.37 the files were
    # deleted here, before the grouping stage, so a failure there (a worksheet
    # mistake, say) left a working directory holding only the zip.
    #
    # The VCF files that were read are archived to vcf_starting_files.zip, so
    # leaving the loose copies behind just clutters the working directory --
    # deleting them is the useful default.  What must never happen is deleting
    # someone else's files: wd_vcf_list is the caller's own list unless
    # --output_dir made copies first, which is how `vsnp3_step2.py -wd
    # /path/to/vcf_database` used to empty a shared database.  -k/--keep_ind_vcfs
    # existed to prevent that and was parsed but never read.
    #
    # So: delete by default, but only files this run owns, and only now that the
    # archive holding them is complete.
    if args.keep_ind_vcfs:
        print(f'\n-k/--keep_ind_vcfs: leaving {len(wd_vcf_list)} VCF file(s) in place.')
    else:
        # Ours if we copied them ourselves, or if they sit inside the working
        # directory the user pointed us at.
        owned = bool(args.output_dir)
        if not owned:
            cwd_real = os.path.realpath(os.getcwd())
            owned = all(os.path.realpath(v).startswith(cwd_real + os.sep)
                        for v in wd_vcf_list)
        if owned:
            for each_vcf in wd_vcf_list:
                try:
                    os.remove(each_vcf)
                except FileNotFoundError:
                    # already gone, eg. it was empty and was dropped earlier
                    pass
        else:
            print(f'\nThe {len(wd_vcf_list)} VCF file(s) read are outside the working '
                  'directory and were not copied by this run, so they are left in place '
                  'rather than deleted -- they may be a shared database. A copy is '
                  f'archived in {os.path.basename(starting_files)}.zip. Use --output_dir '
                  'to work on copies, or delete them yourself.')

    try:
        os.remove(notification_file)
        print(f"Deleted file: {notification_file}")
    except FileNotFoundError:
        pass

    # A run in which no group produced a table or tree has not succeeded, whatever
    # it managed to write along the way.  Exiting 0 here is what let a broken
    # install look healthy to every wrapper script and scheduler.  A run where only
    # some groups failed still exits 0: the report names them, and the groups that
    # completed are valid.
    if group.group_failures and not group.raxml_version:
        print(f'\n### No group produced a table or tree. See the {len(group.group_failures)} '
              f'failure(s) reported above.')
        sys.exit(1)
# Created 2021 by Tod Stuber
