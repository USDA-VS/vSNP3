#!/usr/bin/env python

from vsnp3_version import __version__

import os
import sys
import io
import shutil
import subprocess
import re
import pickle
import locale
import argparse
import textwrap
import pandas as pd
import numpy as np
import zipfile
import glob
from datetime import datetime

from collections import defaultdict
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
from vsnp3_file_setup import Setup
from vsnp3_file_setup import redact_paths
from vsnp3_group_on_defining_snps import Group
from vsnp3_reference_options import Ref_Options
from vsnp3_remove_from_analysis import Remove_From_Analysis
from vsnp3_input_validator import InputValidator, VCFValidationResults

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

global_date_stamp=None
global_working_dir='.'

class VCF_to_DF():
    '''
    Enhanced VCF processing with comprehensive validation and error reporting
    '''

    def __init__(self, vcf_list=None, debug=False, assume_gt_only_quality=False,
                 first_sample_only=False): #write_out=False,
        '''
        Start at class call with comprehensive VCF validation
        '''
        self.startTime = datetime.now()
        self.vcf_bad_list = []
        self.vcf_original_count = len(vcf_list)
        cpu_count = int(multiprocessing.cpu_count() / 1.2)
        dataframes = {}

        # Normalization used to rewrite the caller's VCF files in place.  It now
        # writes copies here and the originals are never touched.  Set before the
        # process pool starts, because check_and_fix reads these off self.
        self.assume_gt_only_quality = assume_gt_only_quality
        self.first_sample_only = first_sample_only
        self.normalized_dir = os.path.join(os.getcwd(), 'vcf_normalized')
        os.makedirs(self.normalized_dir, exist_ok=True)

        # Initialize input validator for comprehensive VCF validation
        validator = InputValidator(debug=debug)

        # COMPREHENSIVE VCF VALIDATION - EARLY DETECTION
        print(f"\n🔍 VALIDATING {self.vcf_original_count} VCF FILES...")
        print("="*60)

        # Perform comprehensive validation of all VCF files
        validation_results = validator.validate_vcf_list(vcf_list)
        self.validation_results = validation_results  # Store for HTML reporting

        # Print validation summary to terminal
        validator.print_vcf_validation_summary(validation_results)

        # Write detailed validation log
        log_file = f'vcf_validation_log-{global_date_stamp}.txt'
        validator.write_validation_log(log_file, validation_results)
        print(f"\n📋 Detailed validation log written to: {log_file}")

        # CRITICAL ERROR HANDLING - REFERENCE MISMATCHES
        if validation_results.reference_mismatches:
            print("\n" + "="*80)
            print("🔥💥 CRITICAL ERROR: VCF REFERENCE MISMATCH DETECTED! 💥🔥")
            print("="*80)
            print("\nThe following VCF files use different reference genomes:")
            print("-" * 60)

            for mismatch in validation_results.reference_mismatches:
                file_name = os.path.basename(mismatch['file'])
                print(f"❌ {file_name}")
                print(f"   {mismatch['message']}")

            print("\n💡 SOLUTION:")
            print("   • All VCF files must use the same reference genome")
            print("   • Remove or re-process the mismatched files")
            print("   • Check your vSNP3 step1 reference settings")
            print("\n📋 Full details in validation log: " + log_file)
            print("="*80 + "\n")

            # Also log to validation log before exiting
            with open(log_file, 'a') as f:
                f.write("\n" + "="*50 + "\n")
                f.write("CRITICAL ERROR: REFERENCE MISMATCH\n")
                f.write("="*50 + "\n")
                for mismatch in validation_results.reference_mismatches:
                    f.write(f"File: {mismatch['file']}\n")
                    f.write(f"Error: {mismatch['message']}\n")
                    f.write("-" * 30 + "\n")

            sys.exit(1)

        # Filter out invalid files but continue with valid ones
        valid_vcf_list = validation_results.valid_files
        self.vcf_bad_list = (validation_results.corrupted_files +
                            validation_results.empty_files +
                            validation_results.permission_errors +
                            validation_results.unreadable_files)

        # Update counts after validation
        print(f"\n✅ PROCEEDING WITH {len(valid_vcf_list)} VALID VCF FILES")
        if self.vcf_bad_list:
            print(f"⚠️  EXCLUDED {len(self.vcf_bad_list)} PROBLEMATIC FILES")
        print("="*60 + "\n")

        # Process only valid VCF files
        if debug:
            print(f'Processing {len(valid_vcf_list)} validated VCF files')
            for vcf in valid_vcf_list:
                print(f'Processing: {os.path.basename(vcf)}')
                vcf, df, vcf_bad_list_temp = self.check_and_fix(vcf)
                try:
                    self.chrom = df['CHROM'].iloc[0]
                except (TypeError, AttributeError):
                    pass
                if df is not None:
                    dataframes[os.path.basename(vcf)] = df
                # Pre-validation catches most problems, but not a VCF that is well
                # formed and still unparseable here, so these must not be dropped.
                self.vcf_bad_list.extend(vcf_bad_list_temp)
        else:
            print(f'Processing: Pool processing with {cpu_count} cpus...')
            # Use context manager for process pool to ensure proper cleanup
            with futures.ProcessPoolExecutor(max_workers=cpu_count) as pool:
                for vcf, df, vcf_bad_list_temp in pool.map(self.check_and_fix, valid_vcf_list):
                    try:
                        self.chrom = df['CHROM'].iloc[0]
                    except (TypeError, AttributeError):
                        pass
                    if df is not None:
                        dataframes[os.path.basename(vcf)] = df
                    self.vcf_bad_list.extend(vcf_bad_list_temp)

        self.dataframes = dataframes
        print(f'\n\nDictionary of dataframes to memory runtime: {datetime.now() - self.startTime}\n')

        # Pickle for potential downstream use
        with open('dictionary_of_dataframes.pickle', 'wb') as handle:
            pickle.dump(dataframes, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if not debug:
            try:
                os.remove('dictionary_of_dataframes.pickle')
            except FileNotFoundError:
                pass

    def read_vcf(self, path):
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
        df['POS'] = pd.to_numeric(df['POS'], errors='coerce').fillna(0).astype(int)
        # QUAL and MQ stay nullable.  These were .fillna(0), which turned an ABSENT
        # quality value into the assertion "quality is zero".  For QUAL that means
        # the position falls below n_threshold and is rewritten to the REFERENCE
        # base - so a missing measurement became a positive claim that the sample
        # matches the reference.  A missing value is a no-call; Int64 carries that
        # distinction through to make_groupings, which now tests for it.
        # np.trunc keeps the truncation the previous .astype(int) performed on
        # fractional QUAL values such as 4523.5, so numeric behaviour is unchanged;
        # Int64 is what carries the missing value through as NA rather than 0.
        df['QUAL'] = np.trunc(pd.to_numeric(df['QUAL'], errors='coerce')).astype('Int64')

        # Split the INFO column and extract the AC, DP and MQ fields
        # Updated to handle malformed data more gracefully
        def extract_info_field(info_str, field):
            try:
                info_dict = dict(item.split("=") for item in info_str.split(";") if "=" in item)
                return info_dict.get(field, None)
            except (ValueError, AttributeError):
                return None
                
        df['AC'] = df['INFO'].apply(lambda x: extract_info_field(x, 'AC'))
        df['DP'] = df['INFO'].apply(lambda x: extract_info_field(x, 'DP'))
        df['MQ'] = df['INFO'].apply(lambda x: extract_info_field(x, 'MQ'))
        df['AC'] = pd.to_numeric(df['AC'], errors='coerce').fillna(0).astype(int)
        df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype(int)
        # An absent MQ is unknown, not zero.  As zero it silently failed the
        # MQ >= 56 filter, so those positions were dropped with no indication; and
        # if a whole file lacked MQ every position vanished and the sample went
        # missing from the analysis rather than being reported as unusable.
        df['MQ'] = np.trunc(pd.to_numeric(df['MQ'], errors='coerce')).astype('Int64')
        if df['MQ'].isna().all() and len(df):
            raise ValueError(
                'no MQ field in INFO on any record, so the mapping quality filter '
                'cannot be applied. Sample excluded rather than contributing zero '
                'positions silently. Re-call with a caller that emits MQ (or MQM).')
        df = df.drop(columns=['INFO', 'ID', 'FILTER', 'FORMAT'])

        return df

    def check_and_fix(self, vcf):
        '''
        Normalize a copy of vcf and parse it.  The returned path is always the
        ORIGINAL, so the dataframe key and therefore the sample name are unchanged.
        '''
        vcf_bad_list_temp = []
        df = None
        try:
            normalized = self.vcf_fix(vcf)
            df = self.read_vcf(normalized)
            df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
        except Exception as e:
            # Previously the input file was deleted here.  A parse failure is a
            # reason to exclude a sample and say why, not to destroy the evidence.
            vcf_bad_list_temp.append(
                f'{vcf}  could not be parsed: {type(e).__name__}: {e}')
        if df is not None:
            df = df.drop_duplicates(subset=['abs_pos'])
        return vcf, df, vcf_bad_list_temp

    def vcf_fix(self, vcf):
        '''
        Write a normalized copy of vcf and return its path.  The input is opened
        read-only; earlier versions renamed a temp file over it, which is why an
        os.utime call was needed to hide the fact that the user's files had been
        rewritten.

        Normalization is limited to making the file parseable: freebayes MQM is
        renamed to the VCF-standard MQ, and stray quoting that some tools emit is
        stripped.  It deliberately no longer invents data -- see below.
        '''
        normalized = os.path.join(self.normalized_dir, os.path.basename(vcf))
        try:
            return self._write_normalized(vcf, normalized)
        except Exception:
            # Do not leave a half-written copy behind: it would look like a valid
            # normalized VCF while containing only the records seen before the
            # problem, which is the failure mode this whole stage exists to remove.
            if os.path.exists(normalized):
                os.remove(normalized)
            raise

    def _write_normalized(self, vcf, normalized):
        synthesized_qual = False
        multi_sample = False

        with open(vcf, 'r') as file, open(normalized, 'w') as write_out:
            for line in file:
                line = line.replace('\r\n', '\n')
                if not line.rstrip():
                    continue
                line = line.rstrip()
                line = re.sub(r';MQM=', r';MQ=', line) #Allow Freebayes MQM to be read as MQ.  MQ is VCF standard
                line = re.sub(r'ID=MQM,', r'ID=MQ,', line)
                line = re.sub('"AC=', 'AC=', line)
                line = re.sub('""', '"', line)
                line = re.sub('""', '"', line)
                line = re.sub('""', '"', line)
                line = re.sub('"$', '', line)
                line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                line = re.sub('^"', '', line)
                if line.startswith('##') and line.endswith('"'):
                    line = re.sub('"$', '', line)
                if line.startswith('##'):
                    print(line.split('\t')[0], file=write_out)
                    continue

                line = re.sub('"', '', line)
                line = re.sub(r" +", "\t", line)
                fields = line.split('\t')

                if fields[0] == '#CHROM':
                    if len(fields) > 10:
                        multi_sample = True
                        if not self.first_sample_only:
                            raise ValueError(
                                f'multi-sample VCF ({len(fields) - 9} samples: '
                                f'{", ".join(fields[9:])}). vSNP3 expects one sample per '
                                f'file; split with `bcftools view -s <sample>`, or pass '
                                f'--first_sample_only to analyse {fields[9]} alone.')
                        print(f'  {os.path.basename(vcf)}: multi-sample VCF, keeping '
                              f'only {fields[9]}', flush=True)
                    print("\t".join(fields[0:10]), file=write_out)
                    continue

                if self.assume_gt_only_quality and len(fields) > 9 and fields[8] == 'GT':
                    # Records carrying only a GT have no quality data, so every
                    # threshold downstream rejects them.  On request substitute a
                    # QUAL and leave everything else alone.
                    #
                    # What this replaces: a regex that matched `<digits>\tGT\t<gt>`
                    # and rewrote it to `999\tGT:AD:DP:GQ:PL\t1/1:...`.  Because the
                    # digits it captured were the last numeric value in INFO rather
                    # than QUAL, it overwrote MQ (or DP, or AC) with 999 -- passing
                    # the MQ filter unconditionally -- and forced every genotype to
                    # 1/1.  A 0/1 heterozygous record became AC=999 and 1/1, which
                    # then failed both the AC==2 and AC==1 tests and vanished.
                    fields[5] = str(SYNTHESIZED_QUAL)
                    synthesized_qual = True

                print("\t".join(fields[0:10]), file=write_out)

        if synthesized_qual:
            self._insert_header_note(
                normalized,
                f'##vsnp3_synthesized_quality=QUAL={SYNTHESIZED_QUAL} substituted for '
                f'records whose FORMAT was GT only (--assume_gt_only_quality)')
        if multi_sample:
            self._insert_header_note(
                normalized, '##vsnp3_first_sample_only=other samples in this VCF were dropped')
        return normalized

    @staticmethod
    def _insert_header_note(vcf, note):
        '''Record a normalization decision in the VCF header so it is not invisible.'''
        with open(vcf) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if not line.startswith('##'):
                lines.insert(i, note + '\n')
                break
        with open(vcf, 'w') as f:
            f.writelines(lines)


class HTML_Summary():

    def __init__(self, runtime=None, vcf_to_df=None, reference=None, groupings_dict=None, raxml_version=None, all_vcf_boolen=None, args=None, removed_samples=None, validation_results=None, name_collisions=None, metadata_ambiguous=None):

        htmlfile = open(f'{global_working_dir}/vSNP_step2_summary-{global_date_stamp}.html', 'at', encoding='utf-8')
        
        #MAKE HTML FILE:
        print("<html>\n<head><meta charset=\"UTF-8\"><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)

        print(f"<h2>Script ran using <u>{reference} </u> variables:<br><br>", file=htmlfile)

        print('<div style="font-size:11px; font-weight:normal;">', file=htmlfile)
        if args.metadata:
            print(f"Metadata:  {redact_paths(str(args.metadata))}<br>", file=htmlfile)
        else:
            print("No metadata for describing samples in trees and tables<br>", file=htmlfile)
        if args.defining_snps:
            print(f"Defining SNPs:  {redact_paths(str(args.defining_snps))}<br>", file=htmlfile)
        else:
            print("No defining SNPs files for grouping and filtering<br>", file=htmlfile)
        if args.gbk:
            for each in args.gbk:
                print(f"gbk:  {redact_paths(str(each))}<br>", file=htmlfile)
        else:
            print("No gbk for annotation<br>", file=htmlfile)
        if args.remove_by_name:
            print(f"Remove from analysis:  {redact_paths(str(args.remove_by_name))}<br>", file=htmlfile)
        print('</div><br>', file=htmlfile)

        # Duplicate sample names.  The run completed on one file per name, so the
        # reader has to be told which file that was: the tables and tree carry the
        # name, and nothing else in the report would show that a second VCF claimed
        # it.  Placed high in the summary rather than in a footnote, because it
        # means a sample the user submitted is absent from the analysis.
        if name_collisions:
            print('<div style="font-size:11px; border:1px solid #8a5a00; '
                  'background:#fdf3e0; color:#8a5a00; padding:6px; margin:6px 0;">',
                  file=htmlfile)
            print(f'<b>Duplicate sample names: {len(name_collisions)}</b><br>',
                  file=htmlfile)
            print('More than one VCF resolved to the same sample name, so only one '
                  'was analysed under each. The reason is given per name below; '
                  'correct it and rerun to include the files left out.<br>',
                  file=htmlfile)
            for name, info in sorted(name_collisions.items()):
                dropped = ', '.join(redact_paths(str(d)) for d in info['dropped'])
                print(f'&nbsp;&nbsp;<b>{name}</b>: analysed '
                      f'{redact_paths(str(info["used"]))} &mdash; not analysed: '
                      f'{dropped}<br>', file=htmlfile)
                print(f'&nbsp;&nbsp;&nbsp;&nbsp;<i>because '
                      f'{redact_paths(str(info.get("cause", "")))}</i><br>',
                      file=htmlfile)
            print('</div>', file=htmlfile)

        # A file matching several worksheet rows keeps one name, so no count
        # changes and nothing else in the report would show it.  Separate box from
        # the one above: a different mistake, needing a different correction.
        if metadata_ambiguous:
            print('<div style="font-size:11px; border:1px solid #8a5a00; '
                  'background:#fdf3e0; color:#8a5a00; padding:6px; margin:6px 0;">',
                  file=htmlfile)
            print(f'<b>Ambiguous metadata rows: {len(metadata_ambiguous)}</b><br>',
                  file=htmlfile)
            print('These VCF files match more than one row of the metadata '
                  'worksheet, and those rows give different names. The first '
                  'matching row was used; remove the duplicate rows to make the '
                  'choice explicit.<br>', file=htmlfile)
            for base, info in sorted(metadata_ambiguous.items()):
                targets = ', '.join(str(t) for t in info['targets'])
                print(f'&nbsp;&nbsp;<b>{redact_paths(str(base))}</b>: '
                      f'"{info["key"]}" appears more than once, giving {targets} '
                      f'&mdash; used {info["targets"][0]}<br>', file=htmlfile)
            print('</div>', file=htmlfile)

        print(f'<span style="font-size:11px; font-weight:bold;">SNP calling thresholds:  REF: QUAL <u>&lt;{args.n_threshold}</u>, N: QUAL <u>{args.n_threshold}-{args.qual_threshold}</u>, ALT: QUAL <u>&gt;{args.qual_threshold}</u>, Ambigious: <u>AC=1</u>, MQ: <u>&gt;{args.mq_threshold}</u></span></h4>', file=htmlfile)

        # Add density filtering information to HTML summary
        if args.density_threshold is not None or args.density_window is not None:
            # Set defaults if not provided
            threshold = args.density_threshold if args.density_threshold is not None else 3
            window = args.density_window if args.density_window is not None else 20
            print(f"<h4>Density filtering enabled: SNPs removed when &gt;={threshold} SNPs found within {window} bp window</h4>", file=htmlfile)

        print(f"<h4>{vcf_to_df.vcf_original_count} VCF files initial count<br>", file=htmlfile)
        print(f"{len(vcf_to_df.dataframes)} VCF files in this run<br>", file=htmlfile)
        print(f"{len(vcf_to_df.vcf_bad_list)} VCF files in this run were corrupt and therefore removed</h4>", file=htmlfile)
        
        if all_vcf_boolen:
            print("\n<h4>All_VCFs is available</h4>", file=htmlfile)

        #TIME
        # Format runtime to show hours, minutes, seconds without decimals
        total_seconds = int(runtime.total_seconds())
        hours = total_seconds // 3600
        minutes = (total_seconds % 3600) // 60
        seconds = total_seconds % 60

        runtime_parts = []
        if hours > 0:
            runtime_parts.append(f"{hours} hours")
        if minutes > 0:
            runtime_parts.append(f"{minutes} minutes")
        runtime_parts.append(f"{seconds} seconds")

        runtime_formatted = " ".join(runtime_parts)
        print(f"Total run time: {runtime_formatted}: </h4>", file=htmlfile)

        # ENHANCED VCF VALIDATION RESULTS SECTION
        print("\n<h2>🔍 VCF File Validation Results</h2>", file=htmlfile)

        if validation_results:
            clean_count = validation_results.total_valid - len(validation_results.fixed_files)
            print("<table style='width: 80%;'>", file=htmlfile)
            print("<tr style='background-color: #f2f2f2;'><th>Validation Category</th><th>Count</th><th>Status</th></tr>", file=htmlfile)

            # Total files processed
            print(f"<tr><td><strong>Total VCF Files</strong></td><td>{validation_results.total_files}</td><td>📊 Processed</td></tr>", file=htmlfile)
            print(f"<tr><td><strong>Valid Files</strong></td><td>{validation_results.total_valid}</td><td style='color: green;'>✅ Passed</td></tr>", file=htmlfile)
            if validation_results.fixed_files:
                print(f"<tr><td>&nbsp;&nbsp;Clean (no issues)</td><td>{clean_count}</td><td style='color: green;'>✅ Clean</td></tr>", file=htmlfile)
                print(f"<tr><td>&nbsp;&nbsp;Auto-fixed (encoding corrected)</td><td>{len(validation_results.fixed_files)}</td><td style='color: darkorange;'>🔧 Fixed &amp; included</td></tr>", file=htmlfile)
            print(f"<tr><td><strong>Invalid Files</strong></td><td>{validation_results.total_invalid}</td><td style='color: red;'>❌ Failed</td></tr>", file=htmlfile)

            print("</table><br>", file=htmlfile)

            # Auto-fixed files section
            if validation_results.fixed_files:
                print(f"<h3 style='color: darkorange;'>🔧 Auto-Fixed Files ({len(validation_results.fixed_files)}) - Included in Analysis</h3>", file=htmlfile)
                print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                for entry in validation_results.fixed_files:
                    print(f"🔧 {os.path.basename(entry['file'])} &mdash; {entry['fix_applied']}<br>", file=htmlfile)
                print("</div><br>", file=htmlfile)

            # Detailed validation issues
            if validation_results.total_invalid > 0:
                print("<h3>📋 Detailed Validation Issues</h3>", file=htmlfile)

                # Corrupted files
                if validation_results.corrupted_files:
                    print(f"<h4 style='color: red;'>💥 Corrupted Files ({len(validation_results.corrupted_files)})</h4>", file=htmlfile)
                    print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                    for corrupted_file in validation_results.corrupted_files:
                        print(f"❌ {os.path.basename(corrupted_file)}<br>", file=htmlfile)
                    print("</div><br>", file=htmlfile)

                # Reference mismatches (should be empty if we got this far)
                if validation_results.reference_mismatches:
                    print(f"<h4 style='color: red;'>🔥 Reference Mismatches ({len(validation_results.reference_mismatches)})</h4>", file=htmlfile)
                    print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                    for mismatch in validation_results.reference_mismatches:
                        print(f"❌ {os.path.basename(mismatch['file'])}: {mismatch['message']}<br>", file=htmlfile)
                    print("</div><br>", file=htmlfile)

                # Empty files
                if validation_results.empty_files:
                    print(f"<h4 style='color: orange;'>⚠️ Empty Files ({len(validation_results.empty_files)})</h4>", file=htmlfile)
                    print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                    for empty_file in validation_results.empty_files:
                        print(f"⚠️ {os.path.basename(empty_file)}<br>", file=htmlfile)
                    print("</div><br>", file=htmlfile)

                # Permission errors
                if validation_results.permission_errors:
                    print(f"<h4 style='color: orange;'>🔒 Permission Issues ({len(validation_results.permission_errors)})</h4>", file=htmlfile)
                    print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                    for perm_file in validation_results.permission_errors:
                        print(f"🔒 {os.path.basename(perm_file)}<br>", file=htmlfile)
                    print("</div><br>", file=htmlfile)

                # Unreadable files
                if validation_results.unreadable_files:
                    print(f"<h4 style='color: red;'>❓ Unreadable Files ({len(validation_results.unreadable_files)})</h4>", file=htmlfile)
                    print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                    for unreadable_file in validation_results.unreadable_files:
                        print(f"❓ {os.path.basename(unreadable_file)}<br>", file=htmlfile)
                    print("</div><br>", file=htmlfile)

            elif not validation_results.fixed_files:
                print("<h3 style='color: green;'>✅ All VCF files passed validation!</h3>", file=htmlfile)

        else:
            # Fallback to original error reporting if validation_results not available
            if len(vcf_to_df.vcf_bad_list) < 1:
                print("<h3 style='color: green;'>No corrupt files found</h3>", file=htmlfile)
            else:
                print(f"\n<h3 style='color: red;'>Corrupt files removed ({len(vcf_to_df.vcf_bad_list)})</h3>", file=htmlfile)
                print("<div style='margin-left: 20px; font-size: 10px;'>", file=htmlfile)
                for i in vcf_to_df.vcf_bad_list:
                    print(f"❌ {os.path.basename(i)}<br>", file=htmlfile)
                print("</div><br>", file=htmlfile)

        #GROUPING TABLE
        group_vcfs_dict = defaultdict(list) #invert the key, values
        for group, dataframes in groupings_dict.items():
            for vcf in dataframes.keys():
                group_vcfs_dict[vcf].append(group)
        group_vcfs_dict = dict(sorted(group_vcfs_dict.items())) #sorts on key(vcf name)

        print(f'<h2>Groupings with {len(group_vcfs_dict):,} listed:</h2>', file=htmlfile)
        print("<table>", file=htmlfile)
        print("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)

        for key, value in group_vcfs_dict.items():
            print("<tr>", file=htmlfile)
            print(f"<td>{key}</td>", end='\t', file=htmlfile)
            for group in value:
                if group == "Group Not Found":
                    print(f'<td><span style="color: red;">{group}</span></td>', end='\t', file=htmlfile)
                else:
                    print(f"<td>{group}</td>", end='\t', file=htmlfile)
            print("</tr>", file=htmlfile)
        print("</table>", file=htmlfile)

        # Removed from analysis
        if removed_samples:
            if len(removed_samples) < 1:
                print("<h2>No samples purposely removed from initial VCF file dataset</h2>", file=htmlfile)
            else:
                print('\n<h2>VCF files removed from dataset using "remove_from_analysis.xlsx".</h2>', file=htmlfile)
                for each in removed_samples:
                    print(f"{os.path.basename(each)} <br>", file=htmlfile)
                print("<br>", file=htmlfile)

        try:
            import platform
            print("\n<h2>System Information:</h2>", file=htmlfile)
            
            # Get OS information
            os_name = platform.system()
            os_version = platform.version()
            os_release = platform.release()
            
            # Get architecture information
            arch = platform.machine()
            processor = platform.processor()
            
            # Print OS information with specific details based on OS type
            print(f"<b>Operating System:</b> {os_name} {os_release} {os_version}<br>", file=htmlfile)
            
            # Get and print detailed OS information based on the OS type
            if os_name == 'Darwin':  # macOS
                # Check if ARM (Apple Silicon) or Intel
                if arch == 'arm64':
                    cpu_type = "ARM (Apple Silicon)"
                else:
                    cpu_type = "Intel"
                
                # Get macOS version name
                mac_ver = platform.mac_ver()
                macos_version = f"macOS {mac_ver[0]}"
                
                print(f"<b>macOS Details:</b> {macos_version}, {cpu_type}<br>", file=htmlfile)
                
            elif os_name == 'Linux':
                # Try to get Linux distribution info
                try:
                    import distro
                    linux_distro = distro.name(pretty=True)
                except ImportError:
                    # Fallback if distro module is not available
                    try:
                        with open('/etc/os-release') as f:
                            lines = f.readlines()
                            for line in lines:
                                if line.startswith('PRETTY_NAME='):
                                    linux_distro = line.split('=')[1].strip().strip('"')
                                    break
                            else:
                                linux_distro = "Unknown Linux Distribution"
                    except:
                        linux_distro = "Unknown Linux Distribution"
                
                # Check for HPC environment
                is_hpc = False
                hpc_info = "Unknown"
                
                # Check for common HPC environment variables or files
                hpc_indicators = {
                    'SLURM_CLUSTER_NAME': 'SLURM',
                    'PBS_HOME': 'PBS',
                    'SGE_ROOT': 'SGE',
                    'LSB_JOBID': 'LSF'
                }
                
                for env_var, hpc_type in hpc_indicators.items():
                    if env_var in os.environ:
                        is_hpc = True
                        hpc_info = f"{hpc_type} HPC environment"
                        break
                
                # If no environment variables found, check for common HPC directories
                if not is_hpc:
                    hpc_paths = [
                        ('/opt/slurm', 'SLURM'),
                        ('/opt/pbs', 'PBS'),
                        ('/opt/sge', 'SGE'),
                        ('/opt/lsf', 'LSF')
                    ]
                    
                    for path, hpc_type in hpc_paths:
                        if os.path.exists(path):
                            is_hpc = True
                            hpc_info = f"{hpc_type} HPC environment"
                            break
                
                if is_hpc:
                    print(f"<b>Linux Details:</b> {linux_distro}, {hpc_info}<br>", file=htmlfile)
                else:
                    print(f"<b>Linux Details:</b> {linux_distro}<br>", file=htmlfile)
                
            elif os_name == 'Windows':
                win_ver = platform.win32_ver()
                win_edition = win_ver[0]
                win_build = win_ver[1]
                print(f"<b>Windows Details:</b> Windows {win_edition} (Build {win_build})<br>", file=htmlfile)
            
            # Print CPU architecture information
            print(f"<b>CPU Architecture:</b> {arch}<br>", file=htmlfile)
            print(f"<b>Processor:</b> {processor}<br>", file=htmlfile)
            
            # Get more detailed CPU information using py-cpuinfo if available
            try:
                import cpuinfo
                cpu_info = cpuinfo.get_cpu_info()
                print(f"<b>CPU Model:</b> {cpu_info['brand_raw']}<br>", file=htmlfile)
                print(f"<b>CPU Cores:</b> {cpu_info['count']}<br>", file=htmlfile)
            except (ImportError, Exception) as e:
                # Fall back to less detailed information
                pass
            
            print("<hr>", file=htmlfile)
            print("\n<h2>Program versions:</h2>", file=htmlfile)
            print(f'vSNP3: {__version__} <br>', file=htmlfile)
            print(f'Python: {sys.version} <br>', file=htmlfile)
            
            # Define the list of programs to check
            program_list = [
                'biopython', 'dask', 'humanize', 'numpy', 'pandas', 'openpyxl', 
                'xlsxwriter', 'parallel', 'pigz', 'regex', 'py-cpuinfo', 'raxml', 'plotly'
            ]
            
            # Dictionary for module name mapping (conda/pip name to import name)
            module_mapping = {
                'python': 'Python',
                'biopython': 'Bio',
                'numpy': 'numpy',
                'pandas': 'pandas',
                'openpyxl': 'openpyxl',
                'xlsxwriter': 'xlsxwriter',
                'regex': 're',
                'py-cpuinfo': 'cpuinfo',
                'plotly': 'plotly',
                'dask': 'dask',
                'humanize': 'humanize',
            }
            
            # Check if running in a conda environment
            conda_env = os.environ.get('CONDA_DEFAULT_ENV')
            if conda_env:
                print(f'<b>Versions from conda environment: {conda_env}</b> <br>', file=htmlfile)
            else:
                print('<b>Versions from system installation</b> <br>', file=htmlfile)
            
            for program in program_list:
                version = "nd"  # Default to "nd" (no data)
                
                # First try to get the version from conda
                try:
                    conda_output = subprocess.check_output(["conda", "list", program], 
                                                        stderr=subprocess.STDOUT, 
                                                        universal_newlines=True)
                    # Extract version from conda output
                    lines = conda_output.strip().split('\n')
                    for line in lines:
                        if program in line and not line.startswith('#'):
                            parts = line.split()
                            if len(parts) >= 2:
                                version = parts[1]  # Version is typically the second column
                                break
                except (subprocess.CalledProcessError, FileNotFoundError):
                    # If conda check fails, try to check via module if applicable
                    if program in module_mapping:
                        module_name = module_mapping[program]
                        try:
                            if module_name in sys.modules and hasattr(sys.modules[module_name], '__version__'):
                                version = sys.modules[module_name].__version__
                            elif module_name not in sys.modules:
                                # Try to import the module
                                module = __import__(module_name)
                                if hasattr(module, '__version__'):
                                    version = module.__version__
                        except (ImportError, AttributeError):
                            pass
                    
                    # If still no version, try system command
                    if version == "nd" and shutil.which(program):
                        try:
                            # Different programs report versions differently
                            if program == 'bcftools' or program == 'samtools' or program == 'bwa':
                                cmd_output = subprocess.check_output([program, "--version"], 
                                                                stderr=subprocess.STDOUT, 
                                                                universal_newlines=True)
                                version = cmd_output.strip().split('\n')[0].split(' ')[1]
                            elif program == 'raxml':
                                cmd_output = subprocess.check_output([program, "-v"], 
                                                                stderr=subprocess.STDOUT, 
                                                                universal_newlines=True)
                                version = cmd_output.strip()
                                # Process raxml version string as in the original code
                                version = re.sub('This is ', '', version)
                                version = re.sub(' released by.*', '', version)
                            elif program == 'minimap2' or program == 'spades.py' or program == 'seqkit':
                                cmd_output = subprocess.check_output([program, "--version"], 
                                                                stderr=subprocess.STDOUT, 
                                                                universal_newlines=True)
                                version = cmd_output.strip()
                            elif program == 'freebayes':
                                cmd_output = subprocess.check_output([program, "--version"], 
                                                                stderr=subprocess.STDOUT, 
                                                                universal_newlines=True)
                                version = cmd_output.strip().split(',')[0].split(' ')[-1]
                            elif program == 'pigz' or program == 'parallel':
                                cmd_output = subprocess.check_output([program, "--version"], 
                                                                stderr=subprocess.STDOUT, 
                                                                universal_newlines=True)
                                version = cmd_output.strip().split('\n')[0]
                        except (subprocess.CalledProcessError, FileNotFoundError):
                            pass
                
                # Print the version information
                source = ""
                if version != "nd":
                    # Determine if the version came from conda or system
                    try:
                        conda_output = subprocess.check_output(["conda", "list", program], 
                                                            stderr=subprocess.STDOUT, 
                                                            universal_newlines=True)
                        if program in conda_output:
                            source = "(conda)"
                        else:
                            source = "(system)"
                    except:
                        source = "(system)"
                
                print(f'{program}: {version} {source} <br>', file=htmlfile)
            
        except Exception as e:
            print(f"Error checking versions: {str(e)} <br>", file=htmlfile)

        print("</body>\n</html>", file=htmlfile)
        htmlfile.close()


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
    parser.add_argument('-w', '--qual_threshold', action='store', dest='qual_threshold', default=150, required=False, help='Optional: Minimum QUAL threshold for calling a SNP')
    parser.add_argument('-x', '--n_threshold', action='store', dest='n_threshold', default=50, required=False, help='Optional: Minimum N threshold.  SNPs between this and qual_threshold are reported as N')
    parser.add_argument('-y', '--mq_threshold', action='store', dest='mq_threshold', default=56, required=False, help='Optional: At least one position per group must have this minimum MQ threshold to be called.')
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
    parser.add_argument('--density_threshold', nargs='?', const=3, type=int, dest='density_threshold', help='Optional: Minimum number of SNPs required to trigger density filtering (default: 3)')
    parser.add_argument('--density_window', nargs='?', const=20, type=int, dest='density_window', help='Optional: Window size in base pairs for density filtering (default: 20)')
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
        print(f"   Pattern: *.vcf")
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

    def zipit(src, dst):
        zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
        abs_src = os.path.abspath(src)
        for dirname, subdirs, files in os.walk(src):
            for filename in files:
                absname = os.path.abspath(os.path.join(dirname, filename))
                arcname = absname[len(abs_src) + 1:]
                zf.write(absname, arcname)
        zf.close()
        try:
            shutil.rmtree(src)
        except (FileNotFoundError, PermissionError) as e:
            print(f"Warning: Could not remove directory {src}: {e}")

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

    starting_files = os.path.join(setup.cwd, 'vcf_starting_files')
    os.makedirs(starting_files, exist_ok=True)
    for each_vcf in wd_vcf_list:
        shutil.copy(each_vcf, starting_files)

    vcf_to_df = VCF_to_DF(vcf_list=wd_vcf_list, debug=args.debug,
                          assume_gt_only_quality=args.assume_gt_only_quality,
                          first_sample_only=args.first_sample_only) #write_out=args.write_out
    if args.fix_vcfs:
        print(f'\nNormalized copies written to {vcf_to_df.normalized_dir}')
        print('The input VCF files were not modified.')
        if vcf_to_df.vcf_bad_list:
            print(f'\n{len(vcf_to_df.vcf_bad_list)} file(s) could not be normalized:')
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

    # The VCF files that were read are archived to vcf_starting_files.zip, so
    # leaving the loose copies behind just clutters the working directory --
    # deleting them is the useful default.  What must never happen is deleting
    # someone else's files: wd_vcf_list is the caller's own list unless
    # --output_dir made copies first, which is how `vsnp3_step2.py -wd
    # /path/to/vcf_database` used to empty a shared database.  -k/--keep_ind_vcfs
    # existed to prevent that and was parsed but never read.
    #
    # So: delete by default, but only files this run owns.
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

    print('\nvcf_bad_list')
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
        
    if not args.remove_by_name and ro.remove:
        args.remove_by_name = ro.remove
        print(f"Using reference type remove list: {args.remove_by_name}")
    elif args.remove_by_name:
        print(f"Using explicitly provided remove list: {args.remove_by_name}")
    
    print(f'Before sample filter: {len(vcf_to_df.dataframes)}')
    if args.remove_by_name:
        remove_from_analysis = Remove_From_Analysis(working_directory=global_working_dir, excel_remove=args.remove_by_name, extension="vcf")
        actually_removed = []
        for each in remove_from_analysis.remove_list:
            sample_basename = os.path.basename(each)
            # Only track samples that were actually present in the dataframes
            if sample_basename in vcf_to_df.dataframes:
                vcf_to_df.dataframes.pop(sample_basename, None)
                actually_removed.append(each)
        remove_list = actually_removed  # Only samples that were actually found and removed
    else:
        remove_list = None
    print(f'After sample filter: {len(vcf_to_df.dataframes)}')

    if args.defining_snps:
        shutil.copy(args.defining_snps, starting_files) #package with starting files for the record
    zipit(starting_files, starting_files) # zip starting files directory

    # Before creating the Group instance, check if html_tree is True and set dp accordingly
    if args.html_tree:
        args.dp = True

    group = Group(cwd=global_working_dir, metadata=args.metadata, defining_snps=args.defining_snps, excel_remove=args.remove_by_name, gbk_list=args.gbk, dataframes=vcf_to_df.dataframes, all_vcf=args.all_vcf, find_new_filters=args.find_new_filters, no_filters=args.no_filters, qual_threshold=int(args.qual_threshold), n_threshold=int(args.n_threshold), mq_threshold=int(args.mq_threshold), abs_pos=args.abs_pos, group=args.group, show_groups=args.show_groups, hash_groups=args.hash_groups, html_tree=args.html_tree, dp=args.dp, filter_density=filter_density, density_threshold=args.density_threshold, density_window=args.density_window, debug=args.debug)
    vcf_to_df.vcf_bad_list = vcf_to_df.vcf_bad_list + group.vcf_bad_list

    setup.print_time()
    HTML_Summary(runtime=setup.run_time, vcf_to_df=vcf_to_df, reference=ro.select_ref, groupings_dict=group.groupings_dict, raxml_version=group.raxml_version, all_vcf_boolen=args.all_vcf, args=args, removed_samples=remove_list, validation_results=vcf_to_df.validation_results, name_collisions=group.name_collisions, metadata_ambiguous=group.metadata_ambiguous) 
    
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