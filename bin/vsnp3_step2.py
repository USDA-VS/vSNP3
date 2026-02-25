#!/usr/bin/env python

__version__ = "3.34"

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
import zipfile
import glob
from datetime import datetime

from collections import defaultdict
from concurrent import futures
import multiprocessing

# Move set_start_method inside if __name__ == "__main__" to avoid issues with Python 3.12
# multiprocessing.set_start_method('spawn', True)

import warnings
warnings.filterwarnings('ignore')

from vsnp3_file_setup import Setup
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

    def __init__(self, vcf_list=None, debug=False): #write_out=False,
        '''
        Start at class call with comprehensive VCF validation
        '''
        self.startTime = datetime.now()
        self.vcf_bad_list = []
        self.vcf_original_count = len(vcf_list)
        cpu_count = int(multiprocessing.cpu_count() / 1.2)
        dataframes = {}

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
                # Note: vcf_bad_list_temp is now redundant due to pre-validation
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
                    # Note: vcf_bad_list_temp should be minimal due to pre-validation

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
        df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce').fillna(0).astype(int)
        
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
        df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce').fillna(0).astype(int)
        df = df.drop(columns=['INFO', 'ID', 'FILTER', 'FORMAT'])

        return df

    def check_and_fix(self, vcf):
        vcf_bad_list_temp = []
        try:
            self.vcf_fix(vcf)
            df = self.read_vcf(vcf)
            df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
        except RuntimeError:
            self.vcf_fix(vcf)
            try:
                df = self.read_vcf(vcf)
                df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
            except Exception as e:
                vcf_bad_list_temp.append(vcf)
                try:
                    os.remove(vcf)
                except FileNotFoundError:
                    pass
                df = None
        try:
            if df is not None:
                df = df.drop_duplicates(subset=['abs_pos'])
        except AttributeError:
            # pass if df is empty, NoneType
            pass
        return vcf, df, vcf_bad_list_temp

    def vcf_fix(self, vcf):
        temp_file = vcf + ".temp"
        write_out = open(temp_file, 'w') #r+ used for reading and writing to the same file
        initial_file_time_stats = os.stat(vcf)
        with open(vcf, 'r') as file:
            for line in file:
                line = line.replace('\r\n', '\n')
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip()  #remove right white space
                    line = re.sub(r';MQM=', r';MQ=', line) #Allow Freebayes MQM to be read as MQ.  MQ is VCF standard
                    line = re.sub(r'ID=MQM,', r'ID=MQ,', line)
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub(r'[0-9]+\tGT\t.\/.+', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##') and line.endswith('"'):
                        line = re.sub('"$', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = re.sub(r" +", "\t", line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        write_out.close()  # Explicitly close the file before renaming
        os.rename(temp_file, vcf)
        os.utime(vcf, times=(initial_file_time_stats.st_mtime, initial_file_time_stats.st_mtime))


class HTML_Summary():

    def __init__(self, runtime=None, vcf_to_df=None, reference=None, groupings_dict=None, raxml_version=None, all_vcf_boolen=None, args=None, removed_samples=None, validation_results=None):

        htmlfile = open(f'{global_working_dir}/vSNP_step2_summary-{global_date_stamp}.html', 'at', encoding='utf-8')
        
        #MAKE HTML FILE:
        print("<html>\n<head><meta charset=\"UTF-8\"><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)

        print(f"<h2>Script ran using <u>{reference} </u> variables:<br><br>", file=htmlfile)

        print('<div style="font-size:11px; font-weight:normal;">', file=htmlfile)
        if args.metadata:
            print(f"Metadata:  {args.metadata}<br>", file=htmlfile)
        else:
            print("No metadata for describing samples in trees and tables<br>", file=htmlfile)
        if args.defining_snps:
            print(f"Defining SNPs:  {args.defining_snps}<br>", file=htmlfile)
        else:
            print("No defining SNPs files for grouping and filtering<br>", file=htmlfile)
        if args.gbk:
            for each in args.gbk:
                print(f"gbk:  {each}<br>", file=htmlfile)
        else:
            print("No gbk for annotation<br>", file=htmlfile)
        if args.remove_by_name:
            print(f"Remove from analysis:  {args.remove_by_name}<br>", file=htmlfile)
        print('</div><br>', file=htmlfile)

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
                'cairosvg': 'cairosvg',
                'dask': 'dask',
                'humanize': 'humanize',
                'svgwrite': 'svgwrite'
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
    parser.add_argument('-f', '--fix_vcfs', action='store_true', dest='fix_vcfs', help='Optional: Just fix VCF files and exit')
    parser.add_argument('-k', '--keep_ind_vcfs', action='store_true', dest='keep_ind_vcfs', default=False, help='Optional: Keep VCF files in current working directory when VCF files in current working director are used, VCF files are always saved and zipped in "vcf_starting_files.zip".')
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
        vcf_list = glob.glob('*vcf')
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
        vcf_list = glob.glob(vcf_pattern)
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

    vcf_to_df = VCF_to_DF(vcf_list=wd_vcf_list, debug=args.debug) #write_out=args.write_out
    if args.fix_vcfs:
        sys.exit(0)

    # Create the file to indicate the script is running
    notification_file = "step2_is_running__individual_folders_may_be_complete"   
    with open(notification_file, 'w') as f:
        f.write("Script is still running.")
    print(f"Created file: {notification_file}")

    #rm move vcfs from working directory
    for each_vcf in wd_vcf_list:
        try:
            os.remove(each_vcf)
        except FileNotFoundError:
            # if file was previously removed such as it was empty
            pass

    print('\nvcf_bad_list')
    for each in vcf_to_df.vcf_bad_list:
        print(f'\t {each}')

    if args.reference_type:
        ro = Ref_Options(args.reference_type)
    else:
        ro = Ref_Options(vcf_to_df.chrom)

    if args.abs_pos and not args.group:
        print('\n### -abs_pos must be used with -group option\n')
        sys.exit()
    if args.group and not args.abs_pos:
        print('\n### -group must be used with -abs_pos option\n')
        sys.exit()

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
    HTML_Summary(runtime=setup.run_time, vcf_to_df=vcf_to_df, reference=ro.select_ref, groupings_dict=group.groupings_dict, raxml_version=group.raxml_version, all_vcf_boolen=args.all_vcf, args=args, removed_samples=remove_list, validation_results=vcf_to_df.validation_results) 
    
    try:
        os.remove(notification_file)
        print(f"Deleted file: {notification_file}")
    except FileNotFoundError:
        pass
# Created 2021 by Tod Stuber