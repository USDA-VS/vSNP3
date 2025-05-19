#!/usr/bin/env python

__version__ = "3.29"

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

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

global_date_stamp=None
global_working_dir='.'

class VCF_to_DF():
    ''' 
    '''

    def __init__(self, vcf_list=None, debug=False): #write_out=False, 
        '''
        Start at class call
        '''
        self.startTime = datetime.now()
        self.vcf_bad_list=[]
        cpu_count = int(multiprocessing.cpu_count() / 1.2)
        dataframes={}
        self.vcf_original_count = len(vcf_list)
        if debug:
            print(f'VCF file count {self.vcf_original_count}')
            for vcf in vcf_list:
                print(vcf)
                vcf, df, vcf_bad_list_temp = self.check_and_fix(vcf)
                try:
                    self.chrom = df['CHROM'].iloc[0]
                except (TypeError, AttributeError):
                    pass
                if df is not None:
                    dataframes[os.path.basename(vcf)] = df
                self.vcf_bad_list = self.vcf_bad_list + vcf_bad_list_temp
        else:
            print(f'Fixing: Pool processing with {cpu_count} cpus...')
            # Use context manager for process pool to ensure proper cleanup
            with futures.ProcessPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## process works best for calling on multiple files
                for vcf, df, vcf_bad_list_temp in pool.map(self.check_and_fix, vcf_list):
                    try:
                        self.chrom = df['CHROM'].iloc[0]
                    except (TypeError, AttributeError):
                        pass
                    if df is not None:
                        dataframes[os.path.basename(vcf)] = df
                    self.vcf_bad_list = self.vcf_bad_list + vcf_bad_list_temp
        self.dataframes = dataframes
        print(f'\n\nDictionary of dataframes to memory runtime: {datetime.now() - self.startTime}\n')
        # if write_out: # write out pickle file can be used in downstream applications
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

    def __init__(self, runtime=None, vcf_to_df=None, reference=None, groupings_dict=None, raxml_version=None, all_vcf_boolen=None, args=None, removed_samples=None):

        htmlfile = open(f'{global_working_dir}/vSNP_step2_summary-{global_date_stamp}.html', 'at')
        
        #MAKE HTML FILE:
        print("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)

        print(f"<h2>Script ran using <u>{reference} </u> variables<br>", file=htmlfile)

        if args.metadata:
            print(f"<h4>Metadata:  {args.metadata}<br>", file=htmlfile)
        else:
            print(f"No metadata for describing samples in trees and tables<br>", file=htmlfile)
        if args.defining_snps:
            print(f"Defining SNPs:  {args.defining_snps}<br>", file=htmlfile)
        else:
            print(f"No defining SNPs files for grouping and filtering<br>", file=htmlfile)
        if args.gbk:
            for each in args.gbk:
                print(f"gbk:  {each}<br>", file=htmlfile)
        else:
            print(f"No gbk for annotation<br>", file=htmlfile)

        print(f"SNP calling thresholds:  REF: QUAL <u><{args.n_threshold}</u>, N: QUAL <u>{args.n_threshold}-{args.qual_threshold}</u>, ALT: QUAL <u>>{args.qual_threshold}</u>, Ambigious: <u>AC=1</u>, MQ: <u>>{args.mq_threshold}</u></h4>", file=htmlfile)

        print(f"<h4>{vcf_to_df.vcf_original_count} VCF files initial count<br>", file=htmlfile)
        print(f"{len(vcf_to_df.dataframes)} VCF files in this run<br>", file=htmlfile)
        print(f"{len(vcf_to_df.vcf_bad_list)} VCF files in this run were corrupt and therefore removed</h4>", file=htmlfile)
        
        if all_vcf_boolen:
            print("\n<h4>All_VCFs is available</h4>", file=htmlfile)

        #TIME
        print(f"Total run time: {runtime}: </h4>", file=htmlfile)

        # ERROR LIST
        if len(vcf_to_df.vcf_bad_list) < 1:
            print("<h2>No corrupt files found</h2>", file=htmlfile)
        else:
            print("\n<h2>Corrupt files removed</h2>", file=htmlfile)
            for i in vcf_to_df.vcf_bad_list:
                print(f"{i} <br>", file=htmlfile)
            print("<br>", file=htmlfile)

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
                print('\n<h2>VCF files in "remove_from_analysis.xlsx" and not included from dataset.</h2>', file=htmlfile)
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
                print(f'<b>Versions from system installation</b> <br>', file=htmlfile)
            
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
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling.  A pickle file will be kept for troubleshooting to be used directly in vsnp3_group_on_defining_snps.py.  This saves processing time')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    setup = Setup(debug=args.debug)
    global_date_stamp = setup.date_stamp
    global_working_dir = setup.cwd
    cwd_test = False
    if args.wd == '.':
        cwd_test = True
        vcf_list = glob.glob('*vcf')
        wd_vcf_list = vcf_list
    else:
        #get VCFs from a directory
        wd = os.path.expanduser(args.wd)
        wd = os.path.abspath(wd)
        vcf_list = glob.glob(f'{wd}/*vcf')
        wd_vcf_list = vcf_list

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
            wd_vcf_list.append(f'{output_dir}/{os.path.basename(each_vcf)}')
        os.chdir(output_dir)
        setup.cwd = os.getcwd()
        global_working_dir = setup.cwd

    starting_files = f'{setup.cwd}/vcf_starting_files'
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

    print(f'\nvcf_bad_list')
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
        for each in remove_from_analysis.remove_list:
            vcf_to_df.dataframes.pop(os.path.basename(each), None)
        remove_list = remove_from_analysis.remove_list
    else:
        remove_list = None
    print(f'After sample filter: {len(vcf_to_df.dataframes)}')

    if args.defining_snps:
        shutil.copy(args.defining_snps, starting_files) #package with starting files for the record
    zipit(starting_files, starting_files) # zip starting files directory

    # Before creating the Group instance, check if html_tree is True and set dp accordingly
    if args.html_tree:
        args.dp = True

    group = Group(cwd=global_working_dir, metadata=args.metadata, defining_snps=args.defining_snps, excel_remove=args.remove_by_name, gbk_list=args.gbk, dataframes=vcf_to_df.dataframes, all_vcf=args.all_vcf, find_new_filters=args.find_new_filters, no_filters=args.no_filters, qual_threshold=int(args.qual_threshold), n_threshold=int(args.n_threshold), mq_threshold=int(args.mq_threshold), abs_pos=args.abs_pos, group=args.group, show_groups=args.show_groups, hash_groups=args.hash_groups, html_tree=args.html_tree, dp=args.dp, debug=args.debug)
    vcf_to_df.vcf_bad_list = vcf_to_df.vcf_bad_list + group.vcf_bad_list

    setup.print_time()
    HTML_Summary(runtime=setup.run_time, vcf_to_df=vcf_to_df, reference=ro.select_ref, groupings_dict=group.groupings_dict, raxml_version=group.raxml_version, all_vcf_boolen=args.all_vcf, args=args, removed_samples=remove_list) 
    
    try:
        os.remove(notification_file)
        print(f"Deleted file: {notification_file}")
    except FileNotFoundError:
        pass
# Created 2021 by Tod Stuber