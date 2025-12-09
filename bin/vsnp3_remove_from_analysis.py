#!/usr/bin/env python

__version__ = "3.32"

import os
import sys
import pandas as pd
import glob
import argparse
import textwrap
from pathlib import Path


class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'


class Remove_From_Analysis:
    '''
    Class to handle removal of samples from analysis based on an Excel file list
    ''' 
    def __init__(self, working_directory='.', excel_remove=None, extension=None):
        if working_directory == '.':
            working_directory = os.getcwd()
        
        # Use Path for better path handling
        working_dir_path = Path(working_directory)
        if working_dir_path.is_absolute():
            print(f'working directory: {working_directory}')
        else:
            print('##### PROVIDE A FULL PATH')
            print(f'directory given: "{working_directory}"')
            sys.exit(0)

        # Validate Excel file exists
        if not os.path.isfile(excel_remove):
            print(f'{bcolors.RED}ERROR: Excel file not found: {excel_remove}{bcolors.ENDC}')
            sys.exit(1)

        # Updated pandas read_excel approach with error handling
        try:
            # Read the first column with no header
            df = pd.read_excel(excel_remove, index_col=0, usecols=[0], header=None)
        except Exception as e:
            print(f'{bcolors.RED}ERROR: Failed to read Excel file {excel_remove}: {str(e)}{bcolors.ENDC}')
            sys.exit(1)
        
        if df.empty:
            print(f'{bcolors.YELLOW}WARNING: Excel file is empty or has no valid data{bcolors.ENDC}')
            self.remove_list = []
        else:
            remove_list = []
            for each_sample in df.index:
                # Convert to string in case of numeric values
                sample_name = str(each_sample).strip()
                if not sample_name:  # Skip empty sample names
                    continue
                
                # Use Path.joinpath for robust cross-platform path handling
                base_path = working_dir_path / sample_name
                base_with_ext = working_dir_path / f"{sample_name}.{extension}"
                base_with_zc_ext = working_dir_path / f"{sample_name}_zc.{extension}"
                
                # Convert to string for glob compatibility
                remove_list.append(str(base_path))  # if .vcf is supplied with the sample name
                remove_list.append(str(base_with_ext))  # most common behavior
                remove_list.append(str(base_with_zc_ext))  # allow _zc to not be specified
            
            self.remove_list = remove_list
        
        self.excel_remove = excel_remove
        self.working_directory = working_directory

    def remove_files(self):
        num_files_removed = 0
        files_not_found = []
        files_failed_to_remove = []
        
        print(f'Removing samples listed in {self.excel_remove}')
        
        for each_sample in self.remove_list:
            try:
                # Escape special characters in the path for glob
                escaped_sample = each_sample.replace('[', r'\[').replace(']', r'\]')
                glob_list = glob.glob(escaped_sample)
                
                if not glob_list:
                    files_not_found.append(each_sample)
                    continue
                
                for item in glob_list:
                    try:
                        if os.path.isfile(item):
                            os.remove(item)
                            num_files_removed += 1
                            # print(f'\tRemoving: {item}')
                        else:
                            print(f'{bcolors.YELLOW}WARNING: Not a file, skipping: {item}{bcolors.ENDC}')
                    except OSError as e:
                        files_failed_to_remove.append(f"{item}: {str(e)}")
                        print(f'{bcolors.RED}ERROR: Failed to remove {item}: {str(e)}{bcolors.ENDC}')
                        
            except Exception as e:
                print(f'{bcolors.RED}ERROR: Failed to process pattern {each_sample}: {str(e)}{bcolors.ENDC}')
        
        # Summary output
        print(f'\n{bcolors.RED}{num_files_removed:,} {bcolors.ENDC}{bcolors.WHITE}files removed from the analysis{bcolors.ENDC}\n')
        
        if files_not_found:
            print(f'{bcolors.YELLOW}Note: {len(files_not_found)} file patterns had no matches{bcolors.ENDC}')
        
        if files_failed_to_remove:
            print(f'{bcolors.RED}WARNING: {len(files_failed_to_remove)} files could not be removed{bcolors.ENDC}')
            for failed_file in files_failed_to_remove:
                print(f'  {failed_file}')
        
        self.removed_file_count = num_files_removed
        self.files_not_found = files_not_found
        self.files_failed_to_remove = files_failed_to_remove


if __name__ == "__main__":  # execute if directly accessed by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Script is called by vSNP step 2

    Provide an Excel file with single column listing VCF file to remove from analysis.  This is convenient for controlling samples that may or may not be helpful in a comparison.  See it as a way to add or remove isolates without disrupting your VCF file database directory

    File names listed in the first column of the remove_from_analysis.xlsx can have extension, or the extension will be added (default .vcf).  Also _zc.vcf will be looked for and removed.
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r', '--excel_remove', action='store', dest='excel_remove', required=True, help='Excel file containing samples to remove from analysis. Column 1: to match sample name minus extension. No header allowed')
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Optional: path to VCF files')
    parser.add_argument('-e', '--extension', action='store', dest='extension', required=False, default="vcf", help='File extension type to be renamed')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()

    try:
        remove_from_analysis = Remove_From_Analysis(working_directory=args.working_directory, excel_remove=args.excel_remove, extension=args.extension)
        remove_from_analysis.remove_files()
    except KeyboardInterrupt:
        print(f'\n{bcolors.YELLOW}Process interrupted by user{bcolors.ENDC}')
        sys.exit(1)
    except Exception as e:
        print(f'{bcolors.RED}FATAL ERROR: {str(e)}{bcolors.ENDC}')
        sys.exit(1)