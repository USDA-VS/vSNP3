#!/usr/bin/env python

__version__ = "3.15"

import os
import sys
import pandas as pd
import glob
import argparse
import textwrap

class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE='\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

class Remove_From_Analysis:
    '''
    ''' 
    def __init__(self, working_directory='.', excel_remove=None, extension=None):
        if working_directory == '.':
            working_directory = os.getcwd()
        if list(working_directory)[0] == '/':
            print(f'working directory: {working_directory}')
        else:
            print(f'##### PROVIDE A FULL PATH')
            print(f'directory given: "{working_directory}"')
            sys.exit(0)

        df = pd.read_excel(excel_remove, index_col=0, usecols=[0], header=None)
        remove_list=[]
        for each_sample in df.index:
            remove_list.append(f'{working_directory}/{each_sample}') #if .vcf is supplied with the sample name in remove_from_analysis.xlsx
            remove_list.append(f'{working_directory}/{each_sample}.{extension}') #most common behavior if sample name without extension is provided in remove_from_analysis.xlsx
            remove_list.append(f'{working_directory}/{each_sample}_zc.{extension}') #allow _zc to not be specified in remove_from_analysis.xlsx
        self.excel_remove = excel_remove
        self.remove_list = remove_list

    def remove_files(self):
        num_files_removed = 0
        print(f'Removing samples listed in {self.excel_remove}')
        for each_sample in self.remove_list:
            glob_list = glob.glob(each_sample)
            for item in glob_list:
                num_files_removed += 1
                # print(f'\tRemoving: {item}')
                os.remove(item)
        print(f'\n{bcolors.RED}{num_files_removed:,} {bcolors.ENDC}{bcolors.WHITE}files removed from the analysis{bcolors.ENDC}\n')
        self.removed_file_count = num_files_removed


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Script is called by vSNP step 2

    Provide an Excel file with single column listing VCF file to remove from analysis.  This is convenient for controlling samples that may or may not be helpful in a comparison.  See it as a way to add or remove isolates without distrupting your VCF file database directory

    File names listed in the first column of the remove_from_analysis.xlsx can have extension, or the extension will be added (default .vcf).  Also _zc.vcf will be looked for and removed.
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r', '--excel_remove', action='store', dest='excel_remove', required=True, help='Excel file containing samples to remove from analysisColumn 1: to match sample name minus extension.No header allowed')
    parser.add_argument('-w', '--cwd', action='store', dest='working_directory', required=False, default='.', help='Optional: path to VCF files')
    parser.add_argument('-e', '--extension', action='store', dest='extension', required=False, default="vcf", help='File extension type to be renamed')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()

    remove_from_analysis = Remove_From_Analysis(working_directory=args.working_directory, excel_remove=args.excel_remove, extension=args.extension)
    remove_from_analysis.remove_files()