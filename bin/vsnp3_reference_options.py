#!/usr/bin/env python

__version__ = "3.14"

import os
import sys
import glob
import re
import argparse
import textwrap
from collections import defaultdict


class Ref_Options():
    '''
    Directory locations can be added to reference_options_path.txt which is located in the dependencies folder of vSNP.  Paths of parent directory of dependency types are to be used.  Therefore when the following dependency types "species1" and "species2" are added the full working directory up to parent is provided as path.  Paths can be added manually or with path_adder.
    Example:
        /home/user/parent/species1
        /home/user/parent/species2
        Add: /home/user/parent to include species1 and species2 as reference options

    vsnp3_path_adder.py can be used to add paths.
    Usage:
        Change working diretory to folder containing reference files.
        vsnp3_path_adder.py -d `pwd`
    To see all available paths/reference:
        vsnp3_path_adder.py -s

    '''

    def __init__(self, select_ref=None):
        
        self.select_ref = select_ref
        self.defining_snps = None
        self.metadata = None
        self.gbk = None
        self.remove = None
        self.fasta = None
        all_ref_options = []
        script_path = os.path.dirname(os.path.realpath(__file__))
        self.script_path = script_path
        #don't use just the script path dependencies, but gather external dependency paths too
        ref_options_file = os.path.abspath(f'{script_path}/../dependencies/reference_options_paths.txt')
        self.ref_options_file = ref_options_file
        with open(f'{ref_options_file}', 'r') as dep_paths:
            dependency_paths = [line.strip() for line in dep_paths]
        #the additional dependency paths point to more reference options
        for path in dependency_paths:
            ref_options = glob.glob(f'{path}/*')
            all_ref_options = all_ref_options + ref_options
        all_ref_options = [x for x in all_ref_options if os.path.isdir(x)] #only capture directories
        self.all_ref_options = all_ref_options

        #if select_ref is a directory name then use files in it
        for directory in all_ref_options:
            if select_ref == directory.split('/')[-1]: #find reference by directory name
                self.select_ref = select_ref
                #if the asked for reference is a match grab files
                self.path = directory
                self.metadata_gather(directory)
            else:
                reference_header_capture = defaultdict(list) # find reference by fasta headers
                # for directory in all_ref_options:
                for each_fasta in glob.glob(f'{directory}/*fasta'):
                    with open(each_fasta, 'r') as each_fasta:
                        for each_line in each_fasta:
                            if each_line.startswith(">"):
                                reference_header_capture[directory].append(each_line.strip())
                for directory, header in reference_header_capture.items():
                    if select_ref in " ".join(header):
                        self.select_ref = select_ref
                        self.path = directory
                        self.metadata_gather(directory)

    def metadata_gather(self, directory):
        defining_snps = glob.glob(f'{directory}/*xlsx')
        defining_snps = [efile for efile in defining_snps if not re.search('~\$*', efile)] #ignore opened files
        # there are 3 excel files.  only 1 excel file for variable "excel".  it must be the non-*meta* and non-*remove* file
        defining_snps = [efile for efile in defining_snps if not re.search('.*remove.*', efile)] # just incase it is still in directory, remove it so only one excel is found below
        defining_snps = [efile for efile in defining_snps if not re.search('.*meta.*', efile)]
        #check that multiple files are not found for a single variable.  Each variable must point to just one file.
        if len(defining_snps) > 1:
            print(f'\n\n##### Exiting script {self.__eq__select_ref} contains more than one an Excel at {directory}\n')
            sys.exit(0)
        else:
            try:
                self.defining_snps = defining_snps[0]
            except IndexError:
                self.defining_snps = None
        # remove from analysis
        remove = glob.glob(f'{directory}/*remove*xlsx')
        remove = [efile for efile in remove if not re.search('~\$*', efile)] #ignore opened files
        if len(remove) > 1:
            print(f'\n\n##### Exiting script {self.select_ref} contains more than one remove file at {directory}\n')
            sys.exit(0)
        else:
            try:
                self.remove = remove[0]
            except IndexError:
                self.remove =  None
        #metadata
        metadata = glob.glob(f'{directory}/*meta*xlsx')
        metadata = [efile for efile in metadata if not re.search('~\$*', efile)] #ignore opened files
        if len(metadata) > 1:
            print(f'\n\n##### Exiting script {self.select_ref} contains more than one metadata file at {directory}\n')
            sys.exit(0)
        else:
            try:
                self.metadata = metadata[0]
            except IndexError:
                self.metadata =  None
        self.fasta = glob.glob(f'{directory}/*fasta')
        self.gbk = glob.glob(f'{directory}/*gbk')

    def files_in_directory(self):
        all_files = glob.glob(f'{self.path}/*')
        for each_file in all_files:
            print(f'{each_file}')
        print("")

    def print_options(self):
        each_reference_option=[]
        print("\nReference option files available:")
        print(f'Path are listed here: {self.ref_options_file}')
        for option in self.all_ref_options:
            each_reference_option.append((f'\t{os.path.split(option)[-1]}'))
        for each in sorted(each_reference_option):
            print(each)
        print("\nSee vsnp3_path_adder.py -h for more information\n")
        
        
if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Used by vSNP to get reference options by reference type

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-s', '--select_ref', action='store', dest='select_ref', required=True, help='Required: reference name')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    select_ref = args.select_ref
    ro = Ref_Options(select_ref)
    