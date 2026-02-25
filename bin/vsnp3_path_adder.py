#!/usr/bin/env python

__version__ = "3.34"

import os
import glob
import argparse
import textwrap


class Add_Path:
    '''
    Dependencies:
        Defining SNP/Filter file, xlsx - REQUIRED
        FASTA - REQUIRED
        GBK - OPTIONAL
        GFF - OPTIONAL
        metadata.xlsx (2 column: filename and updated name) - OPTIONAL
    '''

    def __init__(self, added_directory):
        if added_directory == '.':
            added_directory = os.getcwd()
        else:
            added_directory = added_directory
        try: 
            if list(added_directory)[0] == '/':
                print(f'\nDirectory Given: {added_directory}\n')
                self.valid_path = True
            else:
                self.valid_path = False
        except (TypeError, IndexError):
            # If None is passed error because not iterable, or empty string
            self.valid_path = False
        script_path = os.path.dirname(os.path.realpath(__file__))
        self.script_path = script_path
        abspath = os.path.abspath(os.path.join(script_path, '../dependencies/reference_options_paths.txt'))
        self.abspath = abspath
        self.added_directory = added_directory

    #insure no blank lines are in file
    def remove_blank_lines(self):
        try:
            with open(self.abspath, "r") as open_file:
                all_lines = open_file.readlines()
            with open(self.abspath, "w") as write_back:  
                [write_back.write(line) for line in all_lines if line.strip()]
        except (IOError, OSError) as e:
            print(f"Error accessing file {self.abspath}: {e}")

    def add_to_path(self):
        #append new path to file
        try:
            with open(self.abspath, 'a') as dep_paths:
                print(self.added_directory, file=dep_paths)

            # Remove duplicate lines from file
            with open(self.abspath, "r") as infile:
                all_lines = [line.rstrip() for line in infile]
            unique_lines = set(all_lines)
            with open(self.abspath, "w") as outfile:
                for each_line in unique_lines:
                    print(each_line, file=outfile)
        except (IOError, OSError) as e:
            print(f"Error writing to file {self.abspath}: {e}")

    def show_options(self):
        # Show reference options
        try:
            with open(self.abspath, 'r') as dep_paths:
                reference_options_paths = [line.rstrip() for line in dep_paths]
        except (IOError, OSError) as e:
            print(f"Error reading file {self.abspath}: {e}")
            return
            
        for path in reference_options_paths:
            print(f'Path: {path}')
            try:
                ref_options = glob.glob(os.path.join(path, '*'))
                ref_options = [x for x in ref_options if os.path.isdir(x)]
                if ref_options:
                    for option in sorted(ref_options):
                        print(f'\t{os.path.basename(option)}')
                    print('\n')
                elif not os.path.exists(path):
                    print('\tPath does not exist\n')
                else:
                    print('\tPath is empty\n')
            except (OSError, IOError) as e:
                print(f'\tError accessing path {path}: {e}\n')
             

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''
    ---------------------------------------------------------

    vsnp3_path_adder.py is used to add reference types.
    Usage:
        Change working diretory to directory one up from directory containing dependencies.
        For example if AF2122 and H37 need to be added and are located here
        ~/Mycobacterium/AF2122
        ~/Mycobacterium/H37
        Change to ~/Mycobacterium and run
        vsnp3_path_adder.py -d $(pwd)
        This will add both AF2122 and H37 as reference types.
    To see all available paths/reference:
        vsnp3_path_adder.py -s

    Files to have in each reference type:
    - FASTA (required)
    - GBK (optional)
    - define_filter.xlsx (required)
    - remove_from_analysis.xlsx (optional)

    See vsnp3_download_fasta_gbk_gff_by_acc.py to download FASTA and GBK
    define_filter.xls and remove_from_analysis.xlsx templates are found in vsnp3 dependencies directory

    ---------------------------------------------------------
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-d', '--cwd', action='store', dest='directory', required=False, help='Absolute directory path to be added to find reference option files.')
    parser.add_argument('-s', '--show', action='store_true', dest='show', required=False, help='Show available directories.')
    
    # Fixed version argument to avoid potential f-string issues
    version_info = f'{os.path.abspath(__file__)}: version {__version__}'
    parser.add_argument('-v', '--version', action='version', version=version_info)
    args = parser.parse_args()

    if args.directory:
        directory = args.directory
        add_path = Add_Path(directory)
        add_path.remove_blank_lines()
        if add_path.valid_path:
            add_path.add_to_path()
            add_path.show_options()
            print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')
        else:
            print('##### PROVIDE AN ABSOLUTE PATH')
            print(f'Directory Given: "{add_path.added_directory}"')
            print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')
    else:
        print('\nNO NEW DIRECTORY ADDED\nONLY SHOWING WHAT IS CURRENTLY AVAILABLE')
        add_path = Add_Path(None)
        add_path.remove_blank_lines()
        try:
            if os.stat(add_path.abspath).st_size == 0:
                print(f'\nPaths have not been added.  File is empty. \nAdd path using -d option.\nSee -h for more options\n\t{add_path.abspath}\n')
            else:
                add_path.show_options()
                print(f'\nPaths files can be manually edited here: {add_path.abspath}\n')
        except (OSError, IOError) as e:
            print(f'Error accessing file {add_path.abspath}: {e}')