#!/usr/bin/env python

__version__ = "3.34"

import os
import re
import pandas as pd
import argparse
import textwrap
import sys


parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
Place description

'''), epilog='''---------------------------------------------------------''')

parser.add_argument('-t1', '--table1', action='store', dest='table1', required=True, help='Provide a vSNP SNP Table')
parser.add_argument('-t2', '--table2', action='store', dest='table2', required=True, default=None, help='Provide a vSNP SNP Table')
parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
args = parser.parse_args()

print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
print(args)
print("\n")

table_bwa = args.table1
table_minimap2 = args.table2

# More modern pandas approach with error handling
try:
    df_bwa = pd.read_excel(table_bwa, index_col=0)
    df_mp = pd.read_excel(table_minimap2, index_col=0)
except FileNotFoundError as e:
    print(f"Error: File not found - {e}", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"Error reading Excel files: {e}", file=sys.stderr)
    sys.exit(1)

# Get samples lists with more robust error handling
samples_bwa = df_bwa.index.tolist()
try:
    samples_bwa.remove('annotations')
except ValueError:
    pass

samples_bwa = [s.replace('_zc', '') for s in samples_bwa]

samples_mp = df_mp.index.tolist()
try:
    samples_mp.remove('annotation')
except ValueError:
    pass

samples_mp = [s.replace('_zc', '') for s in samples_mp]

# Get position lists
pos_bwa = df_bwa.columns.tolist()  # Use columns instead of iloc for clarity
pos_mp = df_mp.columns.tolist()

# Use set operations for better efficiency
sample_bwa_table_only = list(set(samples_bwa) - set(samples_mp))
sample_mp_table_only = list(set(samples_mp) - set(samples_bwa))
positions_bwa_table_only = list(set(pos_bwa) - set(pos_mp))
positions_mp_table_only = list(set(pos_mp) - set(pos_bwa))

# Prepare filenames
bwa_name = re.sub(r'[.].*', '', os.path.basename(table_bwa))
mp_name = re.sub(r'[.].*', '', os.path.basename(table_minimap2))
out_name = f'{bwa_name}--{mp_name}.tab'

# Use context manager for file handling with explicit encoding
try:
    with open(out_name, 'w', encoding='utf-8') as write_out:
        print('Samples found only in table1', file=write_out)
        if not sample_bwa_table_only:
            print('List is empty, All samples in table1 are in table2', file=write_out)
        else:
            for sample in sample_bwa_table_only:
                print(sample, file=write_out)
        print('----', file=write_out)

        print('Samples found only in table2', file=write_out)
        if not sample_mp_table_only:
            print('List is empty, All samples in table2 are in table1', file=write_out)
        else:
            for sample in sample_mp_table_only:
                print(sample, file=write_out)
        print('----', file=write_out)

        print('Positions found only in table1', file=write_out)
        if not positions_bwa_table_only:
            print('List is empty, All positions in table1 are in table2', file=write_out)
        else:
            for position in positions_bwa_table_only:
                print(position, file=write_out)
        print('----', file=write_out)

        print('Positions found only in table2', file=write_out)
        if not positions_mp_table_only:
            print('List is empty, All positions in table2 are in table1', file=write_out)
        else:
            for position in positions_mp_table_only:
                print(position, file=write_out)

except IOError as e:
    print(f"Error writing output file {out_name}: {e}", file=sys.stderr)
    sys.exit(1)

# Convert tab file to Excel with better error handling
try:
    # Read the tab file we just created with proper delimiter
    df = pd.read_csv(out_name, sep='\t', header=None, names=['Content'], encoding='utf-8')
    excel_filename = f'{bwa_name}--{mp_name}.xlsx'
    df.to_excel(excel_filename, index=False, engine='openpyxl')
    print(f"Output files created: {out_name} and {excel_filename}")
except Exception as e:
    print(f"Error converting tab file to Excel: {e}", file=sys.stderr)
    print(f"Tab file {out_name} was created successfully, but Excel conversion failed.")

# Created 2022 by Tod Stuber