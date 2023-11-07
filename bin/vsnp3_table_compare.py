#!/usr/bin/env python

__version__ = "3.17"

import os
import re
import pandas as pd
import argparse
import textwrap


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
df_bwa=pd.read_excel(table_bwa, index_col=0)
df_mp=pd.read_excel(table_minimap2, index_col=0)
samples_bwa = df_bwa.index.to_list()
try:
    samples_bwa.remove('annotations')
except ValueError:
    pass
samples_bwa = [s.replace('_zc', '') for s in samples_bwa]
samples_mp = df_mp.index.to_list()
try:
    samples_mp.remove('annotation')
except ValueError:
    pass
samples_mp = [s.replace('_zc', '') for s in samples_mp]
pos_bwa = df_bwa.iloc[0].index.to_list()
pos_mp = df_mp.iloc[0].index.to_list()

sample_bwa_table_only=[x for x in samples_bwa if x not in set(samples_mp)]
sample_mp_table_only=[x for x in samples_mp if x not in set(samples_bwa)]
positions_bwa_table_only=[x for x in pos_bwa if x not in set(pos_mp)]
positions_mp_table_only=[x for x in pos_mp if x not in set(pos_bwa)]

bwa_name = re.sub('[.].*', '', os.path.basename(table_bwa))
mp_name = re.sub('[.].*', '', os.path.basename(table_minimap2))
out_name = f'{bwa_name}--{mp_name}.tab'
write_out = open(out_name, 'w')

print(f'Samples found only in table1', file=write_out)
if not sample_bwa_table_only:
    print(f'List is empty, All samples in table1 are in table2', file=write_out)
else:
    for i in sample_bwa_table_only:
        print(f'{i}', file=write_out)
print(f'----', file=write_out)


print(f'Samples found only in table2', file=write_out)
if not sample_mp_table_only:
    print(f'List is empty, All samples in table2 are in table1', file=write_out)
else:
    for i in sample_mp_table_only:
        print(f'{i}', file=write_out)
print(f'----', file=write_out)

print(f'Positions found only in table1', file=write_out)
if not positions_bwa_table_only:
    print(f'List is empty, All postions in table1 are in table2', file=write_out)
else:
    for i in positions_bwa_table_only:
        print(f'{i}', file=write_out)
print(f'----', file=write_out)

print(f'Positions found only in table2', file=write_out)
if not positions_mp_table_only:
    print(f'List is empty, All postions in table2 are in table1', file=write_out)
else:
    for i in positions_mp_table_only:
        print(f'{i}', file=write_out)

write_out.close()
df = pd.read_csv(out_name, on_bad_lines='skip')
df.to_excel(f'{bwa_name}--{mp_name}.xlsx', index=False)

# Created 2022 by Tod Stuber
