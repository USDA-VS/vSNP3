#!/usr/bin/env python

__version__ = "3.28"

import os
import io
import sys
import pandas as pd
import argparse
import textwrap

from vsnp3_reference_options import Ref_Options


class GroupReporter:

    def __init__(self, vcf, reference_type=None):
        if reference_type and vcf:
            reference_options = Ref_Options(reference_type)
            self.vcf = vcf
        else:
            print(f'VCF file and reference option must be provided')
            sys.exit(0)
                
        excel_path = reference_options.defining_snps
        xl = pd.ExcelFile(excel_path)
        sheet_names = xl.sheet_names
        ws = pd.read_excel(excel_path, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_dict = dict(zip(defining_snps.index, defining_snps.to_numpy()))
        defining_snps={}
        inverted_defining_snps={}
        for abs_pos, group in defsnp_dict.items():
            if '!' in abs_pos:
                inverted_defining_snps[abs_pos.replace('!', '')] = group
            else:
                defining_snps[abs_pos.replace('###', '')] = group #capture groups blocked out by ###
        self.defining_snps = defining_snps
        self.inverted_defining_snps = inverted_defining_snps

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
        # Using a more efficient approach to extract values from INFO column
        def extract_info_field(info_str, field):
            try:
                info_dict = {item.split('=')[0]: item.split('=')[1] for item in info_str.split(';') if '=' in item}
                return info_dict.get(field, None)
            except (IndexError, AttributeError):
                return None
                
        df['AC'] = df['INFO'].apply(lambda x: extract_info_field(x, 'AC'))
        df['DP'] = df['INFO'].apply(lambda x: extract_info_field(x, 'DP'))
        df['MQ'] = df['INFO'].apply(lambda x: extract_info_field(x, 'MQ'))
        df['AC'] = pd.to_numeric(df['AC'], errors='coerce').fillna(0).astype(int)
        df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype(int)
        df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce').fillna(0).astype(int)
        df = df.drop(columns=['INFO', 'ID', 'FILTER', 'FORMAT'])

        return df

    def find_initial_positions(self, filename):
        found_positions = {}
        found_positions_mix = {}
        # Values must be hardcoded.  This is used in step1 where options are not used.  Defining SNPs should be of relatively high quality.
        AC = 2
        qual_threshold = 150
        MQ = 56
        try:
            df = self.read_vcf(filename)
            try:
                for index, record in df.iterrows():
                    try:
                        record_qual = int(record.QUAL)
                    except TypeError:
                        record_qual = 0 
                    chrom = record.CHROM
                    position = record.POS
                    absolute_positon = str(chrom) + ":" + str(position)
                    try:
                        if str(record.ALT[0]) != "None" and record.AC == AC and len(record.REF) == 1 and record_qual > qual_threshold and record.MQ > MQ:
                            found_positions[absolute_positon] = record.REF
                        if str(record.ALT[0]) != "None" and record.AC == 1 and len(record.REF) == 1 and record_qual > qual_threshold and record.MQ > MQ:
                            found_positions_mix[absolute_positon] = record.REF
                    except KeyError:
                        pass
                return filename, found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                return filename, f'see error', {}, {}
        except (SyntaxError, AttributeError) as e:
            # print(type(e)(str(e) + f'\n### VCF SyntaxError {filename} File Removed'))
            os.remove(filename)
            return filename, f'see error', {}, {}
    
    def bin_and_html_table(self, filename, found_positions, found_positions_mix):
        sample_groups_list = []
        tablename = os.path.basename(filename)
        defining_snps = self.defining_snps
        inverted_defining_snps = self.inverted_defining_snps
        try:
            sample_groups_list = []
            defining_snp = False
            # Using set operations more efficiently
            defining_positions = set(defining_snps.keys())
            found_positions_set = set(found_positions.keys())
            found_positions_mix_set = set(found_positions_mix.keys())
            
            # Check intersections with both found position sets
            matching_positions = defining_positions.intersection(found_positions_set.union(found_positions_mix_set))
            
            for abs_position in matching_positions:
                group = defining_snps[abs_position]
                sample_groups_list.append(group)
                if defining_positions.intersection(found_positions_mix_set):
                    tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'
                defining_snp = True
                
            # Check for inverted defining SNPs
            inverted_positions = set(inverted_defining_snps.keys())
            if not inverted_positions.intersection(found_positions_set.union(found_positions_mix_set)):
                for abs_position in inverted_positions:
                    group = inverted_defining_snps[abs_position]
                    sample_groups_list.append(group)
                    defining_snp = True

            if not defining_snp:  # extra step to get the group name when there are multiple defining snps for a group
                for abs_position in defining_positions:
                    set_abs_position = set(abs_position.split(", "))
                    is_subset = set_abs_position.issubset(found_positions_set)
                    if is_subset:
                        group = defining_snps[abs_position]
                        sample_groups_list.append(group)

            if not sample_groups_list:
                sample_groups_list = ['No defining SNPs']
            else:
                sample_groups_list = sorted(sample_groups_list)

        except TypeError:
            message = f'File TypeError'
            print(f'{message}: {filename}')
            sample_groups_list = [f'{message}: {filename}']
            
        return sample_groups_list

    def get_groups(self):
        filename, found_positions, found_positions_mix = self.find_initial_positions(self.vcf)
        sample_groups_list = self.bin_and_html_table(filename, found_positions, found_positions_mix)
        return sample_groups_list

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    
    Used by vsnp3 step 1 to report on group

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-vcf', '--vcf', action='store', dest='vcf', required=True, help='Required: vcf file')
    parser.add_argument('-r', '--reference_type', action='store', dest='reference_type', required=True, default=None, help='Provide reference option.  See vsnp3_path_adder.py -s for options')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
 
    group_reporter = GroupReporter(args.vcf, args.reference_type)
    sample_groups_list = group_reporter.get_groups()
    print(sample_groups_list)
    
    print("Done")