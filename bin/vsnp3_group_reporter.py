#!/usr/bin/env python

__version__ = "3.34"

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
            print('VCF file and reference option must be provided')
            sys.exit(0)
                
        excel_path = reference_options.defining_snps
        xl = pd.ExcelFile(excel_path)
        sheet_names = xl.sheet_names
        ws = pd.read_excel(excel_path, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_dict = dict(zip(defining_snps.index, defining_snps.to_numpy()))
        
        defining_snps = {}
        inverted_defining_snps = {}
        mixed_defining_snps = {}  # New dictionary for mixed definitions
        
        for abs_pos, group in defsnp_dict.items():
            # Clean up the position string
            clean_abs_pos = abs_pos.replace('###', '')
            
            # Check if this is a mixed definition (contains both regular and inverted positions)
            positions = [pos.strip() for pos in clean_abs_pos.split(", ")]
            has_regular = any(not pos.endswith('!') for pos in positions)
            has_inverted = any(pos.endswith('!') for pos in positions)
            
            if has_regular and has_inverted:
                # Mixed definition - keep as is for special handling
                mixed_defining_snps[clean_abs_pos] = group
            elif has_inverted and not has_regular:
                # All positions are inverted
                inverted_pos = clean_abs_pos.replace('!', '')
                inverted_defining_snps[inverted_pos] = group
            else:
                # All positions are regular (no inversions)
                defining_snps[clean_abs_pos] = group
        
        self.defining_snps = defining_snps
        self.inverted_defining_snps = inverted_defining_snps
        self.mixed_defining_snps = mixed_defining_snps

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
                    except (TypeError, ValueError):
                        record_qual = 0 
                    chrom = record.CHROM
                    position = record.POS
                    absolute_position = f"{chrom}:{position}"
                    try:
                        # Safely access ALT column with proper checks
                        alt_value = record.ALT
                        if pd.notna(alt_value) and str(alt_value) != "None" and alt_value != "." and record.AC == AC and len(str(record.REF)) == 1 and record_qual > qual_threshold and record.MQ > MQ:
                            found_positions[absolute_position] = record.REF
                        if pd.notna(alt_value) and str(alt_value) != "None" and alt_value != "." and record.AC == 1 and len(str(record.REF)) == 1 and record_qual > qual_threshold and record.MQ > MQ:
                            found_positions_mix[absolute_position] = record.REF
                    except (KeyError, AttributeError, TypeError):
                        pass
                return filename, found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                return filename, 'see error', {}, {}
        except (SyntaxError, AttributeError) as e:
            # print(type(e)(str(e) + f'\n### VCF SyntaxError {filename} File Removed'))
            os.remove(filename)
            return filename, 'see error', {}, {}
    
    def bin_and_html_table(self, filename, found_positions, found_positions_mix):
        sample_groups_list = []
        tablename = os.path.basename(filename)
        defining_snps = self.defining_snps
        inverted_defining_snps = self.inverted_defining_snps
        mixed_defining_snps = self.mixed_defining_snps
        
        try:
            sample_groups_list = []
            defining_snp = False
            
            # Convert to sets for efficient operations
            found_positions_set = set(found_positions.keys())
            found_positions_mix_set = set(found_positions_mix.keys())
            all_found_positions = found_positions_set.union(found_positions_mix_set)
            
            # Check regular defining SNPs (may contain multiple positions)
            for abs_position, group in defining_snps.items():
                # Split positions by comma and strip whitespace
                required_positions = [pos.strip() for pos in abs_position.split(", ")]
                
                # Check if this is a multi-position requirement
                if len(required_positions) > 1:
                    # For multiple positions, all must be found
                    required_positions_set = set(required_positions)
                    if required_positions_set.issubset(all_found_positions):
                        sample_groups_list.append(group)
                        defining_snp = True
                        if found_positions_mix_set.intersection(required_positions_set):
                            tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'
                else:
                    # Single position check
                    if abs_position in all_found_positions:
                        sample_groups_list.append(group)
                        defining_snp = True
                        if abs_position in found_positions_mix_set:
                            tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'
            
            # Check inverted defining SNPs (positions that should NOT be present)
            for abs_position, group in inverted_defining_snps.items():
                # Split positions by comma and strip whitespace
                positions_to_check = [pos.strip() for pos in abs_position.split(", ")]
                
                # For inverted positions, NONE of them should be found
                positions_set = set(positions_to_check)
                if not positions_set.intersection(all_found_positions):
                    sample_groups_list.append(group)
                    defining_snp = True
            
            # Handle mixed regular and inverted positions in the same group definition
            # This handles cases like: "pos1, pos2, pos3!" where pos1 and pos2 must be present, pos3 must be absent
            for abs_position, group in mixed_defining_snps.items():
                positions = [pos.strip() for pos in abs_position.split(", ")]
                
                # Separate regular and inverted positions
                regular_positions = []
                inverted_positions = []
                
                for pos in positions:
                    if pos.endswith('!'):
                        inverted_positions.append(pos[:-1])  # Remove the '!' for checking
                    else:
                        regular_positions.append(pos)
                
                regular_positions_set = set(regular_positions) if regular_positions else set()
                inverted_positions_set = set(inverted_positions) if inverted_positions else set()
                
                # Check conditions based on what we have
                regular_match = True  # Default to True if no regular positions
                inverted_match = True  # Default to True if no inverted positions
                
                if regular_positions:
                    regular_match = regular_positions_set.issubset(all_found_positions)
                
                if inverted_positions:
                    inverted_match = not inverted_positions_set.intersection(all_found_positions)
                
                if regular_match and inverted_match:
                    sample_groups_list.append(group)
                    defining_snp = True
                    if regular_positions and found_positions_mix_set.intersection(regular_positions_set):
                        tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'

            if not defining_snp:  # extra step to get the group name when there are multiple defining snps for a group
                for abs_position in defining_snps.keys():
                    set_abs_position = set(abs_position.split(", "))
                    is_subset = set_abs_position.issubset(found_positions_set)
                    if is_subset:
                        group = defining_snps[abs_position]
                        sample_groups_list.append(group)

            if not sample_groups_list:
                sample_groups_list = ['No defining SNPs']
            else:
                sample_groups_list = sorted(list(set(sample_groups_list)))  # Remove duplicates and sort

        except TypeError:
            message = 'File TypeError'
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