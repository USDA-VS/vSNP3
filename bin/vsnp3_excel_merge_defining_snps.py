#!/usr/bin/env python

__version__ = "0.0.1"

import os
import re
import sys
import math
import pandas as pd
from openpyxl.styles import Font
from openpyxl import Workbook
from datetime import datetime
import argparse
import textwrap

import locale
# Set the locale to United States
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

class Expand_Range():
    ''' 
    '''
    def ex_range(self, file=None):

        file_sample_name = re.sub('[.].*', '', file)
        # date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        expanded_file_name = f'{file_sample_name}_expanded.xlsx'

        df = pd.read_excel(file)
        # dictionary for incoming lists
        column_lists = {}
        # Iterate over each column and store values in lists
        for column in df.columns:
            column_lists[column] = df[column].tolist()

        # Print the lists
        update_ranges=[]
        for column, values in column_lists.items():
            expansion=[]
            group_name = values[0]
            reference = values[1].split(':')[0]
            string_strip = f'{reference}:'
            values = values[1:]
            values = [x for x in values if x is not None and not (isinstance(x, float) and math.isnan(x))]
            items_without_matching_prefix = [item for item in values if not item.startswith(string_strip)]
            if items_without_matching_prefix:
                print("### Position entries have wrong reference prefix.\n### There may be mixed reference types.\n### Must check worksheet and remove additional references\n")
                print(f'### The correct reference type appears to be: {reference}')
                print('### Positions not matching:')
                for item in items_without_matching_prefix:
                    print(item, end=', ')
                print("\n\n### ERROR:  Must Fix to Continue\n")
                sys.exit(1)
            values = [item.replace(string_strip, '') for item in values]
            for value in values:
                if "-" in value:
                    start, end = map(int, value.split('-'))
                    expanded_range = list(range(start, end + 1))
                    expansion = list(set(expansion + expanded_range))
                else:
                    expansion.append(int(value))
            expansion = list(set(expansion))
            expansion.sort()
            expansion = [f'{string_strip}{item}' for item in expansion]
            # print(f"\ndefining snp '{column}', group {group_name}\n{expansion}")
            expansion.insert(0, group_name)
            if column == 'Unnamed: 0':
                expansion.insert(0, '')
            else:
                expansion.insert(0, column)
            update_ranges.append(expansion)

        # Create a DataFrame from the list of lists
        df = pd.DataFrame(update_ranges).transpose()

        # Write the DataFrame to an Excel file
        df.to_excel(expanded_file_name, index=False, header=False)
        return expanded_file_name

class Merge_Defining_SNPs():
    ''' 
    '''
    def __init__(self, file1=None, file2=None):

        def numbers_to_ranges(numbers):
            numbers = [int(x) for x in numbers]
            numbers.sort()
            ranges = []
            start = numbers[0]

            for i in range(1, len(numbers)):
                if numbers[i] != numbers[i-1] + 1:
                    if start == numbers[i-1]:
                        ranges.append(str(start))
                    else:
                        ranges.append(f"{start}-{numbers[i-1]}")
                    start = numbers[i]

            # Handle the last number or range
            if start == numbers[-1]:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{numbers[-1]}")
            ranges = sorted(ranges, key=lambda x: int(x.split('-')[0]))
            return ranges

        def sort_df(df):
            def search_list_for_pattern(lst, pattern):
                for idx, item in enumerate(lst):
                    if pattern in item:
                        break
                return idx
            
            second_row = df.iloc[1, :].tolist()
            all_index = search_list_for_pattern(second_row, "-All") # search for "All" column to move left
            sorted_columns = sorted(range(len(second_row)), key=lambda k: second_row[k])
            # Move the target column to the first position
            sorted_columns.remove(df.columns.get_loc(all_index))
            sorted_columns.insert(0, df.columns.get_loc(all_index))
            # Reorder the columns of the DataFrame based on the sorted list
            df = df.iloc[:, sorted_columns]
            return df
        
        file_sample_name1 = re.sub('_expanded.xlsx', '', file1)
        file_sample_name2 = re.sub('_expanded.xlsx', '', file2)
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

        auto_merge=[]

        df = pd.read_excel(file1)
        # dictionary for incoming lists
        column_dict1 = {}
        # Iterate over each column and store values in lists
        for column in df.columns:
            column_dict1[column] = df[column].tolist()

        df = pd.read_excel(file2)
        # dictionary for incoming lists
        column_dict2 = {}
        # Iterate over each column and store values in lists
        for column in df.columns:
            column_dict2[column] = df[column].tolist()

        keys_set1 = set(column_dict1.keys())
        keys_set2 = set(column_dict2.keys())

        defining_in_both = keys_set1.intersection(keys_set2)
        defining_only_in1_not2 = keys_set1.difference(keys_set2)
        defining_only_in2_not1 = keys_set2.difference(keys_set1)

        defining_same_with_diff_positions=[]
        for each in defining_in_both:
            values1 = column_dict1[each]
            group_name1 = values1[0]
            reference1 = values1[1].split(':')[0]
            string_strip = f'{reference1}:'
            values1 = [x for x in values1 if x is not None and not (isinstance(x, float) and math.isnan(x))]
            values1 = [item.replace(string_strip, '') for item in values1]
            values1 = set(values1[1:])

            values2 = column_dict2[each]
            group_name2 = values2[0]
            values2 = [x for x in values2 if x is not None and not (isinstance(x, float) and math.isnan(x))]
            values2 = [item.replace(string_strip, '') for item in values2]
            values2 = set(values2[1:])

            in_list1_only = list(values1.difference(values2))
            if in_list1_only:
                in_list1_only = numbers_to_ranges(in_list1_only)
                in_list1_only = [f'{string_strip}{item}' for item in in_list1_only]
                in_list1_only.insert(0, f'Only in {file_sample_name1}: {group_name1}')
                if each == 'Unnamed: 0':
                    in_list1_only.insert(0, '')
                else:
                    in_list1_only.insert(0, each)
                defining_same_with_diff_positions.append(in_list1_only)

            in_list2_only = list(values2.difference(values1))
            if in_list2_only:
                in_list2_only = numbers_to_ranges(in_list2_only)
                in_list2_only = [f'{string_strip}{item}' for item in in_list2_only]
                in_list2_only.insert(0, f'Only in {file_sample_name2}: {group_name2}')
                if each == 'Unnamed: 0':
                    in_list2_only.insert(0, '')
                else:
                    in_list2_only.insert(0, each)
                defining_same_with_diff_positions.append(in_list2_only)

            combined_set = list(values1.union(values2))
            if combined_set:
                combined_set = numbers_to_ranges(combined_set)
                combined_set = [f'{string_strip}{item}' for item in combined_set]
                combined_set.insert(0, f'{group_name1}')
                if each == 'Unnamed: 0':
                    combined_set.insert(0, '')
                else:
                    combined_set.insert(0, each)
                auto_merge.append(combined_set)
        # Create a DataFrame from the list of lists
        dfall = pd.DataFrame(defining_same_with_diff_positions,).transpose()
        dfall = sort_df(dfall)

        pos_only_in1=[]
        for each in defining_only_in1_not2:
            values1 = column_dict1[each]
            values1 = [x for x in values1 if x is not None and not (isinstance(x, float) and math.isnan(x))]
            values1 = [item.replace(string_strip, '') for item in values1]
            group_name1 = values1[0]
            values1 = list(set(values1[1:]))
            values1 = numbers_to_ranges(values1)
            values1 = [f'{string_strip}{item}' for item in values1]
            values1.insert(0, group_name1)
            if each == 'Unnamed: 0':
                values1.insert(0, '')
            else:
                values1.insert(0, each)
            pos_only_in1.append(values1)
            auto_merge.append(values1)
            df1 = pd.DataFrame(pos_only_in1).transpose()
            df1 = sort_df(df1)

        pos_only_in2=[]
        for each in defining_only_in2_not1:
            values2 = column_dict2[each]
            values2 = [x for x in values2 if x is not None and not (isinstance(x, float) and math.isnan(x))]
            values2 = [item.replace(string_strip, '') for item in values2]
            group_name2 = values2[0]
            values2 = list(set(values2[1:]))
            values2 = numbers_to_ranges(values2)
            values2 = [f'{string_strip}{item}' for item in values2]
            values2.insert(0, group_name2)
            if each == 'Unnamed: 0':
                values2.insert(0, '')
            else:
                values2.insert(0, each)
            pos_only_in2.append(values2)
            auto_merge.append(values2)
            df2 = pd.DataFrame(pos_only_in2).transpose()
            df2 = sort_df(df2)

        merged_file_name = f'merge_differences_{file_sample_name1}_{file_sample_name2}_{date_stamp}.xlsx'
        with pd.ExcelWriter(merged_file_name) as diff_writer:
            # Write each DataFrame to a separate sheet with a specific name
            dfall.to_excel(diff_writer, sheet_name='Defin_SNP_in_both_but_diff_pos', index=False, header=False)
            df1.to_excel(diff_writer, sheet_name=f'Defin_SNP_only_{file_sample_name1}'[:31], index=False, header=False)
            df2.to_excel(diff_writer, sheet_name=f'Defin_SNP_only_{file_sample_name2}'[:31], index=False, header=False)

        dfauto = pd.DataFrame(auto_merge).transpose()
        dfauto = sort_df(dfauto)

        merged_auto = f'merge_auto_{file_sample_name1}_{file_sample_name2}_{date_stamp}.xlsx'
        with pd.ExcelWriter(merged_auto, engine='openpyxl') as write_auto:
            # Write an updated defining snps Excel file with all combined from both worksheets
            dfauto.to_excel(write_auto, sheet_name='Sheet1', index=False, header=False)
            # Access the worksheet
            worksheet = write_auto.sheets['Sheet1']

            # Apply bold font to the specified rows
            bold_font = Font(bold=True)
            for cell in worksheet[2]: #bold row 2
                cell.font = bold_font

# sort merge auto dataframe by second row

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f1', '--file1', action='store', dest='file1', required=False, help='merge 1 file')
    parser.add_argument('-f2', '--file2', action='store', dest='file2', required=False, default=None, help='merge 2 file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    #file1
    expand_range = Expand_Range()
    file1_expanded = expand_range.ex_range(args.file1)

    #file2
    expand_range = Expand_Range()
    file2_expanded = expand_range.ex_range(args.file2)

    merge_defining_snps = Merge_Defining_SNPs(file1_expanded, file2_expanded)

    try:
        os.remove(file1_expanded)
    except OSError as e:
        print(f"Error occurred while deleting the file: {e}")
    try:
        os.remove(file2_expanded)
    except OSError as e:
        print(f"Error occurred while deleting the file: {e}")