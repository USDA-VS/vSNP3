#!/usr/bin/env python

__version__ = "3.32"

import os
import subprocess
import sys
import re
import glob
import time
import random
import locale
from datetime import datetime
import pandas as pd
import multiprocessing
from concurrent import futures
import argparse
import textwrap
import itertools
from Bio import SeqIO
from collections import defaultdict

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class Tree:
    ''' 
    '''

    def __init__(self, fasta_alignments=None, write_path=None, tree_name=None, debug=False,):
        # Find an optimal compiled version of RAxML in conda
        try:
            subprocess.call("raxml", stdout=subprocess.DEVNULL)
            raxml = "raxml"
        except OSError:
            try:
                subprocess.call("raxmlHPC-PTHREADS-AVX2", stdout=subprocess.DEVNULL)
                raxml = "raxmlHPC-PTHREADS-AVX2"
            except OSError:
                try:
                    subprocess.call("raxmlHPC-PTHREADS-AVX", stdout=subprocess.DEVNULL)
                    raxml = "raxmlHPC-PTHREADS-AVX"
                except OSError:
                    try:
                        subprocess.call("raxmlHPC-PTHREADS", stdout=subprocess.DEVNULL)
                        raxml = "raxmlHPC-PTHREADS"
                    except OSError:
                        try:
                            subprocess.call("raxmlHPC-SSE3", stdout=subprocess.DEVNULL)
                            raxml = "raxmlHPC-SSE3"
                        except OSError:
                            try:
                                subprocess.call("raxmlHPC", stdout=subprocess.DEVNULL)
                                raxml = "raxmlHPC"
                            except OSError:
                                raxml = 'raxmlHPC'
            # print('set RAxML to {}'.format(raxml))
        self.raxml = raxml
        self.cpu_count = int(multiprocessing.cpu_count() / 1.2)
        # self.hash_names = hash_names
        self.debug = debug
        
        self.beginTime = datetime.now()
        self.startTime = datetime.now()
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')

        if not write_path:
            write_path = os.getcwd()
        if not tree_name:
            tree_name = "raxml"
        self.write_path = write_path
        self.tree_name = tree_name

        raxml_version = subprocess.check_output('{} -v'.format(raxml), shell=True, text=True)
        self.raxml_version = raxml_version.split('\n')[2]

        os.system('{} -s {} -n raxml -m GTRCATI -o root -w {} -p 456123 -T 4 > /dev/null 2>&1'.format(raxml, fasta_alignments, write_path)) #> /dev/null 2>&1
        # os.system('{} -s {} -n raxml -m ASC_GTRGAMMA --asc-corr lewis  -o root -w {} -p 456123 -T 4 > /dev/null 2>&1'.format(raxml, fasta_alignments, write_path)) #> /dev/null 2>&1
        try:
            newick = os.path.join(write_path, '{}_{}.tre'.format(tree_name, self.st))
            os.rename(os.path.join(write_path, 'RAxML_bestTree.raxml'), newick)
            raxml_to_remove = glob.glob(os.path.join(write_path, 'RAxML*'))
            for each in raxml_to_remove:
                os.remove(each)
            try:
                reduced_file = glob.glob(os.path.join(write_path, '*.reduced'))
                os.remove(reduced_file[0])
            except (FileNotFoundError, IndexError) as e:
                pass
        except FileNotFoundError:
            with open(os.path.join(write_path, 'SEE_RAXML_INFO'), 'w') as message_out:
                print('check sample numbers', file=message_out)
        self.newick = newick
    
    def checksum_match_to_text(self, tree):
        # read entire tree into variable as string obj
        with open(tree, 'rt') as open_tree:
            entire_file = open_tree.read()
            print (entire_file)
        with open("idtable.txt" , 'rt') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split("\t")
                org = str(line[1])
                org = org.rstrip()
                check = str(line[0])
                check = check.rstrip()
                print (org)
                print (check)
                entire_file = re.sub(check, "'" + org + "'", entire_file)
            f.close()
        outfile = "NAMES-UPDATED-" + os.path.basename(tree)
        write_out = open(outfile , 'wt')
        write_out.write(entire_file)
        write_out.close()
        return outfile

    def df_to_fasta(self, df, alignment_file):
        test_duplicates=[] # if duplicate name in alignment fasta raxml with error and exit
        with open(alignment_file, 'w') as write_out:
            for index, row in df.iterrows():
                test_duplicates.append(row.name)
                if test_duplicates.count(row.name) < 2:
                    print('>{}'.format(row.name), file=write_out)
                    for pos in row:
                        print(pos, end='', file=write_out)
                    print("", file=write_out)

    def get_parsimonious_pos(self, in_df):
        try:
            ref_series = in_df.loc['reference_seq']
            in_df = in_df.drop(['reference_seq']) #in all_vcf reference_seq needs to be removed
        except KeyError:
            print('Check that there is a "reference_seq" nameed')
            sys.exit(0)
        # print('in_df size: {}'.format(in_df.shape))
        parsimony = in_df.loc[:, (in_df != in_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = in_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        return out_df


class Tables:

    def __init__(self, fasta_alignments=None, df_alignments=None, tree=None, gbk=None, mq=None, write_path=None, groupings_dict=None, show_groups=False, table_name=None, debug=False, sample_coverage_dict=None):
        self.fasta_alignments = fasta_alignments
        self.df_alignments = df_alignments
        self.tree = tree
        self.gbk = gbk
        self.mq = mq
        self.show_groups = show_groups
        self.debug = debug
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
        self.groupings_dict = groupings_dict
        self.sample_coverage_dict = sample_coverage_dict or {}
                
        if self.groupings_dict:
            group_vcfs_dict = defaultdict(list) #invert the key, values
            for group, dataframes in self.groupings_dict.items():
                for vcf in dataframes.keys():
                    group_vcfs_dict[vcf].append(group)
            self.group_vcfs_dict = dict(sorted(group_vcfs_dict.items())) #sorts on key(vcf name)

        if not write_path:
            write_path = os.getcwd()
        if not table_name:
            table_name = "generic"
        self.write_path = write_path
        self.table_name = table_name

    def build_tables(self, ):
        if self.df_alignments is None:
            reformated = os.path.join(self.write_path, 'reformated.fasta')
            with open(reformated, 'w') as reformat:
                sequence = SeqIO.parse(self.fasta_alignments, "fasta")
                for each in sequence:
                    print('>{}\n{}'.format(each.description, each.seq), file=reformat)
            fasta_df = pd.read_csv(reformated, header=None, sep='^')
            os.remove(reformated)
            seq = fasta_df.iloc[1::2].reset_index(drop=True)
            header = fasta_df.iloc[0::2].reset_index(drop=True).replace(to_replace=r'>', value='', regex=True)
            seq = seq.rename({0:"seq"}, axis='columns')
            header = header.rename({0:"header"}, axis='columns')
            fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
            fasta_df = fasta_df.rename({0: 'header', 1: 'seq'}, axis='columns')
            seq = fasta_df['seq'].apply(lambda x: pd.Series(list(x)))
            fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
            fasta_df = fasta_df.set_index(0)
        else:
            fasta_df = self.df_alignments
        
        # Clean Unicode characters in the FASTA DataFrame index to match the main script processing
        import unicodedata
        
        def clean_unicode(text):
            # Normalize Unicode and convert to ASCII
            normalized = unicodedata.normalize('NFKD', str(text))
            ascii_text = normalized.encode('ascii', 'ignore').decode('ascii')
            return ascii_text
        
        # Apply the same Unicode cleaning to FASTA sample names that was applied in main script
        fasta_df.index = fasta_df.index.map(clean_unicode)
        
        with open(self.tree, 'rt') as tree_file:
            for line in tree_file:
                line = re.sub('[:,]', '\n', line)
                line = re.sub('[)(]', '', line)
                line = re.sub('[0-9].*\.[0-9].*\n', '', line)
                line = re.sub("'", '', line)
                line = re.sub('root\n', '', line)
        
        sample_order = line.split('\n')
        sample_order = list(filter(None, sample_order))
        sample_order.insert(0, 'root')
        
        # Clean Unicode characters in sample_order to match the processing done in main script
        cleaned_sample_order = []
        for sample in sample_order:
            if sample == 'root':
                cleaned_sample_order.append(sample)
            else:
                cleaned_sample_order.append(clean_unicode(sample))
        
        sample_order = cleaned_sample_order
        
        # DEBUG: Print available samples and requested samples
        if self.debug:
            print("[{}] Available samples in fasta_df:".format(self.table_name))
            print(sorted(fasta_df.index.tolist()))
            print("[{}] Requested samples in sample_order:".format(self.table_name))
            print(sorted(sample_order))
        
        # Filter sample_order to only include samples that exist in fasta_df
        available_samples = set(fasta_df.index)
        filtered_sample_order = [sample for sample in sample_order if sample in available_samples]
        
        # Check for missing samples
        missing_samples = set(sample_order) - available_samples
        if missing_samples:
            print("Warning [{}]: The following samples from the tree are not found in the FASTA data: {}".format(self.table_name, missing_samples))
            print("This may be due to Unicode character differences or name mismatches.")
            
            # Try to find close matches for missing samples
            for missing_sample in missing_samples:
                # Look for samples that might be close matches
                possible_matches = []
                for available_sample in available_samples:
                    # Check if the missing sample is a substring of available sample or vice versa
                    if missing_sample.lower() in available_sample.lower() or available_sample.lower() in missing_sample.lower():
                        possible_matches.append(available_sample)
                
                if possible_matches:
                    print("  [{}] Possible matches for '{}': {}".format(self.table_name, missing_sample, possible_matches))
                    # Use the first possible match
                    best_match = possible_matches[0]
                    filtered_sample_order.append(best_match)
                    print("  [{}] Using '{}' as substitute for '{}'".format(self.table_name, best_match, missing_sample))
        
        # Remove duplicates while preserving order
        seen = set()
        final_sample_order = []
        for sample in filtered_sample_order:
            if sample not in seen:
                final_sample_order.append(sample)
                seen.add(sample)
        
        if self.debug:
            print("[{}] Final sample order: {}".format(self.table_name, final_sample_order))
        
        # Use the filtered sample order
        tree_order = fasta_df.loc[final_sample_order]
        tree_order2 = fasta_df.loc[final_sample_order]
        
        # Continue with the rest of the existing method...
        self.write_out_table(tree_order, 'sorted')
        
        ## Sort bias to total number of SNPs per column
        # count number of SNPs in each column
        snp_per_column = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            for element in column:
                if element != column[0] and element != '-':
                    count = count + 1
            snp_per_column.append(count)
        row1 = pd.Series(snp_per_column, tree_order.columns, name="snp_per_column")

        # get the snp count per column
        # for each column in the table
        snp_from_top = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            # skip the first element
            for element in column[1:]:
                if element == column[0] or element == '-':
                    count = count + 1
                else:
                    break
            snp_from_top.append(count)
        row2 = pd.Series(snp_from_top, tree_order.columns, name="snp_from_top")
        tree_order = pd.concat([tree_order, pd.DataFrame([row1])])
        tree_order = pd.concat([tree_order, pd.DataFrame([row2])])
        tree_order = tree_order.T
        tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order = tree_order.T

        # remove snp_per_column and snp_from_top rows
        cascade_order_df = tree_order[:-2]
        self.write_out_table(cascade_order_df, 'cascade1')
        del column
        del row1
        del row2
        del snp_from_top
        del cascade_order_df
        
        ## Start 2nd cascading table: sort bias to longest continues vertical SNP count per column
        row1 = pd.Series(snp_per_column, tree_order2.columns, name="snp_per_column")
        # get the snp count per column
        # for each column in the table
        snp_from_top = []
        for column_header in tree_order2:
            count = 0
            column = tree_order2[column_header]
            index_list_of_ref_differences=[]
            for ind, list_item in enumerate(column[1:].to_list()):
                if list_item != column[0] and list_item != '-':
                    index_list_of_ref_differences.append(ind)
            if index_list_of_ref_differences:  # Check if list is not empty
                c = itertools.count()
                val = max((list(g) for _, g in itertools.groupby(index_list_of_ref_differences, lambda x: x-next(c))), key=len)
                snp_from_top.append(val[0])
            else:
                snp_from_top.append(0)  # Default value if no differences found
        row2 = pd.Series(snp_from_top, tree_order2.columns, name="snp_from_top")
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row1])])
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row2])])
        tree_order2 = tree_order2.T
        tree_order2 = tree_order2.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order2 = tree_order2.T

        # remove snp_per_column and snp_from_top rows
        cascade_order_df = tree_order2[:-2]
        self.write_out_table(cascade_order_df, 'cascade2')


    ###Break and write out table
    # In Tables class, modify write_out_table method:
    def write_out_table(self, df, table_type=None):
        if hasattr(self.mq, 'abs_pos'):
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.mq)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T

        if hasattr(self.gbk, 'abs_pos') and not self.gbk.empty:
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.gbk)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T
        else:
            df = pd.concat([df, pd.Series(name='no annotations').to_frame().T])

        # Add coverage data (only once)
        if self.sample_coverage_dict:
            coverage_data = {
                sample: self.sample_coverage_dict.get(sample, 0.0)
                for sample in df.index 
                if sample not in ['root', 'no annotations', 'Average Coverage']
            }
            if coverage_data:  # Only add if we have coverage data
                df['Average Coverage'] = pd.Series(coverage_data)

        if self.show_groups:
            # Add grouping information
            joined_data = {key: '; '.join(map(str, value)) for key, value in self.group_vcfs_dict.items()}
            new_series = pd.Series(joined_data)
            if 'Grouping' not in df.columns:
                df.insert(0, 'Grouping', new_series)

        # Continue with existing chunking code...
        max_size=10000
        count=0
        chunk_start=0
        chunck_end=0
        column_count = df.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                chunck_end += max_size
                df_split = df.iloc[:, chunk_start:chunck_end]
                if 'Grouping' not in df.columns and self.show_groups:
                    df_split.insert(0, 'Grouping', new_series)
                df_split.to_json(os.path.join(self.write_path, 'df{}.json'.format(count)), orient='split')
                self.excel_formatter(os.path.join(self.write_path, 'df{}.json'.format(count)), 
                                os.path.join(self.write_path, '{}_{}table{}-{}.xlsx'.format(self.table_name, table_type + '_', count, self.st)))
                os.remove(os.path.join(self.write_path, 'df{}.json'.format(count)))
                chunk_start += max_size
                column_count -= max_size
            count += 1
            df_split = df.iloc[:, chunk_start:]
            if 'Grouping' not in df.columns and self.show_groups:
                df_split.insert(0, 'Grouping', new_series)
            df_split.to_json(os.path.join(self.write_path, 'df{}.json'.format(count)), orient='split')
            self.excel_formatter(os.path.join(self.write_path, 'df{}.json'.format(count)), 
                            os.path.join(self.write_path, '{}_{}table{}-{}.xlsx'.format(self.table_name, table_type + '_', count, self.st)))
            os.remove(os.path.join(self.write_path, 'df{}.json'.format(count)))
            self.multiple_excel_files = True
        else:
            if 'Grouping' not in df.columns and self.show_groups:
                df.insert(0, 'Grouping', new_series)
            df.to_json(os.path.join(self.write_path, 'df.json'), orient='split')
            self.excel_formatter(os.path.join(self.write_path, 'df.json'), 
                            os.path.join(self.write_path, '{}_{}_table-{}.xlsx'.format(self.table_name, table_type, self.st)))
            os.remove(os.path.join(self.write_path, 'df.json'))
            self.multiple_excel_files = False
            if table_type == 'cascade1':
                self.table_to_tree_path = os.path.join(self.write_path, '{}_{}_table-{}.xlsx'.format(self.table_name, table_type, self.st))

    def excel_formatter(self, df_json, write_to, group=None):
        """Format Excel output with proper styling and conditional formatting."""
        import pandas.io.formats.excel
        pandas.io.formats.excel.header_style = None
        st = self.st
        table_df = pd.read_json(df_json, orient='split')
        writer = pd.ExcelWriter(write_to, engine='xlsxwriter')
        writer.book.use_zip64()
        table_df.to_excel(writer, sheet_name='Sheet1')
        wb = writer.book
        ws = writer.sheets['Sheet1']

        # Define formats
        formatA = wb.add_format({'bg_color': '#58FA82'})
        formatG = wb.add_format({'bg_color': '#F7FE2E'})
        formatC = wb.add_format({'bg_color': '#0000FF'})
        formatT = wb.add_format({'bg_color': '#FF0000'})
        formatnormal = wb.add_format({'bg_color': '#FDFEFE'})
        formatlowqual = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formathighqual = wb.add_format({'font_color': '#000000', 'bg_color': '#FDFEFE'})
        formatambigous = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formatN = wb.add_format({'bg_color': '#E2CFDD'})
        format_rotation = wb.add_format({})
        format_rotation.set_rotation(90)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        
        # Format for integer display in coverage column
        formatcoverage = wb.add_format({
            'num_format': '0',  # Integer format
            'align': 'right'
        })

        rows, cols = table_df.shape
        end_col = cols

        # Set starting column based on whether groups are shown
        if self.show_groups:
            start_col = 2
        else:
            start_col = 1

        # Basic column formatting
        ws.set_column(0, 0, 30)  # Sample name column
        ws.set_column(1, cols, 2.1)  # Position columns

        # Get Coverage column position if present
        coverage_col = None
        if 'Average Coverage' in table_df.columns:
            coverage_col = table_df.columns.get_loc('Average Coverage') + 1
            ws.set_column(coverage_col, coverage_col, 4.2)  # Double width for coverage column

        # Freeze panes
        ws.freeze_panes(2, start_col)

        # Set annotation row height
        ws.set_row(rows + 1, cols + 1, formatannotation)

        # Apply quality-based formatting (excluding Average Coverage column)
        format_end_col = coverage_col if coverage_col is not None else end_col
        ws.conditional_format(rows - 2, start_col, rows - 1, format_end_col - 1,
                            {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})

        # Apply reference-based formatting (excluding Average Coverage column)
        if self.show_groups:
            ws.conditional_format(2, start_col, rows - 2, format_end_col - 1,
                                {'type': 'cell', 'criteria': '==', 'value': 'C$2', 'format': formatnormal})
        else:
            ws.conditional_format(2, start_col, rows - 2, format_end_col - 1,
                                {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': formatnormal})

        # Apply nucleotide-based formatting
        nucleotides = [
            ('A', formatA), ('G', formatG), ('C', formatC), ('T', formatT),
            ('S', formatambigous), ('Y', formatambigous), ('R', formatambigous),
            ('W', formatambigous), ('K', formatambigous), ('M', formatambigous),
            ('N', formatN), ('-', formatN)
        ]
        
        for nuc, format_obj in nucleotides:
            ws.conditional_format(2, start_col, rows - 2, format_end_col - 1,
                                {'type': 'text',
                                'criteria': 'containing',
                                'value': nuc,
                                'format': format_obj})

        # Write column headers with rotation
        for col_num, col_name in enumerate(table_df.columns):
            ws.write(0, col_num + 1, col_name, format_rotation)

        # Handle Average Coverage column
        coverage_col_name = 'Average Coverage'
        if coverage_col_name in table_df.columns:
            coverage_col = table_df.columns.get_loc(coverage_col_name) + 1
            coverage_data = table_df[coverage_col_name]

            # First, write blank cells for all special rows and the annotation row
            for row_idx in range(1, rows + 2):  # Include both regular rows and annotation row
                try:
                    if row_idx == rows + 1:  # This is the annotation row
                        ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))
                    else:
                        row_name = table_df.index[row_idx-1]
                        if row_name in ['root', 'MQ', 'no annotations']:
                            ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))
                except IndexError:
                    ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))

            # Then write coverage values for regular rows
            for row_idx in range(1, rows):
                row_name = table_df.index[row_idx-1]
                if row_name not in ['root', 'MQ', 'no annotations']:
                    try:
                        value = coverage_data.iloc[row_idx-1]
                        if pd.notna(value):
                            # Convert to integer and write
                            int_value = int(round(float(value)))
                            ws.write(row_idx, coverage_col, int_value, formatcoverage)
                        else:
                            ws.write(row_idx, coverage_col, 'N/A', formatcoverage)
                    except (IndexError, ValueError, TypeError):
                        ws.write(row_idx, coverage_col, 'N/A', formatcoverage)

        # Set annotation row
        ws.set_row(rows, 400, formatannotation)
        
        writer.close()

class Hash_Names:

    def __init__(self, fasta_alignments=None, debug=False,):
        self.fasta_alignments = fasta_alignments
        self.debug = debug

    def hash_fasta(self,):
        unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
        hashed_fasta = self.fasta_alignments.replace('.fasta', '_hashed.fasta')
        self.hashed_fasta = hashed_fasta
        checksum_dict = {}
        record_iterator = SeqIO.parse(self.fasta_alignments, "fasta")
        outfasta = open(hashed_fasta , 'at')
        for fasta_file in record_iterator:
            if fasta_file.description == "root":
                print('>{}\n{}'.format(fasta_file.description, fasta_file.seq), file=outfasta)
                checksum_dict.update({fasta_file.description:fasta_file.description})
            else:
                unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
                if fasta_file.description in checksum_dict.values():
                    dup_header = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(2))
                    checksum_dict.update({unique_number:fasta_file.description + "-DUPLICATE_HEADER_NAME-" + dup_header})
                else:
                    checksum_dict.update({unique_number:fasta_file.description})
                print ('>{}\n{}'.format(unique_number, fasta_file.seq), file=outfasta)
        outfasta.close()

        idtable = open("idtable.txt" , 'wt')
        for key, value in checksum_dict.items():
            print('{}\t{}'.format(key, value), file=idtable)
        idtable.close()
        self.idtable = "idtable.txt"
        return hashed_fasta

    def hash_tree_to_original_name(self, tree):
        # read entire tree into variable as string obj
        with open(tree, 'rt') as open_tree:
            entire_file = open_tree.read()
        with open(self.idtable , 'rt') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split("\t")
                org = str(line[1])
                org = org.rstrip()
                check = str(line[0])
                check = check.rstrip()
                entire_file = re.sub(check, "'" + org + "'", entire_file)
            f.close()
        original_tree_names = tree.replace('.tre', '_original_names.tre')
        write_out = open(original_tree_names , 'wt')
        write_out.write(entire_file)
        write_out.close()
        return original_tree_names


class Parsimonious:
    def __init__(self, fasta_alignments=None, debug=False,):

        with open('reformated.fasta', 'w') as reformat:
            sequence = SeqIO.parse(fasta_alignments, "fasta")
            for each in sequence:
                print('>{}\n{}'.format(each.description, each.seq), file=reformat)

        fasta_df = pd.read_csv('reformated.fasta', header=None, sep='^')
        seq = fasta_df.iloc[1::2].reset_index(drop=True)
        header = fasta_df.iloc[0::2].reset_index(drop=True).replace(to_replace=r'>', value='', regex=True)
        seq = seq.rename({0:"seq"}, axis='columns')
        header = header.rename({0:"header"}, axis='columns')
        fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
        fasta_df = fasta_df.rename({0: 'header', 1: 'seq'}, axis='columns')
        seq = fasta_df['seq'].apply(lambda x: pd.Series(list(x)))
        fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
        fasta_df = fasta_df.set_index(0)
        try:
            ref_series = fasta_df.loc['root']
            fasta_df = fasta_df.drop(['root']) #in all_vcf reference_seq needs to be removed
        except KeyError:
            print('Check that there is a "root" named')
            sys.exit(0)
        parsimony = fasta_df.loc[:, (fasta_df != fasta_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = fasta_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        parsimonious_fasta = fasta_alignments.replace('.fasta', '_parsimonious.fasta')
        with open(parsimonious_fasta, 'w') as write_out:
            for index, row in out_df.iterrows():
                print('>{}'.format(row.name), file=write_out)
                for pos in row:
                    print(pos, end='', file=write_out)
                print("", file=write_out)
        self.parsimonious_fasta = parsimonious_fasta

        if not debug:
            os.remove('reformated.fasta')


if __name__ == '__main__':
    # Set multiprocessing start method for Python 3.12 compatibility
    multiprocessing.set_start_method('spawn', force=True)

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
        ---------------------------------------------------------
        Usage:
        vsnp3_fasta_to_snps_table.py -f *fasta -pn
        ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True, help='Provide an alignment file in FASTA format')
    parser.add_argument('-p', '--parsimonious', action='store_true', dest='parsimonious', help='Only keep parsimonious SNPs from FASTA alignment file.  This is different than the uninformative SNPs removed via vSNP pipeline.  This is to be used when just working with an aligned FASTA file.')
    parser.add_argument('--show_groups', action='store_true', dest='show_groups', help='Show group names in SNP table')
    parser.add_argument('-n', '--hash_names', action='store_true', dest='hash_names', help='Hash FASTA names to rid of any RAxML illegal characters')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling')
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.abspath(__file__), __version__))
    args = parser.parse_args()
    print('\n{} SET ARGUMENTS:'.format(os.path.basename(__file__)))
    print(args)

    initial_fasta = args.fasta

    if args.parsimonious:
        parsimonious = Parsimonious(args.fasta)
        args.fasta = parsimonious.parsimonious_fasta

    if args.hash_names:
        hash_names = Hash_Names(args.fasta)
        args.fasta = hash_names.hash_fasta()

    tree = Tree(fasta_alignments=args.fasta, debug=args.debug)

    if args.hash_names:
        tree.newick = hash_names.hash_tree_to_original_name(tree.newick)
    
    if args.parsimonious:
        fasta_alignments=parsimonious.parsimonious_fasta
    else:
        fasta_alignments=initial_fasta
    tables = Tables(fasta_alignments=fasta_alignments, tree=tree.newick, debug=False)
    tables.build_tables()