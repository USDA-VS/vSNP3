#!/usr/bin/env python

__version__ = "3.18"

import os
import subprocess
import sys
import re
import glob
import time
import random
from datetime import datetime
import pandas as pd
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from concurrent import futures
import argparse
import textwrap
import itertools
from Bio import SeqIO


class Tree:
    ''' 
    '''

    def __init__(self, fasta_alignments=None, write_path=None, tree_name=None, debug=False,):
        # Find an optimal compiled version of RAxML in conda
        try:
            subprocess.call("raxml", stdout=open(os.devnull, 'wb'))
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
            # print(f'set RAxML to {raxml}')
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

        raxml_version = subprocess.check_output(f'{raxml} -v', shell=True, text=True)
        self.raxml_version = raxml_version.split('\n')[2]

        os.system(f'{raxml} -s {fasta_alignments} -n raxml -m GTRCATI -o root -w {write_path} -p 456123 -T 4 > /dev/null 2>&1') #> /dev/null 2>&1
        try:
            newick = f'{write_path}/{tree_name}_{self.st}.tre'
            os.rename(f'{write_path}/RAxML_bestTree.raxml', newick)
            raxml_to_remove = glob.glob(f'{write_path}/RAxML*')
            for each in raxml_to_remove:
                os.remove(each)
            try:
                reduced_file = glob.glob(f'{write_path}/*.reduced')
                os.remove(reduced_file[0])
            except (FileNotFoundError, IndexError) as e:
                pass
        except FileNotFoundError:
            with open(f'{write_path}/SEE_RAXML_INFO', 'w') as message_out:
                print(f'check sample numbers', file=message_out)
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
                    print(f'>{row.name}', file=write_out)
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
        # print(f'in_df size: {in_df.shape}')
        parsimony = in_df.loc[:, (in_df != in_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = in_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        return out_df


class Tables:

    def __init__(self, fasta_alignments=None, df_alignments=None, tree=None, gbk=None, mq=None, write_path=None, table_name=None, debug=False,):
        self.fasta_alignments = fasta_alignments
        self.df_alignments = df_alignments
        self.tree = tree
        self.gbk = gbk
        self.mq = mq
        self.debug = debug
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')

        if not write_path:
            write_path = os.getcwd()
        if not table_name:
            table_name = "generic"
        self.write_path = write_path
        self.table_name = table_name

    def build_tables(self, ):
        if self.df_alignments is None:
            reformated = f'{self.write_path}/reformated.fasta'
            with open(reformated, 'w') as reformat:
                # try:
                sequence = SeqIO.parse(self.fasta_alignments, "fasta")
                # except AttributeError as e:
                #     print(f'\n\t#####\n\t##### {e}, Reformated FASTA {reformated} \n\t#####\n')
                for each in sequence:
                    print(f'>{each.description}\n{each.seq}', file=reformat)
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
        with open(self.tree, 'rt') as tree_file: #must be the single line newick format.  Not Nexus which will be mutliline often with formating
            for line in tree_file:
                line = re.sub('[:,]', '\n', line)
                line = re.sub('[)(]', '', line)
                line = re.sub('[0-9].*\.[0-9].*\n', '', line)
                line = re.sub("'", '', line)
                line = re.sub('root\n', '', line)
        sample_order = line.split('\n')
        sample_order = list(filter(None, sample_order))
        sample_order.insert(0, 'root')
        tree_order = fasta_df.loc[sample_order] #cascading1 table
        tree_order2 = fasta_df.loc[sample_order] #cascading2 table
        self.write_out_table(tree_order, 'sorted') #table_type, sorted or cascade, type is labeled on the output Excel file
        
        ## Sort bias to total number of SNPs per column
        # count number of SNPs in each column
        snp_per_column = []
        for column_header in tree_order:
            count = 0
            column = tree_order[column_header]
            # for each element in the column
            for element in column:
                if element != column[0] and element != '-': #column[0] is top row/root/reference, element is everything below it.
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
                if element == column[0] or element == '-': # when - keep count, essentially skip -
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
        self.write_out_table(cascade_order_df, 'cascade1') #table_type, sorted or cascade, type is labeled on the output Excel file
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
            c = itertools.count()
            val = max((list(g) for _, g in itertools.groupby(index_list_of_ref_differences, lambda x: x-next(c))), key=len)
            snp_from_top.append(val[0]) #starting row number with longest continous SNPs in column
        row2 = pd.Series(snp_from_top, tree_order2.columns, name="snp_from_top")
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row1])])
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row2])])
        tree_order2 = tree_order2.T
        tree_order2 = tree_order2.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order2 = tree_order2.T

        # remove snp_per_column and snp_from_top rows
        cascade_order_df = tree_order2[:-2]
        self.write_out_table(cascade_order_df, 'cascade2') #table_type, sorted or cascade, type is labeled on the output Excel file

    ###Break and write out table
    def write_out_table(self, df, table_type=None):
        if hasattr(self.mq, 'abs_pos'):
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.mq)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T
        # else:
            # df = df.append(pd.Series(name='no map quality'))
        if hasattr(self.gbk, 'abs_pos') and not self.gbk.empty:
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.gbk)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T
        else:
            # df = df.append(pd.Series(name='no annotations'))
            df = pd.concat([df, pd.Series(name='no annotations').to_frame().T])
        max_size=10000 #max columns allowed in tables
        count=0
        chunk_start=0
        chunck_end=0
        column_count = df.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                # print(f'{column_count} columns > {max_size}, cascade table break {count}')
                chunck_end += max_size
                df_split = df.iloc[:, chunk_start:chunck_end]
                df_split.to_json(f'{self.write_path}/df{count}.json', orient='split')
                self.excel_formatter(f'{self.write_path}/df{count}.json', f'{self.write_path}/{self.table_name}_{table_type}_table{count}-{self.st}.xlsx')
                os.remove(f'{self.write_path}/df{count}.json')
                chunk_start += max_size
                column_count -= max_size
            count += 1
            # print(f'Last break {column_count} columns, cascade table break {count}')
            df_split = df.iloc[:, chunk_start:]
            df_split.to_json(f'{self.write_path}/df{count}.json', orient='split')
            self.excel_formatter(f'{self.write_path}/df{count}.json', f'{self.write_path}/{self.table_name}_{table_type}_table{count}-{self.st}.xlsx')
            os.remove(f'{self.write_path}/df{count}.json')
        else: # no break needed
            df.to_json(f'{self.write_path}/df.json', orient='split')
            self.excel_formatter(f'{self.write_path}/df.json', f'{self.write_path}/{self.table_name}_{table_type}_table-{self.st}.xlsx')
            os.remove(f'{self.write_path}/df.json')

    def excel_formatter(self, df_json, write_to, group=None):
        import pandas.io.formats.excel
        pandas.io.formats.excel.header_style = None
        # sample_path_name = self.sample_path_name
        st = self.st
        table_df = pd.read_json(df_json, orient='split')
        writer = pd.ExcelWriter(write_to, engine='xlsxwriter')
        table_df.to_excel(writer, sheet_name='Sheet1')
        wb = writer.book
        ws = writer.sheets['Sheet1']
        formatA = wb.add_format({'bg_color': '#58FA82'})
        formatG = wb.add_format({'bg_color': '#F7FE2E'})
        formatC = wb.add_format({'bg_color': '#0000FF'})
        formatT = wb.add_format({'bg_color': '#FF0000'})
        formatnormal = wb.add_format({'bg_color': '#FDFEFE'})
        formatlowqual = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formathighqual = wb.add_format({'font_color': '#000000', 'bg_color': '#FDFEFE'})
        formatambigous = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formatN = wb.add_format({'bg_color': '#E2CFDD'})
        rows, cols = table_df.shape

        ws.set_column(0, 0, 30)
        ws.set_column(1, cols, 2.1)
        ws.freeze_panes(2, 1)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        #set last row
        ws.set_row(rows + 1, cols + 1, formatannotation)

        #'first_row', 'first_col', 'last_row', and 'last_col'
        # Careful that row/column locations don't overlap
        ws.conditional_format(rows - 2, 1, rows - 1, cols, {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': formatnormal})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'A', 'format': formatA})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'G', 'format': formatG})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'C', 'format': formatC})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'T', 'format': formatT})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'S', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'Y', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'R', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'W', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'K', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'M', 'format': formatambigous})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': 'N', 'format': formatN})
        ws.conditional_format(2, 1, rows - 2, cols, {'type': 'text', 'criteria': 'containing', 'value': '-', 'format': formatN})

        format_rotation = wb.add_format({})
        format_rotation.set_rotation(90)
        # ws.set_row(0, None, format_rotation)
        for columnnum, columnname in enumerate(list(table_df.columns)):
            ws.write(0, columnnum + 1, columnname, format_rotation)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        #set last row
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
                print(f'>{fasta_file.description}\n{fasta_file.seq}', file=outfasta)
                checksum_dict.update({fasta_file.description:fasta_file.description})
            else:
                unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
                if fasta_file.description in checksum_dict.values():
                    dup_header = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(2))
                    checksum_dict.update({unique_number:fasta_file.description + "-DUPLICATE_HEADER_NAME-" + dup_header})
                else:
                    checksum_dict.update({unique_number:fasta_file.description})
                print (f'>{unique_number}\n{fasta_file.seq}', file=outfasta)
        outfasta.close()

        idtable = open("idtable.txt" , 'wt')
        for key, value in checksum_dict.items():
            print(f'{key}\t{value}', file=idtable)
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
        self.original_tree_names = original_tree_names
        return original_tree_names

        if not self.debug:
            os.remove(self.idtable)
            os.remove(self.hashed_fasta)


class Parsimonious:
    def __init__(self, fasta_alignments=None, debug=False,):

        with open('reformated.fasta', 'w') as reformat:
            sequence = SeqIO.parse(fasta_alignments, "fasta")
            for each in sequence:
                print(f'>{each.description}\n{each.seq}', file=reformat)

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
                print(f'>{row.name}', file=write_out)
                for pos in row:
                    print(pos, end='', file=write_out)
                print("", file=write_out)
        self.parsimonious_fasta = parsimonious_fasta

        if not debug:
            os.remove('reformated.fasta')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
        ---------------------------------------------------------
        Usage:
        vsnp3_fasta_to_snps_table.py -f *fasta -pn
        ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True, help='Provide an alignment file in FASTA format')
    parser.add_argument('-p', '--parsimonious', action='store_true', dest='parsimonious', help='Only keep parsimonious SNPs from FASTA alignment file.  This is different than the uninformative SNPs removed via vSNP pipeline.  This is to be used when just working with an aligned FASTA file.')
    parser.add_argument('-n', '--hash_names', action='store_true', dest='hash_names', help='Hash FASTA names to rid of any RAxML illegal characters')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')
    args = parser.parse_args()
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
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