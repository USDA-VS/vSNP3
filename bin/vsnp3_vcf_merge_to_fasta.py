#!/usr/bin/env python

__version__ = "0.0.1"

import os
import re
import io
import argparse
import textwrap
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from vsnp3_file_setup import Setup


class Merge_VCF(Setup):
    ''' 
    '''

    def __init__(self, vcf=None, fasta=None, qual_threshold=300, mq_threshold=56, ac_threshold=1, indels=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, FASTA=fasta)
        self.fasta = fasta
        self.fasta_name = re.sub('[.].*', '', os.path.basename(fasta))
        self.vcf = vcf
        self.vcf_name = re.sub('[.].*', '', os.path.basename(vcf)).strip('_zc')
        self.output_name = f'{self.vcf_name}_vcf_merged_{self.fasta_name}'
        self.fasta = fasta
        self.qual_threshold = qual_threshold
        self.mq_threshold = mq_threshold
        self.ac_threshold = ac_threshold
        self.indels = indels

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
        df['AC'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x.split(";") if "=" in item).get('AC', None))
        df['DP'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x.split(";") if "=" in item).get('DP', None))
        df['MQ'] = df['INFO'].apply(lambda x: dict(item.split("=") for item in x.split(";") if "=" in item).get('MQ', None))
        df['AC'] = pd.to_numeric(df['AC'], errors='coerce').fillna(0).astype(int)
        df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype(int)
        df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce').fillna(0).astype(int)
        df = df.drop(columns=['INFO', 'ID', 'FILTER', 'FORMAT'])

        return df

    def run(self,):
        '''
        description
        '''
        os.system('awk \'BEGIN { OFS=FS="\t" } $4 !~ /,/\' ' + self.vcf + ' | awk \'BEGIN { OFS=FS="\t" } $5 !~ /,/\' > temp_file.vcf')
        # single_df = single_df[(single_df['QUAL'] > self.qual_threshold) & (single_df['AC'] == 2) & (single_df['REF'].str.len() == 1) & (single_df['ALT'].str.len() == 1) & (single_df['MQ'] >= self.mq_threshold)]
        df = self.read_vcf('temp_file.vcf')
        if self.indels:
            df = df[(df['QUAL'] > self.qual_threshold) & (df['MQ'] >= self.mq_threshold) & (df['AC'] >= self.ac_threshold)]
        else:
            df['REF'] = df['REF'].astype('str')
            df['ALT'] = df['ALT'].astype('str')
            df = df[(df['QUAL'] > self.qual_threshold) & (df['MQ'] >= self.mq_threshold) & (df['AC'] >= self.ac_threshold) & (df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)]
        df = df.drop(['CHROM', 'QUAL', 'REF', 'AC', 'DP', 'MQ'], axis=1)
        sequence = SeqIO.read(self.fasta, "fasta")
        dfdict = dict(zip(df.POS, df.ALT))
        #sequence list will be 0 indexed, so update vcf to be the same
        vcf_dict_zero_index={}
        for key, value in dfdict.items():
            vcf_dict_zero_index[key-1] = value
        sequence_list = list(sequence.seq)
        seq_dict = dict(enumerate(sequence_list))
        merge_dicts = {**seq_dict, **vcf_dict_zero_index}
        updated_seq = "".join(list(merge_dicts.values()))

        record = SeqRecord(
            Seq(updated_seq),
            id = self.output_name,
            name = self.output_name,
            description = ""
        )

        with open(f'{self.output_name}.fasta', "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        os.remove('temp_file.vcf')

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Place description

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-v', '--vcf', action='store', dest='vcf', required=False, help='Required: VCF file will be merged into FASTA file it was made from')
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', help='Required: FASTA used to build the VCF file')

    parser.add_argument('-x', '--qual_threshold', action='store', dest='qual_threshold', default=300, required=False, help='Optional: Minimum QUAL threshold for using a SNP')
    parser.add_argument('-y', '--mq_threshold', action='store', dest='mq_threshold', default=56, required=False, help='Optional: Minimum MQ threshold for using a SNP')
    parser.add_argument('-z', '--ac', action='store', dest='ac_threshold', default=1, required=False, help='Optional: Minimum AC threshold for using a SNP')
    parser.add_argument('-i', '--indels', action='store_true', dest='indels', default=False, help='Optional: Include indels')
    
    parser.add_argument('--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    merge_vcf = Merge_VCF(args.vcf, args.fasta, int(args.qual_threshold), int(args.mq_threshold), int(args.ac_threshold), args.indels,)
    merge_vcf.run()


# Created 2022 by Tod Stuber
