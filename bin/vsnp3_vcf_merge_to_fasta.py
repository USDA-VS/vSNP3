#!/usr/bin/env python

__version__ = "3.28"

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
from collections.abc import Mapping  # For Python 3.12 compatibility

from vsnp3_file_setup import Setup


class Merge_VCF(Setup):
    ''' 
    Merges VCF file SNPs into a FASTA reference sequence
    '''

    def __init__(self, vcf=None, fasta=None, qual_threshold=300, mq_threshold=56, ac_threshold=1, indels=False):
        '''
        Initialize the Merge_VCF class
        
        Parameters:
        -----------
        vcf : str
            Path to the VCF file
        fasta : str
            Path to the FASTA file
        qual_threshold : int
            Quality threshold for filtering SNPs
        mq_threshold : int
            Mapping quality threshold for filtering SNPs
        ac_threshold : int
            Allele count threshold for filtering SNPs
        indels : bool
            Whether to include indels in the output
        '''
        Setup.__init__(self, FASTA=fasta)
        self.fasta = fasta
        self.fasta_name = re.sub('[.].*', '', os.path.basename(fasta))
        self.vcf = vcf
        self.vcf_name = re.sub('[.].*', '', os.path.basename(vcf)).strip('_zc')
        self.output_name = f'{self.vcf_name}_vcf_merged_{self.fasta_name}'
        self.qual_threshold = qual_threshold
        self.mq_threshold = mq_threshold
        self.ac_threshold = ac_threshold
        self.indels = indels

    def read_vcf(self, path):
        '''
        Read a VCF file and return a pandas DataFrame
        
        Parameters:
        -----------
        path : str
            Path to the VCF file
            
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing the VCF data
        '''
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
        
        df = pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})
        
        # Convert columns to appropriate types
        df['POS'] = pd.to_numeric(df['POS'], errors='coerce').fillna(0).astype(int)
        df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce').fillna(0).astype(int)
        
        # Helper function to parse INFO field
        def parse_info(info_str):
            info_dict = {}
            for item in info_str.split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
            return info_dict
        
        # Extract INFO fields
        df['INFO_DICT'] = df['INFO'].apply(parse_info)
        df['AC'] = df['INFO_DICT'].apply(lambda x: x.get('AC', None))
        df['DP'] = df['INFO_DICT'].apply(lambda x: x.get('DP', None))
        df['MQ'] = df['INFO_DICT'].apply(lambda x: x.get('MQ', None))
        
        # Convert to numeric
        df['AC'] = pd.to_numeric(df['AC'], errors='coerce').fillna(0).astype(int)
        df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype(int)
        df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce').fillna(0).astype(int)
        
        # Drop unnecessary columns
        df = df.drop(columns=['INFO', 'INFO_DICT', 'ID', 'FILTER', 'FORMAT'])

        return df

    def run(self):
        '''
        Run the VCF to FASTA merging process
        '''
        # Filter VCF to only include entries with single REF and ALT values
        os.system('awk \'BEGIN { OFS=FS="\t" } $4 !~ /,/\' ' + self.vcf + ' | awk \'BEGIN { OFS=FS="\t" } $5 !~ /,/\' > temp_file.vcf')
        
        # Read and filter VCF
        df = self.read_vcf('temp_file.vcf')
        
        # Apply filters based on parameters
        if self.indels:
            df = df[(df['QUAL'] > self.qual_threshold) & 
                   (df['MQ'] >= self.mq_threshold) & 
                   (df['AC'] >= self.ac_threshold)]
        else:
            df['REF'] = df['REF'].astype('str')
            df['ALT'] = df['ALT'].astype('str')
            df = df[(df['QUAL'] > self.qual_threshold) & 
                   (df['MQ'] >= self.mq_threshold) & 
                   (df['AC'] >= self.ac_threshold) & 
                   (df['REF'].str.len() == 1) & 
                   (df['ALT'].str.len() == 1)]
        
        # Keep only position and alt allele
        df = df[['POS', 'ALT']]
        
        # Read reference sequence
        sequence = SeqIO.read(self.fasta, "fasta")
        
        # Create dictionary of positions to alt alleles
        dfdict = dict(zip(df.POS, df.ALT))
        
        # Convert to 0-indexed positions
        vcf_dict_zero_index = {key-1: value for key, value in dfdict.items()}
        
        # Convert sequence to list and then to dictionary
        sequence_list = list(sequence.seq)
        seq_dict = dict(enumerate(sequence_list))
        
        # Merge dictionaries, with VCF variants overriding reference
        merge_dicts = {**seq_dict, **vcf_dict_zero_index}
        
        # Convert back to sequence
        updated_seq = "".join(list(merge_dicts.values()))

        # Create new sequence record
        record = SeqRecord(
            Seq(updated_seq),
            id=self.output_name,
            name=self.output_name,
            description=""
        )

        # Write to file
        with open(f'{self.output_name}.fasta', "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")
        
        # Clean up temporary file
        os.remove('temp_file.vcf')


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Merge VCF SNPs into a FASTA reference sequence
    This script takes a VCF file and the FASTA file that was used
    to generate it, and creates a new FASTA file with the variants
    from the VCF incorporated into the reference sequence.

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

    merge_vcf = Merge_VCF(args.vcf, args.fasta, int(args.qual_threshold), int(args.mq_threshold), int(args.ac_threshold), args.indels)
    merge_vcf.run()


# Created 2022 by Tod Stuber