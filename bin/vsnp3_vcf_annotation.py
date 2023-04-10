#!/usr/bin/env python

__version__ = "3.14"

import os
import re
import argparse
import textwrap
import pandas as pd

from vsnp3_annotation import Annotation

class VCF_Annotation():
    ''' 
    '''
    def __init__(self, gbk_list=None, vcf_file=None,):

        annotation = Annotation(gbk_list=gbk_list)
        
        header_out = open('v_header.csv', 'w+')
        with open(vcf_file) as fff:
            for line in fff:
                if re.search('^#', line):
                    print(line.strip(), file=header_out)
        header_out.close()

        vcf_df = pd.read_csv(vcf_file, sep='\t', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], comment='#')
        vcf_df['ABS_POS'] = vcf_df['CHROM'].map(str) + ':' + vcf_df['POS'].map(str)
        annotation_dict={}
        for index, row in vcf_df.iterrows():
            annotation.run(row['ABS_POS'], row['ALT'])
            annotation_dict[row['ABS_POS']] = f'cds_nt_start={annotation.cds_nt_start};cds_nt_end={annotation.cds_nt_end};gene={annotation.gene};product={annotation.product};aa_residue_pos={annotation.aa_residue_pos};snp_nt={annotation.snp_nt};aa_pos={annotation.aa_pos};reference base code={annotation.reference_base_code};snp_base_code={annotation.snp_base_code};ref_aa={annotation.ref_aa};snp_aa={annotation.snp_aa};mutation_type={annotation.mutation_type}'
        vcf_df = vcf_df.set_index('ABS_POS')
        vcf_df.drop(['ID'], axis=1, inplace=True)
        annotation_df = pd.DataFrame.from_dict(annotation_dict, orient='index', columns=["ID"])
        annotation_df.index.name = 'ABS_POS'

        vcf_df = vcf_df.merge(annotation_df, how='left', left_index=True, right_index=True)
        vcf_df = vcf_df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]]
        vcf_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
        cat_files = ['v_header.csv', 'v_annotated_body.csv']
        name = vcf_file.replace('.vcf', '')
        with open(f'{name}_annotated.vcf', "wb") as outfile:
            for cf in cat_files:
                with open(cf, "rb") as infile:
                    outfile.write(infile.read())
        os.remove('v_header.csv')
        os.remove('v_annotated_body.csv')
        self.vcf = f'{os.getcwd()}/{name}_annotated.vcf'

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Annotate VCF file

    Usage:
    vsnp3_vcf_annotation.py -b *.gbk -c *.vcf

    Give wildcard when more than one *.gbk

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-b', '--gbk_list', nargs='*', dest='gbk_list', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-c', '--vcf', action='store', dest='vcf', default=None, required=True, help='Provide VCF file, Return VCF annotated')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    vcf_annotation = VCF_Annotation(gbk_list=args.gbk_list, vcf_file=args.vcf)

# Created 2021 by Tod Stuber
