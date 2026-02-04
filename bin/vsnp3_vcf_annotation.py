#!/usr/bin/env python

__version__ = "3.33"

import os
import re
import argparse
import textwrap
import pandas as pd

from vsnp3_annotation import Annotation


class VCF_Annotation():
    ''' 
    VCF Annotation class for processing and annotating VCF files
    '''
    def __init__(self, gbk_list=None, vcf_file=None,):

        annotation = Annotation(gbk_list=gbk_list)
        
        # Extract header from VCF file
        try:
            with open('v_header.csv', 'w+') as header_out:
                with open(vcf_file) as fff:
                    for line in fff:
                        if re.search('^#', line):
                            print(line.strip(), file=header_out)
        except IOError as e:
            raise IOError(f"Error processing VCF file {vcf_file}: {e}")

        # Read VCF data
        try:
            vcf_df = pd.read_csv(vcf_file, sep='\t', header=None, 
                                 names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", 
                                       "INFO", "FORMAT", "Sample"], 
                                 comment='#')
        except Exception as e:
            raise Exception(f"Error reading VCF file: {e}")
            
        vcf_df['ABS_POS'] = vcf_df['CHROM'].astype(str) + ':' + vcf_df['POS'].astype(str)
        
        # Build annotation dictionary
        annotation_dict = {}
        for index, row in vcf_df.iterrows():
            try:
                annotation.run(row['ABS_POS'], row['ALT'])
                
                # Build annotation string in parts to avoid f-string complexity
                annotation_parts = [
                    f'cds_nt_start={getattr(annotation, "cds_nt_start", "")}',
                    f'cds_nt_end={getattr(annotation, "cds_nt_end", "")}',
                    f'gene={getattr(annotation, "gene", "")}',
                    f'product={getattr(annotation, "product", "")}',
                    f'aa_residue_pos={getattr(annotation, "aa_residue_pos", "")}',
                    f'snp_nt={getattr(annotation, "snp_nt", "")}',
                    f'aa_pos={getattr(annotation, "aa_pos", "")}',
                    f'reference base code={getattr(annotation, "reference_base_code", "")}',
                    f'snp_base_code={getattr(annotation, "snp_base_code", "")}',
                    f'ref_aa={getattr(annotation, "ref_aa", "")}',
                    f'snp_aa={getattr(annotation, "snp_aa", "")}',
                    f'mutation_type={getattr(annotation, "mutation_type", "")}'
                ]
                
                annotation_string = ';'.join(annotation_parts)
                annotation_dict[row['ABS_POS']] = annotation_string
                
            except Exception as e:
                print(f"Warning: Error annotating position {row['ABS_POS']}: {e}")
                # Create empty annotation for failed positions
                empty_parts = [
                    'cds_nt_start=',
                    'cds_nt_end=',
                    'gene=',
                    'product=',
                    'aa_residue_pos=',
                    'snp_nt=',
                    'aa_pos=',
                    'reference base code=',
                    'snp_base_code=',
                    'ref_aa=',
                    'snp_aa=',
                    'mutation_type='
                ]
                annotation_dict[row['ABS_POS']] = ';'.join(empty_parts)
        
        # Merge annotations with VCF data
        vcf_df = vcf_df.set_index('ABS_POS')
        vcf_df.drop(['ID'], axis=1, inplace=True)
        annotation_df = pd.DataFrame.from_dict(annotation_dict, orient='index', columns=["ID"])
        annotation_df.index.name = 'ABS_POS'

        vcf_df = vcf_df.merge(annotation_df, how='left', left_index=True, right_index=True)
        vcf_df = vcf_df[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]]
        
        # Write annotated body
        try:
            vcf_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
        except Exception as e:
            raise Exception(f"Error writing annotated body: {e}")
        
        # Combine header and body files
        cat_files = ['v_header.csv', 'v_annotated_body.csv']
        name = vcf_file.replace('.vcf', '')
        output_filename = f'{name}_annotated.vcf'
        
        try:
            with open(output_filename, "wb") as outfile:
                for cf in cat_files:
                    with open(cf, "rb") as infile:
                        outfile.write(infile.read())
        except IOError as e:
            raise IOError(f"Error combining files: {e}")
                    
        # Clean up temporary files
        try:
            for temp_file in ['v_header.csv', 'v_annotated_body.csv']:
                if os.path.exists(temp_file):
                    os.remove(temp_file)
        except OSError as e:
            print(f"Warning: Could not remove temporary file: {e}")
            
        self.vcf = os.path.join(os.getcwd(), output_filename)


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

    try:
        vcf_annotation = VCF_Annotation(gbk_list=args.gbk_list, vcf_file=args.vcf)
        print(f"Successfully created annotated VCF: {vcf_annotation.vcf}")
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

# Created 2021 by Tod Stuber