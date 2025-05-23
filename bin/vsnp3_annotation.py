#!/usr/bin/env python3

__version__ = "3.30"

import os
import shutil
import sys
import argparse
import textwrap
from Bio.Seq import Seq
from Bio import SeqIO
from collections import ChainMap
from Bio.SeqFeature import FeatureLocation, CompoundLocation


class Annotation():
    ''' 
    '''

    def __init__(self, gbk_list=None,):

        self.aa_code = { \
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', \
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', \
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', \
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', \
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', \
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', \
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', \
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', \
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', \
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', \
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', \
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', \
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', \
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', \
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', \
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

        cwd = os.getcwd()
        gbk_dict_list=[]
        for each in gbk_list:
            try: #IF gbk provided as path variable copy local
                gbk_dict_list.append(SeqIO.to_dict(SeqIO.parse(each, "genbank")))
            except shutil.SameFileError:
                pass
        self.gbk_dict = dict(ChainMap(*gbk_dict_list)) #merge dictionaries        


    def run(self, abs_pos, snp_nt):
        '''
        When having both abs_pos and snp_nt, the snp_nt allows updating the amino acid to the mutant call
        # abs_pos='NC_017250.1:264518', snp_nt='T'
        # abs_pos='NC_017251.1:396173', snp_nt='T'
        '''
        #can both abs_pos and SNP nt be provided, which would eliminate snp_dictionary
        chrom, position = abs_pos.split(':')
        position = int(position)
        self.chrom = chrom
        self.position = position
        self.cds_nt_start = "n/a"
        self.cds_nt_end = "n/a"
        self.aa_residue_pos = "n/a"
        self.mutation_type = "n/a"
        self.snp_nt = snp_nt
        rbc = ""
        self.reference_base_code = "n/a"
        self.snp_base_code = "SNP nt not provided"
        self.ref_aa = "n/a"
        self.snp_aa = "n/a"
        self.gene = "unlisted gene"
        self.product = "No annotated product"
        aa_residue_pos = ""
        self.aa_pos = "n/a"
        self.feature_found = False
        self.direction = "n/a"
        try:
            for feature in self.gbk_dict[chrom].features:
                if "CDS" == feature.type or "tRNA" == feature.type or "rRNA" == feature.type or "repeat_region" == feature.type or "mobile_element" == feature.type or "ncRNA" == feature.type:
                    for part in feature.location.parts:
                        if position in range(part.start, part.end):
                            self.feature_found = True
                            self.direction = "forward gene"
                            self.cds_nt_start = part.start
                            self.cds_nt_end = part.end
                            aa_location = (position - feature.location.start) / 3
                            aa_residue, nt_in_aa = str(aa_location).split('.')
                            try:
                                self.gene = feature.qualifiers["gene"][0]
                            except KeyError:
                                pass
                            try:
                                self.product = feature.qualifiers['product'][0]
                            except KeyError:
                                try:
                                    self.product = feature.qualifiers['locus_tag'][0]
                                except KeyError:
                                    try:
                                        self.product = feature.type + ", " + feature.qualifiers['label'][0]
                                    except KeyError:
                                        self.product = f'{feature.type}, product_unknown'
                            if nt_in_aa == '0':
                                nt_index_aa = 2 #set index
                                right = position
                                left = position - 3
                                self.aa_pos = 3
                            elif nt_in_aa[0] == '3':
                                nt_index_aa = 0 #set index
                                right = position + 2
                                left = position - 1
                                self.aa_pos = 1
                            elif nt_in_aa[0] == '6':
                                nt_index_aa = 1 #set index
                                right = position + 1
                                left = position - 2
                                self.aa_pos = 2
                            else:
                                #error out without exception to quit
                                right = ''
                            self.aa_residue_pos = int((left - part.start) / 3)
                            self.aa_residue_pos = self.aa_residue_pos + 1 # aa is zero based so add 1
                            if part.strand == -1: # Reverse complement
                                self.aa_residue_pos = int((part.end - left) / 3)
                            rbc = self.gbk_dict[chrom].seq[left:right]
                            if part.strand == -1: # Reverse complement
                                rbc = rbc.reverse_complement()
                                nucleotide_seq = Seq(snp_nt) if snp_nt is not None else None
                                snp_nt = str(nucleotide_seq.reverse_complement()) if nucleotide_seq is not None else None
                                self.direction = "reverse gene"
                                block = True
                                if nt_index_aa == 0:
                                    nt_index_aa = 2
                                    block = False
                                elif nt_index_aa == 2 and block == True:
                                    nt_index_aa = 0
                            rbc_list = list(str(rbc))  # Convert Bio.Seq to string first
                            self.reference_base_code = "".join(rbc_list)
                            try:
                                self.ref_aa = self.aa_code[self.reference_base_code]
                            except KeyError:
                                self.ref_aa = 'unfound_ref_AA'
                            #change rbc_list to represent SNP
                            if snp_nt is not None:
                                rbc_list[nt_index_aa] = snp_nt
                                # Example snp_dictionary: SNP at abs pos, {'NC_017250.1:264518': 'T', ...}
                                try:
                                    self.snp_base_code = "".join(rbc_list)
                                    try:
                                        self.snp_aa = self.aa_code[self.snp_base_code]
                                    except KeyError:
                                        self.snp_aa = 'ambiguous'
                                except TypeError:
                                    self.snp_aa = "SNP nt not provided"
                            if self.ref_aa == self.snp_aa:
                                self.mutation_type = "silent mutation"
                            elif self.snp_aa == "ambiguous":
                                self.mutation_type = "unsure-ambiguous"
                            elif self.ref_aa == "unfound_ref_AA":
                                self.mutation_type = ""
                            else:
                                self.mutation_type = "nonsynonymous"
                            return
        except KeyError:
            print(f'\n### KeyError: incorrect chrom in dataset\n### grep -l {chrom} *vcf to find file with error\n### File must be removed\n')
            try:
                print(f'Feature Type: {feature.type}')
            except UnboundLocalError:
                pass
            print(f'Position: {position}\n')
            sys.exit(0)
        return


class VCF_Annotation():
    ''' 
    '''
    def __init__(self, vcf_file=None,):
        self.vcf_file = vcf_file
        

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    
    Used by vSNP to annotate VCF files and SNP tables

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-g', '--gbk', nargs='*', dest='gbk', default=None, required=True, help='Full path to .gbk files.  Wildcard can be used')
    parser.add_argument('-p', '--abs_pos', action='store', dest='abs_pos', default=None, required=True, help='Absolute position, example NC_017250.1:264518')
    parser.add_argument('-s', '--snp_nt', action='store', dest='snp_nt', default=None, required=False, help='Position alt call, example T')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    # print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    # print(args)
    print("\n")

    annotation = Annotation(gbk_list=args.gbk)
    annotation.run(args.abs_pos, args.snp_nt)
    print(args.abs_pos)
    print(f'cds_nt_start: {annotation.cds_nt_start}')
    print(f'cds_nt_end: {annotation.cds_nt_end}')
    print(f'gene: {annotation.gene}')
    print(f'product: {annotation.product}')
    print(f'aa_residue_pos: {annotation.aa_residue_pos}')
    print(f'snp_nt: {annotation.snp_nt}')
    print(f'aa_pos: {annotation.aa_pos}')
    print(f'reference base code: {annotation.reference_base_code}')
    print(f'snp_base_code: {annotation.snp_base_code}')
    print(f'ref_aa: {annotation.ref_aa}')
    print(f'snp_aa: {annotation.snp_aa}')
    print(f'mutation_type: {annotation.mutation_type}')
    print(f'Gene direction: {annotation.direction}')
    print(f'{annotation.ref_aa}{annotation.aa_residue_pos}{annotation.snp_aa}')

# Created 2021 by Tod Stuber