#!/usr/bin/env python3

__version__ = "3.32"

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
        
        # Ensure gbk_list is not None and handle file processing safely
        if gbk_list is None:
            gbk_list = []
            
        for each in gbk_list:
            try: #IF gbk provided as path variable copy local
                gbk_dict_list.append(SeqIO.to_dict(SeqIO.parse(each, "genbank")))
            except (FileNotFoundError, PermissionError, OSError) as e:
                print(f"Warning: Could not process file {each}: {e}")
                continue
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
        
        # Validate input parameters
        if not abs_pos or ':' not in abs_pos:
            print("Error: Invalid abs_pos format. Expected format: 'chrom:position'")
            return
            
        try:
            chrom, position = abs_pos.split(':', 1)  # Split only on first ':'
            position = int(position)
        except ValueError:
            print(f"Error: Could not parse position from abs_pos: {abs_pos}")
            return
            
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
        feature = None  # Initialize feature variable
        
        try:
            if chrom not in self.gbk_dict:
                raise KeyError(f"Chromosome {chrom} not found in provided GenBank files")
                
            for feature in self.gbk_dict[chrom].features:
                if "CDS" == feature.type or "tRNA" == feature.type or "rRNA" == feature.type or "repeat_region" == feature.type or "mobile_element" == feature.type or "ncRNA" == feature.type:
                    # Safely handle feature location parts
                    if hasattr(feature.location, 'parts') and feature.location.parts:
                        location_parts = feature.location.parts
                    else:
                        location_parts = [feature.location]
                        
                    for part in location_parts:
                        if position in range(part.start, part.end):
                            self.feature_found = True
                            self.direction = "forward gene"
                            self.cds_nt_start = part.start
                            self.cds_nt_end = part.end
                            aa_location = (position - feature.location.start) / 3
                            aa_residue, nt_in_aa = str(aa_location).split('.')
                            
                            # Safely extract gene information
                            try:
                                if hasattr(feature, 'qualifiers') and 'gene' in feature.qualifiers:
                                    self.gene = feature.qualifiers["gene"][0]
                            except (KeyError, IndexError, TypeError):
                                pass
                                
                            # Safely extract product information with fallback hierarchy
                            try:
                                if hasattr(feature, 'qualifiers') and 'product' in feature.qualifiers:
                                    self.product = feature.qualifiers['product'][0]
                            except (KeyError, IndexError, TypeError):
                                try:
                                    if hasattr(feature, 'qualifiers') and 'locus_tag' in feature.qualifiers:
                                        self.product = feature.qualifiers['locus_tag'][0]
                                except (KeyError, IndexError, TypeError):
                                    try:
                                        if hasattr(feature, 'qualifiers') and 'label' in feature.qualifiers:
                                            feature_type = getattr(feature, 'type', 'unknown')
                                            self.product = f'{feature_type}, {feature.qualifiers["label"][0]}'
                                    except (KeyError, IndexError, TypeError):
                                        feature_type = getattr(feature, 'type', 'unknown')
                                        self.product = f'{feature_type}, product_unknown'
                                        
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
                                print(f"Error: Unexpected nt_in_aa value: {nt_in_aa}")
                                return
                                
                            self.aa_residue_pos = int((left - part.start) / 3)
                            self.aa_residue_pos = self.aa_residue_pos + 1 # aa is zero based so add 1
                            
                            if hasattr(part, 'strand') and part.strand == -1: # Reverse complement
                                self.aa_residue_pos = int((part.end - left) / 3)
                                
                            # Safely extract sequence with bounds checking
                            try:
                                seq_length = len(self.gbk_dict[chrom].seq)
                                if left < 0 or right > seq_length:
                                    print(f"Warning: Position out of sequence bounds: left={left}, right={right}, seq_length={seq_length}")
                                    return
                                rbc = self.gbk_dict[chrom].seq[left:right]
                            except (IndexError, AttributeError) as e:
                                print(f"Error extracting sequence: {e}")
                                return
                                
                            if hasattr(part, 'strand') and part.strand == -1: # Reverse complement
                                try:
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
                                except Exception as e:
                                    print(f"Error processing reverse complement: {e}")
                                    return
                                    
                            rbc_list = list(str(rbc))  # Convert Bio.Seq to string first
                            self.reference_base_code = "".join(rbc_list)
                            
                            # Safely lookup reference amino acid
                            try:
                                self.ref_aa = self.aa_code.get(self.reference_base_code, 'unfound_ref_AA')
                            except (KeyError, TypeError):
                                self.ref_aa = 'unfound_ref_AA'
                                
                            #change rbc_list to represent SNP
                            if snp_nt is not None:
                                try:
                                    if 0 <= nt_index_aa < len(rbc_list):
                                        rbc_list[nt_index_aa] = snp_nt
                                        # Example snp_dictionary: SNP at abs pos, {'NC_017250.1:264518': 'T', ...}
                                        self.snp_base_code = "".join(rbc_list)
                                        self.snp_aa = self.aa_code.get(self.snp_base_code, 'ambiguous')
                                    else:
                                        print(f"Warning: nt_index_aa {nt_index_aa} out of range for sequence length {len(rbc_list)}")
                                        self.snp_aa = "SNP nt not provided"
                                except (TypeError, IndexError) as e:
                                    print(f"Error processing SNP: {e}")
                                    self.snp_aa = "SNP nt not provided"
                            else:
                                self.snp_base_code = "SNP nt not provided"
                                self.snp_aa = "SNP nt not provided"
                                
                            # Determine mutation type
                            if self.ref_aa == self.snp_aa:
                                self.mutation_type = "silent mutation"
                            elif self.snp_aa == "ambiguous":
                                self.mutation_type = "unsure-ambiguous"
                            elif self.ref_aa == "unfound_ref_AA":
                                self.mutation_type = ""
                            elif self.snp_aa == "SNP nt not provided":
                                self.mutation_type = "n/a"
                            else:
                                self.mutation_type = "nonsynonymous"
                            return
                            
        except KeyError as e:
            # Safely handle chromosome error with proper variable access
            chrom_safe = getattr(self, 'chrom', 'unknown')
            position_safe = getattr(self, 'position', 'unknown')
            
            print(f'\n### KeyError: incorrect chrom in dataset\n### grep -l {chrom_safe} *vcf to find file with error\n### File must be removed\n')
            
            # Safely access feature type
            if feature is not None:
                feature_type = getattr(feature, 'type', 'unknown')
                print(f'Feature Type: {feature_type}')
            else:
                print('Feature Type: not processed')
                
            print(f'Position: {position_safe}\n')
            sys.exit(0)
        except Exception as e:
            print(f"Unexpected error in annotation processing: {e}")
            return
            
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

    try:
        annotation = Annotation(gbk_list=args.gbk)
        annotation.run(args.abs_pos, args.snp_nt)
        
        # Safely print results with getattr fallbacks to maintain exact output format
        abs_pos_safe = getattr(args, 'abs_pos', 'unknown')
        print(abs_pos_safe)
        
        print(f'cds_nt_start: {getattr(annotation, "cds_nt_start", "n/a")}')
        print(f'cds_nt_end: {getattr(annotation, "cds_nt_end", "n/a")}')
        print(f'gene: {getattr(annotation, "gene", "unlisted gene")}')
        print(f'product: {getattr(annotation, "product", "No annotated product")}')
        print(f'aa_residue_pos: {getattr(annotation, "aa_residue_pos", "n/a")}')
        print(f'snp_nt: {getattr(annotation, "snp_nt", "n/a")}')
        print(f'aa_pos: {getattr(annotation, "aa_pos", "n/a")}')
        print(f'reference base code: {getattr(annotation, "reference_base_code", "n/a")}')
        print(f'snp_base_code: {getattr(annotation, "snp_base_code", "SNP nt not provided")}')
        print(f'ref_aa: {getattr(annotation, "ref_aa", "n/a")}')
        print(f'snp_aa: {getattr(annotation, "snp_aa", "n/a")}')
        print(f'mutation_type: {getattr(annotation, "mutation_type", "n/a")}')
        print(f'Gene direction: {getattr(annotation, "direction", "n/a")}')
        
        # Final mutation notation with safe attribute access
        ref_aa = getattr(annotation, "ref_aa", "n/a")
        aa_residue_pos = getattr(annotation, "aa_residue_pos", "n/a") 
        snp_aa = getattr(annotation, "snp_aa", "n/a")
        print(f'{ref_aa}{aa_residue_pos}{snp_aa}')
        
    except Exception as e:
        print(f"Error: Failed to create or run annotation: {e}")
        sys.exit(1)

# Created 2021 by Tod Stuber