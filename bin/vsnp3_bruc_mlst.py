#!/usr/bin/env python

__version__ = "3.18"

import os
import io
import shutil
import re
import glob
import pandas as pd
from collections import OrderedDict
import argparse
import textwrap


class Bruc_MLST:

    def __init__(self, read1, read2):
        self.fastq_list = []
        for read in [read1, read2]:
            if read: #append if not None
                self.fastq_list.append(read)
        self.paired = False
        if read2:
            self.paired = True
        self.sample_name = re.sub('[_.].*', '', os.path.basename(read1))

        #https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-7-34
        write_ref = open("ST1-MLST.fasta", 'w')
        print(">ST1-MLST", file=write_ref)
        print("CGTTTCCCCAAGGAAGTGGAGGTTGCAGGCGATACGATCGATGTTGGCTACGGCCCGATCAAGGTTCATGCCGTCCGCAACCCGGCCGAACTGCCGTGGAAGGAAGAAAACGTCGATATCGCCCTTGAATGCACCGGCATTTTCACCTCGCGCGACAAGGCAGCACTTCATCTTGAAGCTGGCGCCAAGCGCGTCATCGTCTCCGCTCCCGCAGACGGTGCCGATCTCACCGTCGTCTATGGTGTCAACAACGACAAGCTGACGAAGGACCATCTGGTCATCTCCAACGCTTCGTGTACCACCAACTGCCTTGCGCCGGTGGCTCAGGTTCTCAACGATACTATCGGTATCGAAAAGGGCTTCATGACCACGATCCACTCCTATACGGGCGACCAGCCGACGCTGGACACCATGCACAAGGATCTCTACCGCGCCCGCGCCGCTGCCCTTTCCATGATCCCGACCTCGACGGGTGCGGCCAAGGCCGTCGGTCTCGTTCTGCCGGAACTGAAAGGCAAGCTCGACGGCGTTGCCATTCGCGTCCCGACCCCAAATGTCTCGGTCGTTGATCTCACCTTCATCGCCAAGCGTGAAACCACCGTTGAAGAAGTCAACAATGCGATCCGCGAAGCCGCCAATGGCCGCCTCAAGGGCATTCTCGGCTATACCGATGAGAAGCTCGTCTCGCACGACTTCAACCACGATTCCCATTCCTCGGTCTTCCACACCGACCAGACCAAGGTTATGGACGGCACCATGGTGCGTATCCTGTCGTGGTACGACAATGAATGGGGCTTCTCCAGCCGCATGAGCGACACCGCCGTCGCTTTGGGCAAGCTGATCTGATAACGGCAACGCCTCTCCTTCACTGGCGAGGCGTTTTCATTTCTTGATAAGGACCGAGAGAAGAAACATGATGTTCCGCACCCTTGACGATGCCAATGTCCAATCCAAGCGCGTGCTGGTCCGTGTTGACCTCAACGTGCCGAAAATCCGCCGTTCTGCTCGCTGGTCTTAACACCCCGGGCGTCACCACCGTGATCGAGCCGGTCATGACGCGCGATCATACGGAAAAGATGCTGCAAGACTTTGGCGCAGACCTGACGGTTGAAACCGATAAGGATGGTGTGCGCCATATCCGTATTGTCGGCCAGGGCAAGCTTACCGGCCAGACCATCGACGTGCCGGGTGATCCCTCGTCAACGGCTTTTCCGCTGGTGGCCGCCCTTCTGGTCGAAGGTTCGGAGGTCACCATCCGCAATGTGCTGATGAACCCGACCCGCACCGGCCTGATCCTGACGTTGCAGGAAATGGGGGCGGATATCGAGATCATCGATCCACGCCTTGCCGGCGGCGAGGATGTCGCCGATCTGCGCGTCAAGGCCTCGAAGCTGAAAGGCGTTGTCGTTCCGCCGGAACGTGCGCCTTCGATGATCGATGAATATCCGGTTCTGGCCATTGCCGCGTCTTTTGCGGAAGGCGAAACCGTGATGGACGGTCTCGATGAACTGCGCGTCAAGGAATCGGATCGTCTGGCGGCCGTTGCGCGCGGCCTTGAAGCCAATGGTGTCGATTGTACCGAAGGCGAGATGTCGCTGACGGTTCGTGGCCGCCCCGGCGGCAAGGGGCTGGGCGGTGGCACGGTTGCAACCCACCTCGACCACCGCATCGCGATGAGTTTCCTCGTCATGGGCCTTGCATCGGAAAAGCCGGTTACGGTGGATGACAGCACCATGATCGCCACCTCTTTCCCGGAATTCATGGGCATGATGGCGGGGCTGGGGGCGAAGATTGCCGAAAGCGGTGCAGAATGAAATCGTTCGTCGTCGCCCCGTTCATTGTCGCCATTGACGGACCGGCCGCCTCGGGCAAGGGAACCCTTGCCCGGCGGATCGCGACACATTACGGGATGCCGCATCTCGATACGGGCCTGACCTATCGCGCGGTCGCCAAGAGCCGCGCTCTGTCATTCTGGCCGTGGCAGGCCCGGTGGACGGCGACGAGATCGACCTCACCAATTGCGACTGGGTCGTGCGTCCTAAAAAGATGATCGCTGATCTGGGCTTTGAAGACGTGACCGTCCTCAATGATTTCGAGGCGCAGGCCCTTGCCGTGGTTTCGCTGGAAGGCCACCATATGGAACAGATCGGCGGCAAACCGGAGGAGGCTGTTGCCACCCGCGTCGTGCTCGGCCCCGGCACGGGCCTTGGCGTGGCAGGTCTGTTTCGCACACGTCATGCATGGGTTCCGGTTCCCGGTGAAGGCGGTCATATCGATATCGGTCCACGCACCGAACGCGACTACCAGATTTTCCCGCATATCGAACGCATCGAAGGGCGTGTCACCGGCGAGCAAATTCTTAGCGGGCGGGGCCTGCGCAACCTCTATCTGGGCATCTGCGCCGCCGACAAGATCACGCCCACCCTTGAGACGCCAGTAGACATTACATCCGCCGGACTGGACGGCAGCAATCCACAAGCCGCAGAAACGCTTGACCTCTTCGCCACCTATCTGGGGCGGCTTGCGGGCGACCTTGCGCTCATTTTCATGGCGCATGGCGGCGTTTATCTTTCGGGTGGCATCCCGGTGCGCATCCTTTCCGCCCTCAAGGCCGGTTCGTTCCGCGCAACCTTCGAGGACAAGGCCCCGCACAAGGCCATCATGCGCGACATACCGGTCCGCGTTATCACATATCAACTGGCGGCCTTAACCGGGCTTTCCGCTTTCGCCCGCACCCCCTCGCGCTTTGAAGTTTCGACCGAGGGCCGCCGCTGGCGCATGCGCCGCTAGAGCATTTCCGAGCCAAAAGTGCGAAGCGGTTCCGTTTCCCAACGAGCCGACCGCGGCTGCGCTTGCCTATGGTCTCGACAAGAGCGAAGGCAAGACCATCGCTGTCTATGACCTTGGCGGCGGTACTTTCGACGTGTCGGTTCTGGAAATCGGCGACGGCGTTTTTGAAGTGAAGTCCACCAATGGCGACACGTTCCTTGGCGGTGAAGACTTCGATATTCGTCTGGTCGAATATCTGGTTGCCGAGTTCAAGAAGGAAAGTGGCATCGACCTGAAGAACGACAAGCTTGCCCTGCAGCGCCTCAAGGAAGCTGCCGAAAAGGCCAAGATCGAACTGTCGTCCTCGCAGCAGACCGAAATCAACCTGCCGTTCATCACGGCTGACCAGACTGGCCCGAAGCATCTGGCGATCAAGCTGTCGCGCGCCAAGTTTGAAAGCCTGGTCGATGATCTCGTGCAGCGCACGGTCGAGCCGTGCAAGGCGGCGCTCAAGGATGCCGGCCTCAAGGCTGGCGAAATTGACGAAGTGGTTCTGGTCGGCGGCATGACCCGCATGCCCAAGATTCAGGAAGTCGTGAAGGCCTTCTTCGGCAAGGAACCGCACAAGGGCGTGAACCCGGATGAAGTCGTGGCCATGGGCGCGGCGATCCAGGGCGGCGTTTTGCAGGGCGACGTCAAGGACGTGCTGCTGCTCGACGTGACCCCGCTTTCGCTCGGCATTGAAACGCTGGGCGGCGTGTTCACCCGCCTGATCGAACGCAACACCACTATCCCGACCAAGAAGTCGCAGACCTTCTCCACGGCTGAGGACAACCAGTCGGCCGTGACGATCCGCGTCTTCCAGGGCGAGCGTGAAATGGCAGCCGATAACAAGCTGCTTGGACAGTTCGACCTCGTTGGCATTCCGCCACGTCCCTGCCCGGAAAGCTTGCCGATTGCCAGGAGCGCGATCCGGCCAAGTCCGAAATCTTCATCGTCGAGGGCGATTCGGCAGGCGGTTCCGCCAAGAGCGGGCGCTCGCGCCAGAATCAGGCCATTCTGCCGCTGCGCGGCAAAATCCTGAACGTGGAACGCGTGCGTTTCGACCGGATGATTTCATCCGATCAGGTGGGCACCCTCATCACGGCGCTTGGCACCTCCATCGGCAAGGATGAAACGCACGGCTTCAACGCCGACAAGCTGCGTTATCACAAGATCATCATCATGACCGACGCCGACGTCGATGGCGCCCATATTCGTACGCTTCTGCTCACCTTCTTCTTCCGGCAGATGCCGGAACTGATCGAACGCGGGCATATCTATATCGCGCAGCCGCCGCTCTATAAGGTGACACGCGGCAAGTCTTCGCAATATATCAAGAACGAAGCCGCCTTTGAGGATTTCCTCATCGAAACCGGCCTTGAAGAAACGACACTGGAACTGGTGACTGGCGAAATGCGCGCCGGGCCGGATTTGCGCTCGGTGGTGGAGGATGCGCGCACGCTGCGTCAGCTTCTGCACGGCCTGCACACCCGCTATGACCGCAGCGTGGTGGAACAGGCGGCAATTGCCGGCCTGCTCAACCCCGATGCCTCAAGGGACAATGCAACGGCACAGCATTCCGCCGATACGGTTGCCAAGCGTCTCGACATGATTTCGGAAGAGACCGAGCGCGGCTGGAGCGGCCATGTGATGGAAGACGGCGGCTATCGCTTCGAGCGTATGGTGCGCGGTGTAAAGGATATCGCCATTCTCGACATGGCCCTGCTCGGCTCGGCCGATGCCCGCCAGGTCGACCGAGATCGAGATGTATTCCCGCCTGATCCATACGGTCGATCATATCGAAGGCCGCCTGCGTGACGGCATGGATGCGTTTGACGGCTTCCTCAGCCATGCATGGGCTGTGACGGTGACAGGCGCGCCGAAGCTGTGGGCAATGCGCTTTCTTGAGGAAAACGAACGCAGCCCGCGCGCATGGTATGGCGGCGCGATCGGCATGATGCATTTCAATGGCGATATGAATACAGGGCTGACGCTGCGCACCATCCGCATCAAGGATGGTGTGGCGGAAATCCGTGCAGGGGCGACGCTTCTGTTCGATTCCAACCCTGACGAGGAAGAAGCCGAGACCGAATTGAAGGCATCGGCCATGATTGCGGCTGTGCGGGACGCACAGAAGAGCAATCAGATCGCGGAAGAAAGTGTGGCGGCAAAGGTGGGTGAGGGGGTTTCGATCCTGCTGGTCGATCACGAGGATTCCTTCGTCCATACGCTTGCCAATTATTTCCGCCAGACGGGCGCCAAGGTTTCCACCGTGCGTTCACCGGTGGCAGAGGAGATATTCGACCGCGTCAATCCCGATCTGGTGGTGTTATCGCCGGGACCGGGCTCGCCGCAGGATTTCGATTGCAAGGCGACCATCGATAAGGCGCGCAAGCGCCAGCTTCCGATTTTTGGCGTCTGCCTCGGCCTTCAGGCACTGGCGGAAGCCTATGGCGGGGCGTTGCGCCAGCTTCGCGTTCCGGTGCATGGCAAGCCTTCACGCATCCGCGTATCAAAGCCGGAGCGCATTTTCTCCGGCTTGCCGGAGGAAGTGACGGTGGGGCGTTATCATTCGATCTTCGCCGATCCTGAACGCCTGCCGGATGATTTTCTCGTCACAGCCGAAACGGAAGACGGGATCATAGCCTGCGGTGGAGGTGGTGATGGTGCCGCCGGGCTCCAGCCTGCCTGCGGATGCGGGGCTTGTCGTGTTGCCCGGCACCAAATCCACGATTGCCGATCTGCTGGCGCTGCGTGAAAACGGCTGGGACCGCGAATTGGTCGCCCATGTGAAGCGGGGCGGGCATGTGCTTGGTATTTGCGGCGGGTTTCAAATGCTTGGACGGCGGATCAGTGACCCGGCGGGTATTGAAGGCAATGTGCGCGATATCGAGGGGCTGGGCCTTCTCGATATCGAGACGATGACGGAGCCGGAAAAAGTGGTTCGCAATGTTGAGGCGGTGTCGCTGCTGCATGATGAGCCGCTGGAGGGCTATGAAATCCACATCGGGCGCACCAGCGGGCCGGATATGGCGCGGCCATTTGCGCGTATCGGCGATCATGATGATGGGGCCGTCTCGCCCGATGGTCGTATCATGGGAACCTATCTCCACGGTATTTTCAGTGCGGATCGTTTCCGCCACCACTTTTTGCGCGCGCTGGGTGTGGAAGGCGGCCAGATGAATTATCGCGAGAGCGTCGAAGAGGCTCTGGGCGAACTGGCTGAAGGGCTGGAAGCCTCGCTGGATATTGATGGCCTGTTTGCGCTGGCATGATTGACGCCGCGAAGCCGAAAGCCTAGTGTCAAACCATGTGACAGGTTTTGCCGGAACGAATCCCCGGCAATACCAAAAGGGAATGCGACGGACGGACCCACGCCGGGCGTCTTTATCGCAGCCGACCCCGCGACTGTAGAGCGGAGAGGGAAGAGGCAAGCCGGGCAACCGGCAGCCACTGGAAATCAGATGCGATAATGCAACATCGCATTTTTGCCATCTTCTCGACAGATTATCTCCACACAATGGGGCATTTCGTGCCGCAATTACCCTCGATATGTCACCCCTGTCAGCGCGGCATGGGCGGTTTACTCCCGATGCTGCCCGCCCGATAAGGGACCGCGCAAAACGTAATTTGTGTAAGGAGAATGCCATGCGCACTCTTAAGTCTCTCGTAATCGTCTCGGCTGCGCTGCTGCCGTTCTCTGCGACCGCTTTTGCTGCCGACGCCATCCAGGAACAGCCTCCGGTTCCGGCTCCGGTTGAAGTAGCTCCCCAGTATAGCTGGGCTGGTGGCTATACCGGTCTTTACCTTGGCTATGGCTGGAACAAGGCCAAGACCAGCACCGTTGGCAGCATCAAGCCTGACGATTGGAAGGCTGGCGCCTTTGCTGGCTGGAACTTCCAGCAGGACCAGATCGTATACGGTGTTGAAGGTGATGCAGGTTATTCCTGGGCCAAGAAGTCCAAGGACGGCCTGGAAGTCAAGCAGGGCTTTGAAGGCTCGCTGCGTGCCCGCGTCGGCTACGACCTGAACCCGGTTATGCCGTACCTCACGGCTGGTATTGCCGGTTCGCAGATCAAGCTTAACAACGGCTTGGACGACGAAAGCAAGTTCCGCGTGGGTTGGACGGCTGGTGCCGGTCTCGAAGCCAAGCTGACGGACAACATCCTCGGCCGCGTTGAGTACCGTTACACCCAGTACGGCAACAAGAACTATGATCTGGCCGGTACGACTGTTCGCAACAAGCTGGACACGCAGGATATCCGCGTCGGCATCGGCTACAAGTTCTAATTATAGCATAATTGGACACGGAAAACCGGACAGGCAACTGTCCGGTTTTTTGTTGTCTGCAAAGGTCGAGAAAGCGCGGCAGAGCAACGGCGGCAGCCTGATTTTCAGGGGAAATGAAGTGGAGGCTTCTGTTGCCAGGTGCCTCCGAACCCCGCCTTAAGGGGCTAACCCTAAGGACTTTAGAGTGGGTTTCCCGCACCGCCATTAGGCAGCGAGAGCATAACCCTGAGCATTGTTGTCATTTGCAACTACTCTGTTGACCCGATAACGGTGGTATCATGCCGAGTAAAAGAGCGATCTTTACACCCTTGTCGATCCTGTTTCGCCCCCGCCACAACACAGCCTGATCGGCAAGCTGTGCTGTGGTGGAGGCGCCGGGTACCGCCCCCGGGTCCAATGGGTTTATTACACCGTCCGTTTATCACCATAGTCGGCTTGCGCCGACAGGACGTATATAGGCGTGGTTTTTACCGATTGGAAGGGGGCTTGTGCGTTTTCGCGCAAGACCGACAGAGGTGGTGCGGCCCTTCCGTTCATTTTCCATTGACAGCTTCCGCGTGCTGGTCAATCCTCACAATATATCGGGATCGGCCTTGAAGAGGCTTGGCGCAGCCGGGGCGGAAACCATGGCTGAAACGGGGACGATATGCCCCAATCGAAGGAGAGTGGATATATGAGTGAATATCTCGCGGATGTCCGTCGCTATGATGCTGCCGCCGATGAGGCCGTTGTCGAGAAAATCGTCAAGCATCTTGGCATTGCGCTTCGCAATCGCGATTCCTCGCTCGTTTCGGCAAGC", file=write_ref)
        write_ref.close()

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

    def find_mlst(self):

        sample_name = self.sample_name + "_mlst"
        fastq_list = self.fastq_list
        mlst_reference = "ST1-MLST.fasta"
        ref = re.sub('.fasta', '', os.path.basename(mlst_reference))
        os.system("samtools faidx {}  2> /dev/null" .format(mlst_reference))
        os.system("picard CreateSequenceDictionary REFERENCE={} OUTPUT={}  2> /dev/null" .format(mlst_reference, ref + ".dict"))
        os.system("bwa index {}  2> /dev/null" .format(mlst_reference))
        samfile_mlst = sample_name + ".sam"
        print(f'### {sample_name} BWA Aligning reads for Brucella MLST...')
        if self.paired:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {mlst_reference} {fastq_list[0]} {fastq_list[1]} > {samfile_mlst}')
        else:
            os.system(f'bwa mem -M -R "@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250" -t 8 {mlst_reference} {fastq_list[0]} > {samfile_mlst}')
        print(f'{sample_name} finding MLST type...')
        all_bam_mlst = sample_name + "_all.bam"
        os.system("picard AddOrReplaceReadGroups INPUT={} OUTPUT={} RGLB=lib1 RGPU=unit1 RGSM={} RGPL=illumina  2> /dev/null" .format(samfile_mlst, all_bam_mlst, sample_name))
        mapbam = sample_name + "_mapped.bam"
        os.system("samtools view -h -F4 -b -T {} {} -o {} 2> /dev/null" .format(mlst_reference, all_bam_mlst, mapbam))
        sortedbam = sample_name + "_sorted.bam"
        os.system("samtools sort {} -o {} 2> /dev/null" .format(mapbam, sortedbam))
        os.system("samtools index {} 2> /dev/null" .format(sortedbam))
        unfiltered_vcf_mlst = sample_name + "_unfiltered_mlst" + ".vcf"
        mapq_fix = sample_name + "_mapq_fix_mlst.vcf"
        vcf_mlst = sample_name + ".vcf"
        os.system("freebayes -E -1 -f {} {} > {}" .format(mlst_reference, sortedbam, unfiltered_vcf_mlst))
        # "fix" MQ notation in VCF to match GATK output
        write_fix = open(mapq_fix, 'w+')
        with open(unfiltered_vcf_mlst, 'r') as unfiltered:
            for line in unfiltered:
                line = line.strip()
                new_line = re.sub(r';MQM=', r';MQ=', line)
                print(new_line, file=write_fix)
            write_fix.close()
        # remove clearly poor positions
        os.system(r'vcffilter -f "QUAL > 20" %s > %s' % (mapq_fix, vcf_mlst))
        pos_call_dict = {}

        df = self.read_vcf(vcf_mlst)
        pos_call_dict = dict(zip(df['POS'], df['ALT']))
        # Position 1629 was too close to the end of glk sequence.  Reads would not assemble properly to call possilbe SNP, therefore 100 bases of the gene were added.  Because of this all positions beyond this point are 100 more.  Same with position 1645 and 2693.

        target_pos_ref = {231: 'C', 297: 'T', 363: 'C', 398: 'C', 429: 'C', 523: 'G', 631: 'G', 730: 'G', 1247: 'G', 1296: 'C', 1342: 'G', 1381: 'A', 1648: 'C', 1685: 'C', 1741: 'C', 1754: 'G', 2165: 'A', 2224: 'T', 2227: 'C', 2297: 'G', 2300: 'A', 2344: 'A', 2352: 'G', 2403: 'C', 2530: 'G', 2557: 'G', 2578: 'G', 2629: 'A', 3045: 'A', 3054: 'G', 3118: 'G', 3295: 'C', 3328: 'C', 3388: 'A', 3966: 'C', 3969: 'G', 4167: 'G', 4271: 'C', 4296: 'G', 4893: 'C', 4996: 'G', 4998: 'T', 5058: 'G', 5248: 'A', 5672: 'G', 5737: 'C', 5928: 'A', 5963: 'G', 5984: 'C', 5987: 'C', 6025: 'G', 6045: 'G', 6498: 'G', 6499: 'C', 6572: 'A', 6627: 'T', 6715: 'C', 6735: 'T', 6745: 'G', 6785: 'T', 6810: 'C', 6828: 'C', 6845: 'C', 6864: 'G', 6875: 'C', 7382: 'G', 7432: 'G', 7464: 'G', 7594: 'G', 7660: 'T', 7756: 'A'}

        #pos_call_dict will replace target_pos_ref
        for key, value in pos_call_dict.items():
            if key in target_pos_ref.keys():
                target_pos_ref[key] = value
        ordered_combined_dict = OrderedDict(sorted(target_pos_ref.items()))
        combined_value_list = list(ordered_combined_dict.values())
        mlst_join = ''.join(combined_value_list)

        mlst_dictionary = {}
        mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 01"
        mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 02"
        mlst_dictionary["CTCCCGTGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 03"
        mlst_dictionary["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCAAGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 04"
        mlst_dictionary["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGAGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 05"
        mlst_dictionary["TTCCTGGGGCAACCCGAGCGAGGCAGGGAGGCCGCGGCTCGTGAGCGGTCGGGCATCTGTCCCGCGGGGTA"] = "MLST type 06"
        mlst_dictionary["CTTCCTGGCCGAGCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 07"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCCGTCTCGCGGTGCT"] = "MLST type 08"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGTGGTGCT"] = "MLST type 09"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGGGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 10"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 11"
        mlst_dictionary["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"] = "MLST type 12"
        mlst_dictionary["CCCCCGGGCCGACTCGAGCGAAGCGAAGAGGCCACGGCGCGTGAGTGACCAGGCACCTATCCCACGGGGTA"] = "MLST type 13"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAGTGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 14"
        mlst_dictionary["CCCCCGGGCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 15"
        mlst_dictionary["CCCCCGGCCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"] = "MLST type 16"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGGTA"] = "MLST type 17"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGCTA"] = "MLST type 18"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGACACGGCGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 19"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGTCCCGCAGGGTA"] = "MLST type 20"
        mlst_dictionary["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGCCCCGCAGGGTA"] = "MLST type 21"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAGGCGGGGAGGCCACGGCGCGGGAGTGGCCAGACACCTGTCCTGCGGGGTA"] = "MLST type 22"
        mlst_dictionary["CCCCCGGGCTGACCCGAGCGAAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23"
        mlst_dictionary["CCCCCGGGCTGACCCGAGCGGAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 23x"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCAAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 24"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 25"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAAGCACCTGTTCCGCGGGGTA"] = "MLST type 26"
        mlst_dictionary["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCACCTGTCCCGCGGGGTA"] = "MLST type 27"
        mlst_dictionary["CCCTCGGGCCGACCTGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCTCCTGTCCCGCGGGGTA"] = "MLST type 28"

        remove_files = glob.glob('ST1-MLST*')
        for i in remove_files:
            os.remove(i)
        remove_files = glob.glob('*-mlst*')
        for i in remove_files:
            os.remove(i)
        remove_files = glob.glob('*_mlst.vcf.idx')
        for i in remove_files:
            os.remove(i)
        os.remove(unfiltered_vcf_mlst)

        write_out = open("mlst.txt", 'w')
        if mlst_join in mlst_dictionary:
            mlst_type = mlst_dictionary[mlst_join]
            print(mlst_type)
            print(mlst_type, file=write_out)
        else:
            print("NO MLST MATCH FOUND\n")
            print("NO MLST MATCH FOUND", file=write_out)
            mlst_type="No Brucella MLST type found"
        write_out.close()
        self.mlst_type = mlst_type

        os.makedirs("mlst")
        shutil.move(vcf_mlst, "mlst")
        shutil.move("mlst.txt", "mlst")
        os.remove(samfile_mlst)
        os.remove(all_bam_mlst)
        os.remove(mapbam)
        os.remove(sortedbam)
        os.remove(f'{sortedbam}.bai')
        os.remove(mapq_fix)

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Based on https://bmcmicrobiol.biomedcentral.com/articles/10.1186/1471-2180-7-34

    Brucella MLST from WGS sequence
    
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--read1', action='store', dest='read1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--read2', action='store', dest='read2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    read1 = args.read1
    read2 = args.read2

    mlst = Bruc_MLST(read1, read2)
    mlst.find_mlst()