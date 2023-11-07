#!/usr/bin/env python

__version__ = "3.17"

import os
import gzip
import re
import glob
import shutil
import regex
import argparse
import textwrap
from collections import OrderedDict
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from concurrent import futures
from dask import delayed
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats


class Spoligo(Setup):

    def __init__(self, SAMPLE_NAME=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.print_run_time('Spoligotype')
        self.cpu_count_half = int(self.cpus / 2)
        real_path = os.path.dirname(os.path.realpath(__file__))
        self.spoligo_db = real_path + "/../dependencies/spoligotype_db.txt" 
        spoligo_dictionary = {}
        spoligo_dictionary["spacer01"] = ["TGATCCAGAGCCGGCGACCCTCTAT", "ATAGAGGGTCGCCGGCTCTGGATCA"]
        spoligo_dictionary["spacer02"] = ["CAAAAGCTGTCGCCCAAGCATGAGG", "CCTCATGCTTGGGCGACAGCTTTTG"]
        spoligo_dictionary["spacer03"] = ["CCGTGCTTCCAGTGATCGCCTTCTA", "TAGAAGGCGATCACTGGAAGCACGG"]
        spoligo_dictionary["spacer04"] = ["ACGTCATACGCCGACCAATCATCAG", "CTGATGATTGGTCGGCGTATGACGT"]
        spoligo_dictionary["spacer05"] = ["TTTTCTGACCACTTGTGCGGGATTA", "TAATCCCGCACAAGTGGTCAGAAAA"]
        spoligo_dictionary["spacer06"] = ["CGTCGTCATTTCCGGCTTCAATTTC", "GAAATTGAAGCCGGAAATGACGACG"]
        spoligo_dictionary["spacer07"] = ["GAGGAGAGCGAGTACTCGGGGCTGC", "GCAGCCCCGAGTACTCGCTCTCCTC"]
        spoligo_dictionary["spacer08"] = ["CGTGAAACCGCCCCCAGCCTCGCCG", "CGGCGAGGCTGGGGGCGGTTTCACG"]
        spoligo_dictionary["spacer09"] = ["ACTCGGAATCCCATGTGCTGACAGC", "GCTGTCAGCACATGGGATTCCGAGT"]
        spoligo_dictionary["spacer10"] = ["TCGACACCCGCTCTAGTTGACTTCC", "GGAAGTCAACTAGAGCGGGTGTCGA"]
        spoligo_dictionary["spacer11"] = ["GTGAGCAACGGCGGCGGCAACCTGG", "CCAGGTTGCCGCCGCCGTTGCTCAC"]
        spoligo_dictionary["spacer12"] = ["ATATCTGCTGCCCGCCCGGGGAGAT", "ATCTCCCCGGGCGGGCAGCAGATAT"]
        spoligo_dictionary["spacer13"] = ["GACCATCATTGCCATTCCCTCTCCC", "GGGAGAGGGAATGGCAATGATGGTC"]
        spoligo_dictionary["spacer14"] = ["GGTGTGATGCGGATGGTCGGCTCGG", "CCGAGCCGACCATCCGCATCACACC"]
        spoligo_dictionary["spacer15"] = ["CTTGAATAACGCGCAGTGAATTTCG", "CGAAATTCACTGCGCGTTATTCAAG"]
        spoligo_dictionary["spacer16"] = ["CGAGTTCCCGTCAGCGTCGTAAATC", "GATTTACGACGCTGACGGGAACTCG"]
        spoligo_dictionary["spacer17"] = ["GCGCCGGCCCGCGCGGATGACTCCG", "CGGAGTCATCCGCGCGGGCCGGCGC"]
        spoligo_dictionary["spacer18"] = ["CATGGACCCGGGCGAGCTGCAGATG", "CATCTGCAGCTCGCCCGGGTCCATG"]
        spoligo_dictionary["spacer19"] = ["TAACTGGCTTGGCGCTGATCCTGGT", "ACCAGGATCAGCGCCAAGCCAGTTA"]
        spoligo_dictionary["spacer20"] = ["TTGACCTCGCCAGGAGAGAAGATCA", "TGATCTTCTCTCCTGGCGAGGTCAA"]
        spoligo_dictionary["spacer21"] = ["TCGATGTCGATGTCCCAATCGTCGA", "TCGACGATTGGGACATCGACATCGA"]
        spoligo_dictionary["spacer22"] = ["ACCGCAGACGGCACGATTGAGACAA", "TTGTCTCAATCGTGCCGTCTGCGGT"]
        spoligo_dictionary["spacer23"] = ["AGCATCGCTGATGCGGTCCAGCTCG", "CGAGCTGGACCGCATCAGCGATGCT"]
        spoligo_dictionary["spacer24"] = ["CCGCCTGCTGGGTGAGACGTGCTCG", "CGAGCACGTCTCACCCAGCAGGCGG"]
        spoligo_dictionary["spacer25"] = ["GATCAGCGACCACCGCACCCTGTCA", "TGACAGGGTGCGGTGGTCGCTGATC"]
        spoligo_dictionary["spacer26"] = ["CTTCAGCACCACCATCATCCGGCGC", "GCGCCGGATGATGGTGGTGCTGAAG"]
        spoligo_dictionary["spacer27"] = ["GGATTCGTGATCTCTTCCCGCGGAT", "ATCCGCGGGAAGAGATCACGAATCC"]
        spoligo_dictionary["spacer28"] = ["TGCCCCGGCGTTTAGCGATCACAAC", "GTTGTGATCGCTAAACGCCGGGGCA"]
        spoligo_dictionary["spacer29"] = ["AAATACAGGCTCCACGACACGACCA", "TGGTCGTGTCGTGGAGCCTGTATTT"]
        spoligo_dictionary["spacer30"] = ["GGTTGCCCCGCGCCCTTTTCCAGCC", "GGCTGGAAAAGGGCGCGGGGCAACC"]
        spoligo_dictionary["spacer31"] = ["TCAGACAGGTTCGCGTCGATCAAGT", "ACTTGATCGACGCGAACCTGTCTGA"]
        spoligo_dictionary["spacer32"] = ["GACCAAATAGGTATCGGCGTGTTCA", "TGAACACGCCGATACCTATTTGGTC"]
        spoligo_dictionary["spacer33"] = ["GACATGACGGCGGTGCCGCACTTGA", "TCAAGTGCGGCACCGCCGTCATGTC"]
        spoligo_dictionary["spacer34"] = ["AAGTCACCTCGCCCACACCGTCGAA", "TTCGACGGTGTGGGCGAGGTGACTT"]
        spoligo_dictionary["spacer35"] = ["TCCGTACGCTCGAAACGCTTCCAAC", "GTTGGAAGCGTTTCGAGCGTACGGA"]
        spoligo_dictionary["spacer36"] = ["CGAAATCCAGCACCACATCCGCAGC", "GCTGCGGATGTGGTGCTGGATTTCG"]
        spoligo_dictionary["spacer37"] = ["CGCGAACTCGTCCACAGTCCCCCTT", "AAGGGGGACTGTGGACGAGTTCGCG"]
        spoligo_dictionary["spacer38"] = ["CGTGGATGGCGGATGCGTTGTGCGC", "GCGCACAACGCATCCGCCATCCACG"]
        spoligo_dictionary["spacer39"] = ["GACGATGGCCAGTAAATCGGCGTGG", "CCACGCCGATTTACTGGCCATCGTC"]
        spoligo_dictionary["spacer40"] = ["CGCCATCTGTGCCTCATACAGGTCC", "GGACCTGTATGAGGCACAGATGGCG"]
        spoligo_dictionary["spacer41"] = ["GGAGCTTTCCGGCTTCTATCAGGTA", "TACCTGATAGAAGCCGGAAAGCTCC"]
        spoligo_dictionary["spacer42"] = ["ATGGTGGGACATGGACGAGCGCGAC", "GTCGCGCTCGTCCATGTCCCACCAT"]
        spoligo_dictionary["spacer43"] = ["CGCAGAATCGCACCGGGTGCGGGAG", "CTCCCGCACCCGGTGCGATTCTGCG"]
        self.spoligo_dictionary = spoligo_dictionary

    def finding_sp(self, spacer_sequence):
        # spacer_id, spacer_sequence = spacer_id_and_spacer_sequence
        total_count = 0
        total_finds = 0
        #if total < 6: # doesn't make a big different.  Might as well get full counts
        #total += sum(seq.count(x) for x in (v)) #v=list of for and rev spacer
        total_finds = [len(regex.findall("(" + spacer + "){s<=1}", self.seq_string)) for spacer in spacer_sequence]
        for number in total_finds:
            total_count += number
        return (total_count)

    def binary_to_octal(self, binary):
        #binary_len = len(binary)
        i = 0
        ie = 1
        octal = ""
        while ie < 43:
            ie = i + 3
            # print(binary[i:ie])
            region = binary[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2: # for the lone spacer 43.  When present needs to be 1 not 4.
                    oct = 1
                else:
                    oct = 4
            else:
                oct = 0
            try:
                if int(region[1]) == 1:
                    oct += 2
                if int(region[2]) == 1:
                    oct += 1
            except IndexError:
                pass
            octal = octal + str(oct)
        return(octal)

    def spoligo(self):

        octal = None
        sbcode = None
        db_binarycode = None
        sample_binary = None

        seq_string = ""
        count_summary = {}
        sequence_list = []
        try:
            for fastq in self.FASTQ_list:
                sum_length = 0
                with gzip.open(fastq, "rt") as in_handle:
                    # all 3, title and seq and qual, were needed
                    count=0
                    for title, seq, qual in FastqGeneralIterator(in_handle):
                        if count < 1000000: #million read max
                            count+=1
                            sequence_list.append(seq)
                            sum_length = sum_length + len(seq)
                    ave_length = sum_length/count
                    ave_length = int(ave_length)
                    print(f'Spoligo calculated average read length: {ave_length}')
        except TypeError:
            # TypeError if not paired
            pass

        if ave_length > 64:
            print("Spoligo check: Average read length >=65, Reads parsed on repeat regions before counting spacers")
            #Three 10bp sequences dispersed across repeat region, forward and reverse
            capture_spacer_sequence = re.compile(".*TTTCCGTCCC.*|.*GGGACGGAAA.*|.*TCTCGGGGTT.*|.*AACCCCGAGA.*|.*TGGGTCTGAC.*|.*GTCAGACCCA.*")
            sequence_list = list(filter(capture_spacer_sequence.match, sequence_list))
            seq_string = "".join(sequence_list)
        else:
            #if <= 70 then search all reads, not just those with repeat regions.
            print("Spoligo check: Average read length < 65.  Looking at all reads for spacers.  Will be very slow. Queue Jepordy theme song.")
            seq_string = "".join(sequence_list)
        self.seq_string = seq_string
        for spacer_id, spacer_sequence in self.spoligo_dictionary.items():
            count = delayed(self.finding_sp)(spacer_sequence)
            count_summary.update({spacer_id: count})
        pull = delayed(count_summary)
        count_summary = pull.compute()
        count_summary = OrderedDict(sorted(count_summary.items()))
        spoligo_binary_dictionary = {}
        self.call_cut_off = 4
        for k, v in count_summary.items():
            if v > self.call_cut_off:
                spoligo_binary_dictionary.update({k: 1})
            else:
                spoligo_binary_dictionary.update({k: 0})
        spoligo_binary_dictionary = OrderedDict(sorted(spoligo_binary_dictionary.items()))
        spoligo_binary_list = []
        for v in spoligo_binary_dictionary.values():
            spoligo_binary_list.append(v)
        sample_binary = ''.join(str(e) for e in spoligo_binary_list)  #sample_binary correct
        self.sample_binary = sample_binary
        self.octal = self.binary_to_octal(sample_binary)
        found = False
        with open(self.spoligo_db) as spoligo_db_file: # put into dictionary or list
            for line in spoligo_db_file:
                line = line.rstrip()
                sbcode = line.split()[1]
                db_binarycode = line.split()[2]
                if sample_binary == db_binarycode:
                    found = True
                    self.sbcode = sbcode
        if not found:
            if sample_binary == '0000000000000000000000000000000000000000000':
                self.sbcode = "spoligo not found, binary all zeros, see spoligo file"
            else:
                self.sbcode = "Not Found"
        self.sample_binary = sample_binary
        self.count_summary_list=[]
        for spacer, count in count_summary.items():
            self.count_summary_list.append(count)

    def latex(self, tex):
        blast_banner = Banner("Spoligotype")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{ l | l | l }', file=tex)
        print(r'\multicolumn{3}{l}{Spacer Counts} \\', file=tex)
        print(r'\hline', file=tex) 
        count_summary = ":".join(map(str, self.count_summary_list))
        print(r'\multicolumn{3}{l}{' + f'{count_summary}' + r' } \\', file=tex)
        print(r'\hline', file=tex)
        print(f'Binary Code, threshold greater than {str(self.call_cut_off)} spacer counts & Octal Code & SB Number \\\\', file=tex)
        print(r'\hline', file=tex)
        print(f'{self.sample_binary} & {self.octal} & {self.sbcode} \\\\', file=tex)
        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        excel_dict['Spoligotype Spacer Counts'] = f'{":".join(map(str, self.count_summary_list))}'
        excel_dict['Spoligotype Binary Code'] = f'binary-{self.sample_binary}'
        excel_dict['Spoligotype Octal Code'] = f'octal-{self.octal}'
        excel_dict['Spoligotype SB Number'] = f'{self.sbcode}'


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    Mycobacterium bovis spoligotype from WGS
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='turn off map.pooling of samples')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    spoligo = Spoligo(SAMPLE_NAME=args.SAMPLE_NAME, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=args.debug)
    spoligo.spoligo()

    #Latex report
    latex_report = Latex_Report(spoligo.sample_name)
    spoligo.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(spoligo.sample_name)
    spoligo.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Updated 2021 by Tod Stuber