#!/usr/bin/env python3

__version__ = "3.28"

import os
import sys
import shutil
import glob
import subprocess
import argparse
import textwrap
import numpy as np
from Bio import SeqIO

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats

from vsnp3_fastq_stats_seqkit import FASTQ_Stats


class Assemble(Setup):
    ''' 
    '''

    def __init__(self, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        '''
        Use file_setup to get the routine done
        '''
        Setup.__init__(self, FASTA=FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        if FASTQ_R1:
            fastq_stats = FASTQ_Stats(FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
            fastq_stats.run()
            self.R1 = fastq_stats.R1
            self.R2 = fastq_stats.R2

    def run(self,):
        '''
        Run SPAdes assembly, single read assumes ion torrent
        '''
        sample_name = self.sample_name
        FASTQ_list = self.FASTQ_list
        cwd = self.cwd
        debug = self.debug

        self.print_run_time('SPAdes')
        try:
            if len(FASTQ_list) == 2:
                #--pe-[1-2]-s assemble unequal read files
                subprocess.run(["spades.py", "--pe1-s", FASTQ_list[0], "--pe2-s", FASTQ_list[1], "-o", "spades_assembly", "--careful"], capture_output=True, check=False) # default -t 16, capture_output=True
            elif len(FASTQ_list) == 1: # assume iontorrent
                subprocess.run(["spades.py", "--iontorrent", "-s", FASTQ_list[0], "-o", "spades_assembly", "--careful"], capture_output=True, check=False)
            else:
                print(f'\n### Must have either single or paired read set.\n')
                sys.exit(0)
        except Exception as e:
            print(f"\n### Error running SPAdes: {str(e)}\n")
            sys.exit(1)

        if os.path.exists(f'{cwd}/spades_assembly/scaffolds.fasta'):
            shutil.copy2(f'{cwd}/spades_assembly/scaffolds.fasta', f'{cwd}/{sample_name}.fasta')
            self.FASTA = f'{cwd}/{sample_name}.fasta'
        else:
            print(f'\n### SPAdes did not complete, see log\n')
            self.FASTA = None
        
        if not debug:
            shutil.rmtree(f'{cwd}/spades_assembly')
        
        try:
            self.spades_version = subprocess.run(['spades.py', '-v'], capture_output=True, text=True, check=False).stdout.splitlines()[0]
        except Exception:
            self.spades_version = "SPAdes version: unable to determine"

    def stats(self, FASTA=None):
        '''
        description
        '''
        if not FASTA or not os.path.exists(FASTA):
            print(f"\n### Error: FASTA file not found or not provided: {FASTA}\n")
            sys.exit(1)
            
        records = list(SeqIO.parse(FASTA, "fasta"))
        if not records:
            print(f"\n### Error: No records found in FASTA file: {FASTA}\n")
            self.contig_count = 0
            self.small_contigs_count = 0
            self.mid_size = 0
            self.greater_one_kb_count = 0
            self.longest_contig = 0
            self.total_contig_lengths = 0
            self.n50 = 0
            self.l50 = 0
            self.mean_coverage = 0.0
            self.coverage_title = 'Coverage'
            return
            
        cov_length_list=[]
        contig_count = 0
        coverage_list=[]
        length_list=[]
        small_contigs=[]
        greater_one_kb=[]
        mid_size = []
        for rec in records:
            header = rec.description
            try:
                coverage_value = header.split('_')[5]
                coverage_value = float(coverage_value)  # Convert to float first
                coverage_value = int(coverage_value)    # Then to int to handle values with decimals
            except (IndexError, ValueError):
                coverage_value = 1
            cov_length_list.append({'name': rec.description, 'cov': coverage_value, 'length': len(rec)})
            coverage_list.append(coverage_value)
            length_list.append(len(rec))
            contig_count += 1
            if len(rec) <= 300:
                small_contigs.append(len(rec))
            elif len(rec) >= 1000:
                greater_one_kb.append(len(rec))
            else:
                mid_size.append(len(rec))
        total_contig_lengths = int(sum(length_list))

        if total_contig_lengths == 0:
            print(f"\n### Error: Total contig length is zero in FASTA file: {FASTA}\n")
            self.contig_count = contig_count
            self.small_contigs_count = len(small_contigs)
            self.mid_size = len(mid_size)
            self.greater_one_kb_count = len(greater_one_kb)
            self.longest_contig = max(length_list) if length_list else 0
            self.total_contig_lengths = total_contig_lengths
            self.n50 = 0
            self.l50 = 0
            self.mean_coverage = 0.0
            self.coverage_title = 'Coverage'
            return

        if hasattr(self, 'FASTQ_R1') and self.FASTQ_R1 and hasattr(self, 'R1'):
            self.coverage_title = 'FASTQ calculated mean coverage' #read count * read size / total assembly length
            try:
                # Handle potential commas in numeric string
                num_seqs = int(self.R1.num_seqs.replace(',', ''))
                mean_coverage = ((num_seqs * float(self.R1.avg_len))*2)/total_contig_lengths
            except (ValueError, AttributeError, ZeroDivisionError):
                # Fallback if there's an issue with the calculation
                mean_coverage = 0.0
        else:
            # No FASTQs so calculate mean coverage via SPAdes reportings.
            normalized_list = []
            for rec in records:
                header = rec.description
                try:
                    coverage_value = header.split('_')[5]
                    coverage_value = float(coverage_value)  # Convert to float first
                    coverage_value = int(coverage_value)    # Then to int
                except (IndexError, ValueError):
                    coverage_value = 1
                normalized_list.append((len(rec) / total_contig_lengths) * coverage_value)
            self.coverage_title = 'SPAdes calculated mean coverage' 
            mean_coverage = sum(normalized_list)

        #N50 calculation
        all_len = sorted(length_list, reverse=True)
        csum = np.cumsum(all_len)
        n2 = int(sum(length_list)/2)
        
        # Handle case where all contigs are smaller than N50 threshold
        if all(c < n2 for c in csum):
            self.n50 = all_len[0] if all_len else 0  # Largest contig is N50
            self.l50 = len(all_len)  # All contigs needed
        else:
            csumn2 = min(csum[csum >= n2])
            ind = np.where(csum == csumn2)
            self.n50 = all_len[int(ind[0])] # n50 smallest size contig which, along with the larger contigs, contain half of sequence of a particular genome
            self.l50 = int(ind[0][0]) + 1 # l50 smallest number of contigs whose length sum makes up half of genome
            
        self.greater_one_kb_count = len(greater_one_kb)
        self.longest_contig = int(max(length_list)) if length_list else 0
        self.small_contigs_count = len(small_contigs)
        self.mid_size = len(mid_size)
        self.contig_count = contig_count
        self.total_contig_lengths = total_contig_lengths
        self.mean_coverage = mean_coverage

        print(f'\t     Contig count: {bcolors.YELLOW}{self.contig_count:,}{bcolors.ENDC}, \n \
            Contig length counts <|301-999bp|>: {bcolors.RED}{self.small_contigs_count:,}{bcolors.ENDC}|{bcolors.BLUE}{self.mid_size:,}{bcolors.ENDC}|{bcolors.GREEN}{self.greater_one_kb_count:,}{bcolors.ENDC}, \n \
            Longest contig: {bcolors.GREEN}{self.longest_contig:,}{bcolors.ENDC}, \n \
            Total length: {bcolors.WHITE}{self.total_contig_lengths:,}{bcolors.ENDC}, \n \
            N50: {bcolors.PURPLE}{self.n50:,}{bcolors.ENDC}, \n \
            {self.coverage_title}: {bcolors.YELLOW}{self.mean_coverage:,.1f}X{bcolors.ENDC}\n')
    
    def latex(self, tex, groups=None):
        blast_banner = Banner("Assembly")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{ l | l | l | l | l | l }', file=tex)
        print(f'Contig count & Contig length counts $<$ | 301-999bp | $>$ & Longest contig & Total length & N50 & {self.coverage_title} \\\\', file=tex)
        print(r'\hline', file=tex)
        print(f'{self.contig_count:,} & {self.small_contigs_count:,} | {self.mid_size:,} | {self.greater_one_kb_count:,} & {self.longest_contig:,} & {self.total_contig_lengths:,} & {self.n50:,} & {self.mean_coverage:,.1f}X \\\\', file=tex)
        print(r'\hline', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\vspace{0.1 mm}', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\\', file=tex)
        # print(r'\begin{flushleft}Results provided by: \href{https://blast.ncbi.nlm.nih.gov/Blast.cgi}{BLAST}\end{flushleft}', file=tex)
        print(r'\end{table}', file=tex)

    def excel(self, excel_dict):
        excel_dict['Contig count'] = f'{self.contig_count:,}'
        excel_dict['Contig length counts <|301-999bp|>'] = f'{self.small_contigs_count:,}|{self.mid_size:,}|{self.greater_one_kb_count:,}'
        excel_dict['Longest contig'] = f'{self.longest_contig:,}'
        excel_dict['Total length'] = f'{self.total_contig_lengths:,}'
        excel_dict['N50'] = f'{self.n50:,}'
        excel_dict[f'{self.coverage_title}'] = f'{self.mean_coverage:,.1f}X'

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Wrapper script for SPAdes assembly.

    Usage:

    Single read:
    vsnp3_assmbly.py -r1 *.fastq.gz

    Paired reads:
    vsnp3_assmbly.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz

    If just stats is needed:
    vsnp3_assembly.py -f *fasta

    Output:  pdf report and excel stat summary

    '''), epilog='''---------------------------------------------------------''')
    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=False, default=None, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-f', '--fasta', action='store', dest='FASTA', default=None, help='provide assembly if just stats are needed.  Assumed SPAdes assembly')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    assemble = Assemble(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2)
    if args.FASTQ_R1:
        assemble.run()
        if assemble.FASTA and os.path.exists(assemble.FASTA):
            assemble.stats(assemble.FASTA)
        else:
            print('### Error: FASTA file not created or not found after assembly')
            sys.exit(1)
    elif args.FASTA:
        if os.path.exists(args.FASTA):
            assemble.stats(args.FASTA)
        else:
            print(f'### Error: FASTA file not found: {args.FASTA}')
            sys.exit(1)
    else:
        print('### Error: Provide FASTQ or FASTA file.  See vsnp3_assembly.py -h for option')
        sys.exit(1)

    #Latex report
    latex_report = Latex_Report(assemble.sample_name)
    assemble.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(assemble.sample_name)
    assemble.excel(excel_stats.excel_dict)
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

# Created 2021 by Tod Stuber