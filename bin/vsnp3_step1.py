#!/usr/bin/env python

__version__ = "3.25"

import os
import sys
import shutil
import glob
import re
import argparse
import textwrap

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats

from vsnp3_fastq_stats_seqkit import FASTQ_Stats
from vsnp3_best_reference_sourmash import Best_Reference
from vsnp3_reference_options import Ref_Options
from vsnp3_alignment_vcf import Alignment
from vsnp3_spoligotype import Spoligo
from vsnp3_group_reporter import GroupReporter
from vsnp3_fasta_to_fastq import Fasta_to_Paired_Fastq


class vSNP3_Step1(Setup):
    ''' 
    '''
    def __init__(self, SAMPLE_NAME=None, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, gbk=None, reference_type=None, nanopore=False, assemble_unmap=None, spoligo=None, debug=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1)
        self.assemble_unmap = assemble_unmap
        self.spoligo = spoligo
        self.latex_report = Latex_Report(self.sample_name)
        self.excel_stats = Excel_Stats(self.sample_name)
        self.nanopore = nanopore
        if FASTA: #IF -f REFERENCE FASTA PROVIDED USE IT, ie .fasta
            concatenated_FASTA = self.concat_fasta(FASTA) # -f option for FASTA will take wildcard for multiple FASTAs
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=gbk, debug=debug)
            self.reference_type = None
            with open(self.FASTA) as f:
                self.top_header_description = f.readline()
        # elif os.path.isdir(reference_type): # -t directory name with full path
        #     reference_options = Ref_Options(reference_type)
        #     concatenated_FASTA = self.concat_fasta(reference_options.fasta)
        #     self.excel_stats.excel_dict["Reference"] = f'{reference_type} Forced'
        #     Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=reference_options.gbk, debug=debug)
        #     self.reference_type = reference_type
        #     with open(self.FASTA) as f:
        #         self.top_header_description = f.readline()
        elif reference_type: # -t directory name
            reference_options = Ref_Options(reference_type)
            concatenated_FASTA = self.concat_fasta(reference_options.fasta)
            self.excel_stats.excel_dict["Reference"] = f'{reference_type} Forced'
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=reference_options.gbk, debug=debug)
            self.reference_type = reference_type
            with open(self.FASTA) as f:
                self.top_header_description = f.readline()
        else: #IF NO REFERENCE PROVIDED SEEK A "BEST REFERENCE", for tb, brucella, paraTB or SARS-CoV-2
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2,)
            self.best_reference = Best_Reference(FASTQ_R1=FASTQ_R1)
            self.best_reference.run()
            self.top_header_description = self.best_reference.top_fasta_header
            if self.best_reference.reference_set:
                reference_options = Ref_Options(self.best_reference.reference_set)
            else:
                print(f'#### {self.sample_name} PROVIDE A REFERENCE, REFERENCE NOT FOUND')
                self.excel_stats.excel_dict["Reference"] = f'PROVIDE A REFERENCE, REFERENCE NOT FOUND'
                print(f'Getting Stats before exiting')
                fastq_stats = FASTQ_Stats(FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
                fastq_stats.run()
                fastq_stats.latex(self.latex_report.tex)
                fastq_stats.excel(self.excel_stats.excel_dict)
                self.latex_report.latex_ending()
                self.excel_stats.post_excel()
                temp_dir = './temp'
                if not os.path.exists(temp_dir):
                    os.makedirs(temp_dir)
                files_grab = []
                for files in ('*.aux', '*.log', '*tex', '*png', '*out', "*_seqkit_stats.txt"):
                    files_grab.extend(glob.glob(files))
                for each in files_grab:
                    shutil.move(each, temp_dir)

                if args.debug is False:
                    shutil.rmtree(temp_dir)
                sys.exit()
            concatenated_FASTA = self.concat_fasta(reference_options.fasta)
            self.gbk = reference_options.gbk
            self.excel_stats.excel_dict["Reference"] = f'{reference_options.select_ref} by Best Reference'
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=reference_options.gbk, debug=debug)
            self.best_reference.latex(self.latex_report.tex)
            self.best_reference.excel(self.excel_stats.excel_dict)
            if reference_options.fasta == None:
                print('Sourmash unable to find a suitable match.  Provide a FASTA (-f) or directory name with dependecies (-n)')
                exit(0)
            self.reference_type = self.best_reference.reference_set

    def concat_fasta(self, FASTAs):
            basenames=[]
            fasta_name=[]
            for each in FASTAs:
                basenames.append(os.path.basename(each))
                fasta_name.append(re.sub('\..*', '', os.path.basename(each)))
            self.excel_stats.excel_dict["FASTA/s"] = ", ".join(basenames)
            #concatenate FASTA list
            concatenated_TEMP = f'{"_".join(fasta_name)}.temp'
            concatenated_FASTA = f'{"_".join(fasta_name)}.fasta'
            with open(concatenated_TEMP,'wb') as wfd:
                for each in FASTAs:
                    with open(each,'rb') as fd:
                        shutil.copyfileobj(fd, wfd)
            shutil.move(concatenated_TEMP, concatenated_FASTA)
            return concatenated_FASTA

    def run(self,):
        '''
        description
        '''
        fastq_stats = FASTQ_Stats(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
        fastq_stats.run()
        fastq_stats.latex(self.latex_report.tex)
        fastq_stats.excel(self.excel_stats.excel_dict)

        MYCO=None # Test for Mycobacterium
        try:
            if any(x in self.top_header_description for x in ['Mycobacterium bovis', 'Mycobacterium tuberculosis', 'Mycobacterium caprae', 'Mycobacterium orygis']):
                MYCO = True
        except AttributeError:
            pass
        self.MYCO = MYCO
        if self.spoligo:
            spoligo = Spoligo(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
            spoligo.spoligo()
            spoligo.latex(self.latex_report.tex)
            spoligo.excel(self.excel_stats.excel_dict)

        if int(float(fastq_stats.R1.max_len.replace(',', ''))) > 701:
            nanopore = True
        elif self.nanopore:
            nanopore = True
        else:
            nanopore = False
            
        alignment = Alignment(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, reference=self.reference, nanopore=nanopore, gbk=self.gbk, assemble_unmap=self.assemble_unmap, debug=self.debug)
        alignment.run()

        groups = "group file not provided"
        if self.reference_type:
            try:
                group_reporter = GroupReporter(alignment.zero_coverage_vcf_file_path, self.reference_type)
                groups = ", ".join(group_reporter.get_groups())
                self.excel_stats.excel_dict['Groups'] = groups
            except ValueError:
                self.excel_stats.excel_dict['Groups'] = groups
        else:
            self.excel_stats.excel_dict['Groups'] = groups
        

        alignment.latex(self.latex_report.tex, groups)
        alignment.excel(self.excel_stats.excel_dict)

        #FASTQ usability
        if alignment.zero_coverage.ave_coverage < 40 or alignment.DUPLICATION_RATIO > .80 or float(fastq_stats.R1.passQ20) < 50.0:
            fastq_usability = 'Poor'
            color_quality = bcolors.RED
        elif alignment.zero_coverage.ave_coverage < 70 or alignment.DUPLICATION_RATIO > .10 or float(fastq_stats.R1.passQ30) < 70.0:
            fastq_usability = 'Questionable'
            color_quality = bcolors.YELLOW
        else:
            fastq_usability = 'Acceptable'
            color_quality = bcolors.GREEN
        # self.excel_stats.excel_dict['FASTQ Usability'] = fastq_usability
        print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}{fastq_usability}{bcolors.ENDC} FASTQ Usability')

        #Reference usability
        if alignment.zero_coverage.genome_coverage < 0.95 or alignment.zero_coverage.percent_ref_with_zero_coverage > 1.0 or alignment.zero_coverage.percent_ref_with_good_snp_count > 0.09:
            reference_usability = 'Poor'
            color_quality = bcolors.RED
            reference_usability_test = False
        elif alignment.zero_coverage.genome_coverage < 0.975 or alignment.zero_coverage.percent_ref_with_zero_coverage > 0.5 or alignment.zero_coverage.percent_ref_with_good_snp_count > 0.03:
            reference_usability = 'Questionable'
            color_quality = bcolors.YELLOW
            reference_usability_test = False
        else:
            reference_usability = 'Acceptable'
            color_quality = bcolors.GREEN
            reference_usability_test = True
        # self.excel_stats.excel_dict['Reference Usability'] = reference_usability
        print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}{reference_usability}{bcolors.ENDC} Reference Usability')
        
        #Contamination determination
        # if reference_usability_test: #if a good reference unmapped reads is relavent for determining posible sample contamination
        #     total_reads = alignment.READ_PAIRS_EXAMINED + alignment.UNPAIRED_READS_EXAMINED + alignment.UNMAPPED_READS
        #     freq_unmapped_reads = alignment.UNMAPPED_READS / total_reads
        #     if freq_unmapped_reads > .1: #10%
        #         sample_contamination = 'High Unmapped Reads'
        #         color_quality = bcolors.RED
        #     elif freq_unmapped_reads > .05: #5%
        #         sample_contamination = 'Little Unmapped Reads'
        #         color_quality = bcolors.YELLOW
        #     else:
        #         sample_contamination = 'Insignificant Unmapped Reads' #<5%
        #         color_quality = bcolors.GREEN
        #     # self.excel_stats.excel_dict['Sample Contamination'] = sample_contamination
        #     print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}{sample_contamination}{bcolors.ENDC} Sample Contamination')
        # else:
        #     # self.excel_stats.excel_dict['Sample Contamination'] = 'N/A - Poor or Questionable Reference'
        #     print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}N/A - Poor or Questionable Reference{bcolors.ENDC} Sample Contamination')

        self.programs = alignment.programs
        self.alignment_vcf_run_summary = alignment.alignment_vcf_run_summary
        self.latex_report.latex_ending()
        self.excel_stats.post_excel()

        if self.remove_copied_gbk: #these were copied from a different location and are not need for any further analyses.
            for each in self.gbk:
                os.remove(each)


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Conda:
    conda install vsnp3=3.13 -c conda-forge -c bioconda
    ---------------------------------------------------------

    When running samples through step1 and 2 of vSNP, or when running a routine analysis, set up dependencies using vsnp3_path_adder.py

    See vsnp3_path_adder.py -h for adding a reference type and for more information

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -t NC_045512_wuhan-hu-1

    For a few reference types (not all) sourmash can be used to find the best reference.  See vsnp3_best_reference_sourmash.py -h for list of references.

    If working or suspect one of these reference types the dependencies will automatically be found

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz

    When running a one-off dependencies can be provided explicitly.

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -f *fasta -b *gbk
    or
    vsnp3_step1.py -r1 ../myfastqs/*_R1*.fastq.gz -r2 ../myfastqs/*_R2*.fastq.gz -f ../myreference/*fasta -b ../myreference/*gbk

    At completion of running script it is recommended to store sample folders by reference.  For each reference there should be a "step1" and "step2" subfolder.  These sample folders should go into step1 subfolders for the reference type used to generate SNPs stored in the VCF file.
    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, help='Provide R1 FASTQ gz file.  A single read file can also be supplied to this option')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: provide R2 FASTQ gz file')
    parser.add_argument('-f', '--FASTAtoFASTQ', action='store', dest='FASTAtoFASTQ', required=False, help='Input a FASTA file, convert to paired FASTQ files, and run.')
    parser.add_argument('-r', '--FASTA', nargs='*', dest='FASTA', required=False, help='FASTA file to be used as reference.  Multiple can be specified with wildcard')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', required=False, default=None, help="Optional: Provide directory name with FASTA and GBK file/s")
    parser.add_argument('-p', '--nanopore', action='store_true', dest='nanopore', default=False, help='if true run alignment optimized for nanopore reads')
    parser.add_argument('-o', '--output_dir', action='store', dest='output_dir', required=False, default=None, help="Optional: Provide a name.  This name will be a directory output files are writen to.  Name can be a directory path, but doesn't have to be.")
    parser.add_argument('-assemble_unmap', '--assemble_unmap', action='store_true', dest='assemble_unmap', help='Optional: skip assembly of unmapped reads.   See also vsnp3_assembly.py')
    parser.add_argument('-spoligo', '--spoligo', action='store_true', dest='spoligo', help='Optional: get spoligotype if TB complex.  See also vsnp3_spoligotype.py')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='keep spades output directory')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    if not args.FASTQ_R1 and not args.FASTAtoFASTQ:
        parser.error('No FASTQ or FASTA provided.  Use -r1 or -f option')

    if args.FASTAtoFASTQ:
        fasta_to_paired_fastq = Fasta_to_Paired_Fastq(args.FASTAtoFASTQ, coverage=100, read_length=300)
        args.FASTQ_R1 = fasta_to_paired_fastq.fastq_r1_file
        args.FASTQ_R2 = fasta_to_paired_fastq.fastq_r2_file
        print(f'Conversion to FASTQ completed')

    program_list = ['Python', 'Bio', 'numpy', 'pandas', 'scipy',]
    python_programs = []
    for name, module in sorted(sys.modules.items()): 
        if hasattr(module, '__version__') and name in program_list: 
            python_programs.append(f'{name}, {module.__version__}')

    if args.output_dir:
        output_dir = args.output_dir
        output_dir = os.path.expanduser(output_dir)
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)
        if args.FASTQ_R1:
            shutil.copy(args.FASTQ_R1, output_dir)
        if args.FASTQ_R2:
            shutil.copy(args.FASTQ_R2, output_dir)
        os.chdir(output_dir)

    vsnp = vSNP3_Step1(SAMPLE_NAME=args.SAMPLE_NAME, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, FASTA=args.FASTA, gbk=args.gbk, reference_type=args.reference_type, nanopore=args.nanopore, assemble_unmap=args.assemble_unmap, spoligo=args.spoligo, debug=args.debug)
    vsnp.run()

    with open(f'{vsnp.sample_name}_run_log.txt', 'w') as run_log:
        print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:', file=run_log)
        print(args, file=run_log)
        if args.FASTAtoFASTQ:
            print(f'Converted FASTA to FASTQ', file=run_log)
        print('\nCall Summary:', file=run_log)
        for each in vsnp.alignment_vcf_run_summary:
            print(each, file=run_log)
        print('\nVersions:', file=run_log)
        print(f'vSNP3: {__version__}', file=run_log)
        for each in python_programs:
            print(each, file=run_log)
        for each in vsnp.programs:
            print(each, file=run_log)

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*_report.out', '*.aux', '*.log', '*tex', '*png', "*_seqkit_stats.txt"):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created 2021 by Tod Stuber
