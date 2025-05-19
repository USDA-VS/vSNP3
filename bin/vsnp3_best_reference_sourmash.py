#!/usr/bin/env python3

__version__ = "3.29"

import os
import subprocess
import shutil
import glob
import pandas as pd
import argparse
import textwrap
from Bio import SeqIO

from vsnp3_file_setup import Setup
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats


class Best_Reference(Setup):
    ''' 
    '''

    def __init__(self, SAMPLE_NAME=None, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1)
        self.script_path = os.path.dirname(os.path.realpath(__file__))
        self.debug = debug

    def run(self,):
        self.print_run_time('Best Reference Finding with Sourmash')
        '''
        # Build database
        for i in *fasta; do sourmash sketch dna $i --name-from-first; done
        sourmash index ref_db.sbt.zip ./*.sig # this .zip will be placed as a dependency file.
        # Prepare read
        sourmash sketch dna *_R1*.fastq.gz
        # Search 
        sourmash search *_R1*.fastq.gz.sig ../sourmash/ref_db.sbt.zip -o sourmash_findings.csv
        '''
        all_ref_options = []
        ref_options_file = os.path.abspath(f'{self.script_path}/../dependencies/reference_options_paths.txt')
        self.ref_options_file = ref_options_file
        
        # Check if reference options file exists
        if not os.path.exists(ref_options_file):
            print(f"Error: Reference options file not found: {ref_options_file}")
            self.top_header_found = "No Reference Options File Found"
            self.top_header = "No Reference Options File Found"
            self.top_fasta = None
            self.reference_set = None
            self.top_fasta_header = "No Reference Options File Found"
            self.sourmash_df = pd.DataFrame()
            return
            
        # Read reference paths
        with open(f'{ref_options_file}', 'r') as dep_paths:
            dependency_paths = [line.strip() for line in dep_paths]
            
        # Collect all reference options from the specified paths
        for path in dependency_paths:
            if os.path.exists(path):
                ref_options = glob.glob(f'{path}/*')
                all_ref_options = all_ref_options + ref_options
            else:
                print(f"Warning: Reference path does not exist: {path}")
                
        # Filter to only include directories
        all_ref_options = [x for x in all_ref_options if os.path.isdir(x)]
        
        # Get FASTA files from each reference directory
        self.fasta_list = []
        for each_path in all_ref_options:
            self.fasta_list.extend(glob.glob(f'{each_path}/*.fasta'))
            
        # Create dictionary mapping FASTA headers to file paths
        header_dict = {}
        for fasta_path in self.fasta_list:  # each available reference/fasta
            try:
                identifiers = [seq_record.description for seq_record in SeqIO.parse(fasta_path, 'fasta')]
                if identifiers:  # Only add if there are identifiers (file not empty)
                    header_dict[identifiers[0]] = fasta_path
            except Exception as e:
                print(f"Warning: Could not parse FASTA file {fasta_path}: {str(e)}")
                
        # Check if FASTQ file exists
        if not os.path.exists(self.FASTQ_R1):
            print(f"Error: FASTQ file not found: {self.FASTQ_R1}")
            self.top_header_found = "FASTQ File Not Found"
            self.top_header = "FASTQ File Not Found"
            self.top_fasta = None
            self.reference_set = None
            self.top_fasta_header = "FASTQ File Not Found"
            self.sourmash_df = pd.DataFrame()
            return
            
        # Check if sourmash database exists
        sourmash_db = f'{self.script_path}/../dependencies/ref_db.sbt.zip'
        if not os.path.exists(sourmash_db):
            print(f"Error: Sourmash database not found: {sourmash_db}")
            self.top_header_found = "Sourmash Database Not Found"
            self.top_header = "Sourmash Database Not Found"
            self.top_fasta = None
            self.reference_set = None
            self.top_fasta_header = "Sourmash Database Not Found"
            self.sourmash_df = pd.DataFrame()
            return
            
        # Run sourmash sketch
        try:
            sketch_result = subprocess.run(
                ["sourmash", "sketch", "dna", self.FASTQ_R1],
                capture_output=True,
                text=True,
                check=False
            )
            if sketch_result.returncode != 0:
                print(f"Warning: Sourmash sketch command failed: {sketch_result.stderr}")
        except Exception as e:
            print(f"Error running sourmash sketch: {str(e)}")
            self.top_header_found = "Sourmash Sketch Failed"
            self.top_header = "Sourmash Sketch Failed"
            self.top_fasta = None
            self.reference_set = None
            self.top_fasta_header = "Sourmash Sketch Failed"
            self.sourmash_df = pd.DataFrame()
            return
            
        # Run sourmash search
        try:
            search_result = subprocess.run(
                [
                    "sourmash", 
                    "search", 
                    f'{self.FASTQ_R1}.sig', 
                    sourmash_db,
                    "-o", 
                    f'{self.sample_name}_search.csv', 
                    '--threshold=0.001'
                ],
                capture_output=True,
                text=True,
                check=False
            )
            if search_result.returncode != 0:
                print(f"Warning: Sourmash search command failed: {search_result.stderr}")
        except Exception as e:
            print(f"Error running sourmash search: {str(e)}")
            self.top_header_found = "Sourmash Search Failed"
            self.top_header = "Sourmash Search Failed"
            self.top_fasta = None
            self.reference_set = None
            self.top_fasta_header = "Sourmash Search Failed"
            self.sourmash_df = pd.DataFrame()
            if os.path.exists(f'{self.FASTQ_R1}.sig'):
                os.remove(f'{self.FASTQ_R1}.sig')
            return
            
        # Read search results
        try:
            if os.path.exists(f'{self.sample_name}_search.csv') and os.path.getsize(f'{self.sample_name}_search.csv') > 0:
                self.sourmash_df = pd.read_csv(f'{self.sample_name}_search.csv')
            else:
                print("Warning: Sourmash search produced no results or empty file")
                self.sourmash_df = pd.DataFrame()
        except Exception as e:
            print(f"Error reading sourmash search results: {str(e)}")
            self.sourmash_df = pd.DataFrame()

        #Force a top hit to a specific reference, ie TB lineages to
        try: 
            self.top_header_found = self.sourmash_df['name'][0].split()[0] # top hit
        except (IndexError, KeyError):
            self.top_header_found = "No Sourmash Findings"
            self.top_fasta_header = "No Sourmash Findings"
            
        self.top_header = self.top_header_found #default
        
        # Force non-bovis/caprae lineages to H37
        force_to_tb = ['NZ_CP041790.1', 'CP023623.1', 'CP023635.1', 'CP063804.1', 'NZ_CP022014.1', 'NZ_CP041790.1', 'NZ_CP041803.1', 'NZ_CP041869.1', 'NZ_CP041875.1']
        if self.top_header_found in force_to_tb:
             self.top_header = 'NC_000962.3'

        # Force caprae to bovis
        force_to_bovis = ['CP016401.1',]
        if self.top_header_found in force_to_bovis:
             self.top_header = 'NC_002945.4'

        self.top_fasta = None #default to no findings
        self.reference_set = None
        
        # Find the FASTA file for the top header
        for header, path in header_dict.items(): #headers in each FASTA set a reference.
            if self.top_header in header:
                self.top_fasta = path
                self.reference_set = os.path.basename(os.path.dirname(self.top_fasta))
                break
                
        print(f'\nSample: {self.sample_name}\nTop Sourmash Finding: {self.top_header_found} \nReference Set: {self.reference_set} \nTop reference that is automatically available: {self.top_fasta}\n')
        
        if not self.top_fasta:
            print(f'Reference not found, Must force reference\n')
            self.top_fasta_header = 'Reference not found: Must force reference'
        else:
            try:
                with open(self.top_fasta) as f:
                    self.top_fasta_header = f.readline().strip()
            except Exception as e:
                print(f"Error reading top FASTA file: {str(e)}")
                self.top_fasta_header = 'Error reading reference file'
        
        # Create sourmash directory and move results
        dir = 'sourmash'
        if not os.path.exists(dir):
            os.makedirs(dir)
            
        if os.path.exists(f'{self.sample_name}_search.csv'):
            shutil.move(f'{self.sample_name}_search.csv', dir)
            
        if os.path.exists(f'{self.FASTQ_R1}.sig'):
            os.remove(f'{self.FASTQ_R1}.sig')
            
        print("#############\n")

    def latex(self, tex):
        blast_banner = Banner("Sourmash Sequence Similarity")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{l|l}', file=tex)
        print(r'Similarity & ID \\', file=tex)
        print(r'\hline', file=tex)
        count=0
        
        if not hasattr(self, 'sourmash_df') or self.sourmash_df.empty:
            print('Sourmash - No Data Output & Sourmash - No Data Output \\\\', file=tex)
        else:
            try:
                for row in self.sourmash_df.itertuples():
                    count+=1
                    if count <= 10:
                        percentage = f'{row[1]:.1%}'
                        name = row[4].replace("_", "\_") if isinstance(row[2], str) else "Invalid Name"
                        print(percentage.replace("%", "\%") + ' & ' + name + ' \\\\', file=tex)
                        print(r'\hline', file=tex)
            except Exception as e:
                print(f'Sourmash - Error processing data: {str(e)} & Error \\\\', file=tex)
                print(r'\hline', file=tex)
                
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)

    def excel(self, excel_dict):
        # Add reference set to Excel dictionary
        excel_dict[f'Found_Reference_Set'] = f'{self.reference_set if self.reference_set else "None"}'


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''

    ---------------------------------------------------------
    Usage:
    vsnp3_best_reference_sourmash.py -r1 *_R1*.fastq.gz

    # Validation Notes:
     * Genomes that will automatically be found as best reference *
    2021-05-06 Sourmash testing.
    Brucella_abortus1 checked good
    Brucella_abortus3 checked good
    Brucella_canis checked good
    Brucella_ceti1 checked good
    Brucella_ceti2 checked good
    Brucella_melitensis-bv1 checked good
    Brucella_melitensis-bv1b checked good
    Brucella_melitensis-bv2 checked good
    Brucella_melitensis-bv3 checked good
    Brucella_neotomae: Finds NC_017251 (suis1), Bceti1Cudo and neotomae at various same or wrong similarity
    Brucella_ovis checked good
    Brucella_suis1 checked good
    Brucella_suis2 checked good
    Brucella_suis3 checked good
    Brucella_suis4: Finds canis and suis4 at the same similarity, sometimes canis more.

    Mycobacterium_AF2122 checked
        # 2021-05-06 Current database does not have caprae, orygis, these were test as what-ifs.
        before adding caprae just barely fell into AF2122, orygis was a tie
        adding caprae finds caprae and still bovis and h37 correct
        adding orygis finds orygis and still bovis and h37 correct
    Mycobacterium_H37 checked good
    para-CP033688 checked good
    para-NC002944 checked good

    NC_045512_wuhan-hu-1 checked good

    # Build database
    for i in *fasta; do sourmash sketch dna $i --name-from-first; done
    sourmash index ref_db.sbt.zip ./*.sig # this .zip will be placed as a dependency file, cp ref_db.sbt.zip ~/git/gitlab/dev_stuber/dependencies
    # Prepare read
    sourmash sketch dna *_R1*.fastq.gz
    # Search 
    sourmash search *_R1*.fastq.gz.sig ../sourmash/ref_db.sbt.zip -o sourmash_findings.csv --threshold=0.001

    # Isolates used to build 2021-08-04 database
    # Isolates saved to /Users/tstuber/OneDrive\ -\ USDA/vSNP/vsnp3/validation/sourmash_vsnp_2021-08-04_genomes

    # TB complex lineage information
    # Best reference will direct to either AF2122 or H37
    NZ_CP041790 Lineage 1
    NZ_CP022014 Lineage 2
    NZ_CP041869 Lineage 2
    NZ_CP041875 Lineage 3
    CP023623 Lineage 4
    CP023635 Lineage 4
    NZ_CP041803 Lineage 7

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=False, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    best_reference = Best_Reference(SAMPLE_NAME=args.SAMPLE_NAME, FASTQ_R1=args.FASTQ_R1, debug=args.debug)
    best_reference.run()

    #Latex report
    latex_report = Latex_Report(best_reference.sample_name)
    best_reference.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(best_reference.sample_name)
    best_reference.excel(excel_stats.excel_dict)
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