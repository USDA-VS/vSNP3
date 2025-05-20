#!/usr/bin/env python

import os
import re
import sys
import shutil
import argparse
import textwrap
import numpy as np
import pandas as pd
import multiprocessing

# Python 3.12 compatibility - move set_start_method inside if __name__ block
# to avoid issues with multiprocessing spawning

from krona_lca_all import force_tax_number

__version__ = "3.30"

class Kraken2_Identification:
    ''' 
    Assemble reads using Spades assembler.
    Paired or single reads
    '''

    def __init__(self, **kwargs):
        '''
        See -h
        '''
        # Internal custom database https://github.com/stuber/instructions/blob/master/kraken2_database.md
        self.db = kwargs.get('database', None) #path to database directory
        # External custom databases https://benlangmead.github.io/aws-indexes/k2
        self.threads = multiprocessing.cpu_count() - 2
        self.FASTA = kwargs.get('FASTA', None)
        self.FASTQ_R1 = kwargs.get('FASTQ_R1', None)
        self.FASTQ_R2 = kwargs.get('FASTQ_R2', None)
        self.directory = kwargs.get('directory', None)
        FASTQ_list = [self.FASTQ_R1, self.FASTQ_R2]
        FASTQ_list = [x for x in FASTQ_list if x is not None]  # remove None when single read
        db_contents = kwargs.get('db_contents', False)
        if db_contents:
            with open(f'{self.db}/README', 'r') as opened_file:
                for line in opened_file:
                    print(f'{line.strip()}')
            print(f'\nDatabase location: {self.db}\n')
        if self.FASTA and self.FASTQ_R1:
            print(f'### Error: Can only provide FASTA or FASTQ, not both file types at the same time')
            sys.exit(0)
        if self.FASTA:
            self.sample_name = re.sub('[._].*', '', self.FASTA)
        else:
            self.sample_name = re.sub('[._].*', '', FASTQ_list[0])
        self.cwd = os.getcwd()
        self.FASTQ_list = FASTQ_list

        
    def kraken2_run(self):
        db = self.db
        threads = self.threads
        sample_name = self.sample_name
        FASTQ_list = self.FASTQ_list
        FASTA = self.FASTA
        cwd = self.cwd
        print(f'Kraken2 Running...')
        if len(FASTQ_list) == 2:
            os.system(f'kraken2 --db {db} --threads {threads} --paired {FASTQ_list[0]} {FASTQ_list[1]} --output {sample_name}-outputkraken.txt --report {sample_name}-reportkraken.txt')
        elif len(FASTQ_list) == 1:
            os.system(f'kraken2 --db {db} --threads {threads} {FASTQ_list[0]} --output {sample_name}-outputkraken.txt --report {sample_name}-reportkraken.txt')
        elif FASTA:
            os.system(f'kraken2 --db {db} --threads {threads} {FASTA} --output {sample_name}-outputkraken.txt --report {sample_name}-reportkraken.txt')
        else:
            print(f'\n### Error: Missing read files\n')
            sys.exit(0)

        if os.path.exists(f'{cwd}/{sample_name}-outputkraken.txt'):
            output = f'{cwd}/{sample_name}-outputkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(0)
        
        if os.path.exists(f'{cwd}/{sample_name}-reportkraken.txt'):
            report = f'{cwd}/{sample_name}-reportkraken.txt'
        else:
            print(f'\n### Error: Kraken report did not complete')
            sys.exit(0)

        if self.directory:
            os.makedirs(self.directory, exist_ok=True)  # Using exist_ok for Python 3.12 compatibility
            shutil.move(report, self.directory)
            shutil.move(output, self.directory)
            report = f'{cwd}/{self.directory}/{sample_name}-reportkraken.txt'
            output = f'{cwd}/{self.directory}/{sample_name}-outputkraken.txt'
            log_file = open("kraken_log.txt", "a")
            try:
                log_file.write(f'DB used: {os.readlink(self.db)}')
            except OSError:
                log_file.write(f'DB used: {self.db}')
            log_file.close()
            shutil.move("kraken_log.txt", self.directory)
            return report, output
        else:
            return report, output

    def krona_make_graph(self, report, output):
        # Output will be: kronaInput.txt
        # Two column file will contain read header and taxid
        force_tax_number(output)
        os.system(f'ktImportTaxonomy kronaInput.txt')
        os.rename(f'taxonomy.krona.html', f'{self.sample_name}-taxonomy.krona.html')
        try:
            shutil.rmtree(f'taxonomy.krona.html.files')
        except FileNotFoundError:
            pass
        os.remove(f'kronaInput.txt')
        if os.path.exists(f'{self.cwd}/{self.sample_name}-taxonomy.krona.html'):
            krona_html = f'{self.cwd}/{self.sample_name}-taxonomy.krona.html'
        else:
            print(f'\n### Error: Krona HTML did not complete')
            sys.exit(0)

        if self.directory:
            shutil.move(krona_html, self.directory)
            krona_html = f'{self.cwd}/{self.directory}/{self.sample_name}-taxonomy.krona.html'
            return krona_html
        else:
            return krona_html

    def bracken(self, report, output):
        os.system(f'bracken -d {self.db} -i {report} -o {self.sample_name}-bracken.txt -r 250')
        
        # Updated pandas usage for better compatibility
        df = pd.read_csv(f'{self.sample_name}-bracken.txt', sep='\t')
        
        # Use pandas' ExcelWriter with engine specification for compatibility
        with pd.ExcelWriter(f'{self.sample_name}-bracken.xlsx', engine='openpyxl') as writer:
            df.to_excel(writer, index=False)
            
        os.remove(f'{self.sample_name}-bracken.txt')
        self.bracken_excel = f'{os.getcwd()}/{self.sample_name}-bracken.xlsx'
        if self.directory:
            shutil.move(f'{self.sample_name}-bracken.xlsx', self.directory)
            self.bracken_excel = f'{os.getcwd()}/{self.directory}/{self.sample_name}-bracken.xlsx'

if __name__ == "__main__": # execute if directly access by the interpreter
    # Set multiprocessing start method here for Python 3.12 compatibility
    multiprocessing.set_start_method('spawn', True)

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        Provide either a single FASTA file, single FASTQ or Paired files.
        Usage:
            kraken2_run.py -r1 *_R1*fastg.gz
            kraken2_run.py -r1 *_R1*fastg.gz -r2 *_R2*fastq.gz -d
            kraken2_run.py -f *fasta

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--FASTA', action='store', dest='FASTA', required=False, help='Provide FASTA file')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, help='Provide R1 FASTQ gz file, or single read')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Provide R2 FASTQ gz file')
    parser.add_argument('-d', '--directory', action='store', dest='directory', required=False, default="kraken2", help='Put output to directory')
    parser.add_argument('-c', '--db_contents', action='store_true', dest='db_contents', help='Show contents of DB by printing README')
    parser.add_argument('--database', required=True, action='store', dest='database', help='Absolute path to database directory')
    args = parser.parse_args()

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    kraken2 = Kraken2_Identification(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, directory=args.directory, db_contents=args.db_contents, database=args.database)
    report, output = kraken2.kraken2_run()
    krona_html = kraken2.krona_make_graph(report, output)
    # kraken2.bracken(report, output)

    print('done')
# Created March 2020 by Tod Stuber