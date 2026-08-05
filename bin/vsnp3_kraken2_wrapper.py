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
import subprocess
from pathlib import Path

# Python 3.12 compatibility - move set_start_method inside if __name__ block
# to avoid issues with multiprocessing spawning

from krona_lca_all import force_tax_number

from vsnp3_version import __version__

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
            readme_path = Path(self.db) / 'README'
            try:
                with open(readme_path, 'r') as opened_file:
                    for line in opened_file:
                        print(line.strip())
                print(f'\nDatabase location: {self.db}\n')
            except FileNotFoundError:
                print(f'Warning: README file not found at {readme_path}')
            except Exception as e:
                print(f'Error reading README file: {e}')
                
        if self.FASTA and self.FASTQ_R1:
            print('### Error: Can only provide FASTA or FASTQ, not both file types at the same time')
            sys.exit(1)
            
        if self.FASTA:
            self.sample_name = re.sub(r'[._].*', '', Path(self.FASTA).name)
        elif FASTQ_list:
            self.sample_name = re.sub(r'[._].*', '', Path(FASTQ_list[0]).name)
        else:
            print('### Error: No input files provided')
            sys.exit(1)
            
        self.cwd = os.getcwd()
        self.FASTQ_list = FASTQ_list

    def _run_command(self, cmd_args, description="Command"):
        """Safely run a command using subprocess instead of os.system"""
        try:
            print(f'{description} running...')
            result = subprocess.run(cmd_args, check=True, capture_output=True, text=True)
            if result.stdout:
                print(result.stdout)
            return True
        except subprocess.CalledProcessError as e:
            print(f'Error running {description}: {e}')
            if e.stderr:
                print(f'Error details: {e.stderr}')
            return False
        except Exception as e:
            print(f'Unexpected error running {description}: {e}')
            return False
        
    def kraken2_run(self):
        db = self.db
        threads = str(self.threads)
        sample_name = self.sample_name
        FASTQ_list = self.FASTQ_list
        FASTA = self.FASTA
        cwd = self.cwd
        
        output_file = f'{sample_name}-outputkraken.txt'
        report_file = f'{sample_name}-reportkraken.txt'
        
        print('Kraken2 Running...')
        
        # Build command safely
        cmd_args = ['kraken2', '--db', db, '--threads', threads]
        
        if len(FASTQ_list) == 2:
            cmd_args.extend(['--paired', FASTQ_list[0], FASTQ_list[1]])
        elif len(FASTQ_list) == 1:
            cmd_args.append(FASTQ_list[0])
        elif FASTA:
            cmd_args.append(FASTA)
        else:
            print('\n### Error: Missing read files\n')
            sys.exit(1)
            
        cmd_args.extend(['--output', output_file, '--report', report_file])
        
        # Run kraken2 command safely
        if not self._run_command(cmd_args, "Kraken2"):
            print('\n### Error: Kraken2 command failed')
            sys.exit(1)

        # Check output files exist
        output_path = Path(cwd) / output_file
        report_path = Path(cwd) / report_file
        
        if not output_path.exists():
            print('\n### Error: Kraken output file was not created')
            sys.exit(1)
            
        if not report_path.exists():
            print('\n### Error: Kraken report file was not created')
            sys.exit(1)

        output = str(output_path)
        report = str(report_path)

        if self.directory:
            try:
                os.makedirs(self.directory, exist_ok=True)
                
                dest_report = Path(self.directory) / report_file
                dest_output = Path(self.directory) / output_file
                
                shutil.move(str(report_path), str(dest_report))
                shutil.move(str(output_path), str(dest_output))
                
                report = str(dest_report)
                output = str(dest_output)
                
                # Create log file
                log_content = f'DB used: '
                try:
                    log_content += os.readlink(self.db)
                except OSError:
                    log_content += self.db
                    
                log_path = Path(self.directory) / "kraken_log.txt"
                try:
                    with open(log_path, "a") as log_file:
                        log_file.write(log_content + '\n')
                except Exception as e:
                    print(f'Warning: Could not write log file: {e}')
                    
            except Exception as e:
                print(f'Error organizing output files: {e}')
                sys.exit(1)

        return report, output

    def krona_make_graph(self, report, output):
        """Create Krona HTML visualization"""
        try:
            # Output will be: kronaInput.txt
            # Two column file will contain read header and taxid
            force_tax_number(output)
            
            # Run ktImportTaxonomy safely
            cmd_args = ['ktImportTaxonomy', 'kronaInput.txt']
            if not self._run_command(cmd_args, "ktImportTaxonomy"):
                print('\n### Error: ktImportTaxonomy command failed')
                sys.exit(1)
                
            # Rename output file
            original_html = 'taxonomy.krona.html'
            new_html = f'{self.sample_name}-taxonomy.krona.html'
            
            if Path(original_html).exists():
                os.rename(original_html, new_html)
            else:
                print('\n### Error: Krona HTML was not created')
                sys.exit(1)
                
            # Clean up files
            files_dir = Path(f'{original_html}.files')
            if files_dir.exists():
                shutil.rmtree(str(files_dir))
                
            krona_input = Path('kronaInput.txt')
            if krona_input.exists():
                os.remove(str(krona_input))

            krona_html_path = Path(self.cwd) / new_html
            if not krona_html_path.exists():
                print('\n### Error: Krona HTML file not found after processing')
                sys.exit(1)
                
            krona_html = str(krona_html_path)

            if self.directory:
                dest_html = Path(self.directory) / new_html
                shutil.move(krona_html, str(dest_html))
                krona_html = str(dest_html)
                
            return krona_html
            
        except Exception as e:
            print(f'Error creating Krona graph: {e}')
            sys.exit(1)

    def bracken(self, report, output):
        """Run Bracken analysis"""
        try:
            bracken_output = f'{self.sample_name}-bracken.txt'
            bracken_excel = f'{self.sample_name}-bracken.xlsx'
            
            # Run bracken safely
            cmd_args = ['bracken', '-d', self.db, '-i', report, '-o', bracken_output, '-r', '250']
            if not self._run_command(cmd_args, "Bracken"):
                print('\n### Error: Bracken command failed')
                return
                
            # Process with pandas
            if Path(bracken_output).exists():
                df = pd.read_csv(bracken_output, sep='\t')
                
                # Use pandas' ExcelWriter with engine specification for compatibility
                with pd.ExcelWriter(bracken_excel, engine='openpyxl') as writer:
                    df.to_excel(writer, index=False)
                    
                os.remove(bracken_output)
                
                self.bracken_excel = str(Path(os.getcwd()) / bracken_excel)
                
                if self.directory:
                    dest_excel = Path(self.directory) / bracken_excel
                    shutil.move(bracken_excel, str(dest_excel))
                    self.bracken_excel = str(dest_excel)
            else:
                print('Warning: Bracken output file not found')
                
        except Exception as e:
            print(f'Error running Bracken analysis: {e}')

if __name__ == "__main__": # execute if directly access by the interpreter
    # Set multiprocessing start method here for Python 3.12 compatibility
    try:
        multiprocessing.set_start_method('spawn', True)
    except RuntimeError:
        # Method already set, continue
        pass

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
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.basename(__file__), __version__))
    
    # parse_args() exits 0 for --help and --version and 2 for a bad argument.  This
    # was wrapped in `except SystemExit: sys.exit(1)`, which turned a successful
    # --help or --version into a failure exit code.
    args = parser.parse_args()

    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)

    try:
        kraken2 = Kraken2_Identification(FASTA=args.FASTA, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, directory=args.directory, db_contents=args.db_contents, database=args.database)
        report, output = kraken2.kraken2_run()
        krona_html = kraken2.krona_make_graph(report, output)
        # kraken2.bracken(report, output)
        
        print('done')
        
    except KeyboardInterrupt:
        print('\nOperation cancelled by user')
        sys.exit(1)
    except Exception as e:
        print(f'Unexpected error: {e}')
        sys.exit(1)

# Created March 2020 by Tod Stuber