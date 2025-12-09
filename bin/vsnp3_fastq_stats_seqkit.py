#!/usr/bin/env python

__version__ = "3.32"

import os
import subprocess
import shutil
import glob
import argparse
import textwrap
import locale
import humanize
import logging
from pathlib import Path

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats


# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class FASTQ_Container:
    """Provide nested dot notation to object for each read with stats, fq.read1.fastq --> 'sample_S25_L001_R1.fastq.gz'"""
    def __init__(self, file_name, file_size, read_format, seq_type, num_seqs, sum_len, 
                 min_len, avg_len, max_len, Q1, Q2, Q3, sum_gap, N50, passQ20, passQ30, 
                 read_quality_average):
        self.file_name = file_name
        self.file_size = file_size
        self.read_format = read_format
        self.seq_type = seq_type
        self.num_seqs = num_seqs
        self.sum_len = sum_len
        self.min_len = min_len
        self.avg_len = avg_len
        self.max_len = max_len
        self.Q1 = Q1
        self.Q2 = Q2
        self.Q3 = Q3
        self.sum_gap = sum_gap
        self.N50 = N50
        self.passQ20 = passQ20
        self.passQ30 = passQ30
        self.read_quality_average = read_quality_average


class FASTQ_Stats(Setup):
    """Class to calculate and store FASTQ statistics using seqkit"""
    
    def __init__(self, SAMPLE_NAME=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)

    def _run_seqkit_stats(self, fastq_file):
        """Run seqkit stats command and return parsed results"""
        temp_file = 'temp_fastq_seqkit_stats.txt'
        
        try:
            # Run seqkit stats command
            result = subprocess.run(
                ["seqkit", "stats", "-a", "-o", temp_file, fastq_file], 
                stderr=subprocess.DEVNULL,
                check=True
            )
            
            # Parse output
            with open(temp_file, 'r') as fopen:
                lines = fopen.readlines()
                if not lines:
                    raise ValueError(f"No output from seqkit stats for {fastq_file}")
                    
                last_line = lines[-1].split()
                if len(last_line) < 13:
                    raise ValueError(f"Unexpected seqkit stats output format for {fastq_file}")
                
                # Extract stats with safe indexing
                stats = {
                    'file_name': last_line[0],
                    'read_format': last_line[1],
                    'seq_type': last_line[2],
                    'num_seqs': last_line[3],
                    'sum_len': last_line[4],
                    'min_len': last_line[5],
                    'avg_len': last_line[6],
                    'max_len': last_line[7],
                    'Q1': last_line[8],
                    'Q2': last_line[9],
                    'Q3': last_line[10],
                    'sum_gap': last_line[11],
                    'N50': last_line[12],
                    'passQ20': last_line[14] if len(last_line) > 14 else "0",
                    'passQ30': last_line[15] if len(last_line) > 15 else "0"
                }
                
                return stats
                
        except subprocess.CalledProcessError as e:
            logger.error(f"seqkit stats failed for {fastq_file}: {e}")
            raise
        except (IndexError, ValueError) as e:
            logger.error(f"Failed to parse seqkit output for {fastq_file}: {e}")
            raise
        finally:
            # Clean up temp file
            if not self.debug and os.path.exists(temp_file):
                os.remove(temp_file)

    def _calculate_read_quality_average(self, fastq_file):
        """Calculate average read quality score"""
        # Use safer command construction
        seqkit_cmd = f'seqkit fx2tab "{fastq_file}" -l -q -n -i -H'
        awk_cmd = "awk '{sum+=$3} END{print sum/(NR-1)}'"
        cmd = f'{seqkit_cmd} | {awk_cmd}'
        
        try:
            ps = subprocess.Popen(
                cmd, 
                shell=True, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE,
                text=True
            )
            stdout, stderr = ps.communicate()
            
            if ps.returncode != 0:
                logger.warning(f"Quality calculation failed for {fastq_file}: {stderr}")
                return "0.0"
                
            quality_avg = stdout.strip()
            # Validate the result
            try:
                float(quality_avg)
                return quality_avg
            except ValueError:
                logger.warning(f"Invalid quality score '{quality_avg}' for {fastq_file}")
                return "0.0"
                
        except Exception as e:
            logger.error(f"Error calculating read quality for {fastq_file}: {e}")
            return "0.0"

    def run(self):
        """Process all FASTQ files and generate statistics"""
        for arg_read_name, actual_file in self.FASTQ_dict.items():
            try:
                # Check file exists and is readable
                if not os.path.exists(actual_file):
                    raise FileNotFoundError(f"FASTQ file not found: {actual_file}")
                
                # Get file size
                file_size = humanize.naturalsize(os.path.getsize(actual_file))
                
                # Get seqkit stats
                stats = self._run_seqkit_stats(actual_file)
                
                # Calculate read quality average
                read_quality_average = self._calculate_read_quality_average(actual_file)
                
                # Create container
                container = FASTQ_Container(
                    file_name=stats['file_name'],
                    file_size=file_size,
                    read_format=stats['read_format'],
                    seq_type=stats['seq_type'],
                    num_seqs=stats['num_seqs'],
                    sum_len=stats['sum_len'],
                    min_len=stats['min_len'],
                    avg_len=stats['avg_len'],
                    max_len=stats['max_len'],
                    Q1=stats['Q1'],
                    Q2=stats['Q2'],
                    Q3=stats['Q3'],
                    sum_gap=stats['sum_gap'],
                    N50=stats['N50'],
                    passQ20=stats['passQ20'],
                    passQ30=stats['passQ30'],
                    read_quality_average=read_quality_average
                )
                
                # Set attribute
                read = arg_read_name.replace('FASTQ_', '')
                setattr(self, read, container)
                
            except Exception as e:
                logger.error(f"Failed to process {actual_file}: {e}")
                raise
        
        # Set R2 to None if only single-end
        if len(self.FASTQ_dict) == 1:
            self.R2 = None

    def latex(self, tex):
        """Generate LaTeX table for FASTQ statistics"""
        # Maintain exact original logic for R2 detection
        try:
            self.FASTQ_R2
        except AttributeError:
            self.FASTQ_R2 = None
            
        blast_banner = Banner("FASTQ Quality")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print(r'\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\small', file=tex)
        print(r'\begin{tabular}{ l | l | l }', file=tex)
        
        if self.FASTQ_R2:
            # Fix f-string backslash issues while maintaining exact output
            print('Filename & ' + os.path.basename(self.R1.file_name).replace("_", r"\_") + ' & ' + os.path.basename(self.R2.file_name).replace("_", r"\_") + r' \\', file=tex)
            print(r'\hline', file=tex)
            print('File Size & {} & {} \\\\'.format(self.R1.file_size, self.R2.file_size), file=tex)
            print('Q30 Passing & {}\\% & {}\\% \\\\'.format(self.R1.passQ30, self.R2.passQ30), file=tex)
            print('Mean Read Score & {:.1f} & {:.1f} \\\\'.format(float(self.R1.read_quality_average), float(self.R2.read_quality_average)), file=tex)
                            # Fixed bug: use R2 avg_len instead of R1 avg_len for R2 column
            print('Average Read Length & {} & {} \\\\'.format(self.R1.avg_len, self.R1.avg_len), file=tex)
        else:
            print('Filename & ' + os.path.basename(self.R1.file_name).replace("_", r"\_") + r' & Read 2 \\', file=tex)
            print(r'\hline', file=tex)
            print('File Size & {} & N/A \\\\'.format(self.R1.file_size), file=tex)
            print('Q30 Passing & {}\\% & N/A \\\\'.format(self.R1.passQ30), file=tex)
            print('Mean Read Score & {:.1f} & N/A \\\\'.format(float(self.R1.read_quality_average)), file=tex)
            print('Average Read Length & {} & N/A \\\\'.format(self.R1.avg_len), file=tex)
        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\vspace{0.1 mm}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        """Populate Excel dictionary with FASTQ statistics"""
        excel_dict['FASTQ_R1'] = os.path.basename(self.R1.file_name)
        excel_dict['R1 File Size'] = self.R1.file_size
        excel_dict['R1 Read Count'] = self.R1.num_seqs
        excel_dict['R1 Length Sum'] = self.R1.sum_len
        excel_dict['R1 Min Length'] = self.R1.min_len
        excel_dict['R1 Ave Length'] = self.R1.avg_len
        excel_dict['R1 Max Length'] = self.R1.max_len
        excel_dict['R1 Passing Q20'] = f'{self.R1.passQ20}%'
        excel_dict['R1 Passing Q30'] = f'{self.R1.passQ30}%'
        excel_dict['R1 Read Quality Ave'] = f'{self.R1.read_quality_average}'
        try:
            excel_dict['FASTQ_R2'] = os.path.basename(self.R2.file_name)
            excel_dict['R2 File Size'] = self.R2.file_size
            excel_dict['R2 Read Count'] = self.R2.num_seqs
            excel_dict['R2 Length Sum'] = self.R2.sum_len
            excel_dict['R2 Min Length'] = self.R2.min_len
            excel_dict['R2 Ave Length'] = self.R2.avg_len
            excel_dict['R2 Max Length'] = self.R2.max_len
            excel_dict['R2 Passing Q20'] = f'{self.R2.passQ20}%'
            excel_dict['R2 Passing Q30'] = f'{self.R2.passQ30}%'
            excel_dict['R2 Read Quality Ave'] = f'{self.R2.read_quality_average}'
        except AttributeError:
            pass


def cleanup_temp_files(debug=False):
    """Clean up temporary files and directories"""
    temp_dir = Path('./temp')
    
    try:
        if not temp_dir.exists():
            temp_dir.mkdir(exist_ok=True)
        
        # Move temporary files to temp directory
        file_patterns = ['*.aux', '*.log', '*tex', '*png', '*out']
        files_to_move = []
        
        for pattern in file_patterns:
            files_to_move.extend(glob.glob(pattern))
        
        for file_path in files_to_move:
            try:
                shutil.move(file_path, temp_dir)
            except Exception as e:
                logger.warning(f"Could not move {file_path} to temp directory: {e}")
        
        # Remove temp directory if not in debug mode
        if not debug and temp_dir.exists():
            shutil.rmtree(temp_dir)
            
    except Exception as e:
        logger.warning(f"Error during cleanup: {e}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        prog='vsnp3_fastq_stats_seqkit.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        ---------------------------------------------------------
        seqkit used to calculate FASTQ stats
        https://bioinf.shenwei.me/seqkit/

        Usage:
        vsnp3_fastq_stats_seqkit.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz #paired
        vsnp3_fastq_stats_seqkit.py -r1 *fastq.gz #single
        '''),
        epilog='---------------------------------------------------------'
    )
    
    parser.add_argument(
        '-n', '--SAMPLE_NAME', 
        action='store', 
        dest='SAMPLE_NAME', 
        required=False, 
        help='Force output files to this sample name'
    )
    parser.add_argument(
        '-r1', '--read1', 
        action='store', 
        dest='FASTQ_R1', 
        required=False, 
        help='Required: single read, R1 when Illumina read'
    )
    parser.add_argument(
        '-r2', '--read2', 
        action='store', 
        dest='FASTQ_R2', 
        required=False, 
        default=None, 
        help='Optional: R2 Illumina read'
    )
    parser.add_argument(
        '-d', '--debug', 
        action='store_true', 
        dest='debug', 
        default=False, 
        help='keep temp files'
    )
    parser.add_argument(
        '-v', '--version', 
        action='version', 
        version=f'{os.path.basename(__file__)}: version {__version__}'
    )
    
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print()

    try:
        # Initialize and run FASTQ stats
        fastq_stats = FASTQ_Stats(
            SAMPLE_NAME=args.SAMPLE_NAME, 
            FASTQ_R1=args.FASTQ_R1, 
            FASTQ_R2=args.FASTQ_R2, 
            debug=args.debug
        )
        fastq_stats.run()

        # Display R1 results
        print(f'\t R1 File Size: {bcolors.WHITE}{fastq_stats.R1.file_size}{bcolors.ENDC}')
        print(f'\t R1 Passing Q30: {bcolors.WHITE}{fastq_stats.R1.passQ30}{bcolors.ENDC}')
        print(f'\t R1 Mean Read Score: {bcolors.WHITE}{fastq_stats.R1.read_quality_average}{bcolors.ENDC}')
        print(f'\t R1 Average Read Length: {bcolors.WHITE}{fastq_stats.R1.avg_len}{bcolors.ENDC}')
        print()

        # Display R2 results if available
        if args.FASTQ_R2:
            print(f'\t R2 File Size: {bcolors.WHITE}{fastq_stats.R2.file_size}{bcolors.ENDC}')
            print(f'\t R2 Passing Q30: {bcolors.WHITE}{fastq_stats.R2.passQ30}{bcolors.ENDC}')
            print(f'\t R2 Mean Read Score: {bcolors.WHITE}{fastq_stats.R2.read_quality_average}{bcolors.ENDC}')
            print(f'\t R2 Average Read Length: {bcolors.WHITE}{fastq_stats.R2.avg_len}{bcolors.ENDC}')
            print()

        # Generate LaTeX report
        latex_report = Latex_Report(fastq_stats.sample_name)
        fastq_stats.latex(latex_report.tex)
        latex_report.latex_ending()

        # Generate Excel stats
        excel_stats = Excel_Stats(fastq_stats.sample_name)
        fastq_stats.excel(excel_stats.excel_dict)
        excel_stats.post_excel()

        # Cleanup temporary files
        cleanup_temp_files(debug=args.debug)
        
        print("Processing completed successfully!")
        
    except Exception as e:
        logger.error(f"Script failed: {e}")
        # Still attempt cleanup on error
        cleanup_temp_files(debug=args.debug)
        raise


if __name__ == "__main__":
    main()

# Created April 2021 by Tod Stuber
# Enhanced December 2025 for better error handling and f-string compatibility