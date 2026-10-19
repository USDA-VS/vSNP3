#!/usr/bin/env python

from vsnp3_version import __version__

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


# seqkit stats -a -T column header -> the name used here.  Looked up by name:
# the Q20 and Q30 columns used to be taken from positions 14 and 15 of the
# whitespace-split line, which assumes the N50_num column seqkit added in 2.4;
# on an older seqkit every sample read AvgQual as its Q30 and was Questionable.
SEQKIT_COLUMNS = {
    'file': 'file_name', 'format': 'read_format', 'type': 'seq_type',
    'num_seqs': 'num_seqs', 'sum_len': 'sum_len', 'min_len': 'min_len',
    'avg_len': 'avg_len', 'max_len': 'max_len', 'Q1': 'Q1', 'Q2': 'Q2', 'Q3': 'Q3',
    'sum_gap': 'sum_gap', 'N50': 'N50', 'Q20(%)': 'passQ20', 'Q30(%)': 'passQ30',
}


def parse_seqkit_stats_table(text, fastq_file=''):
    '''
    The last data row of `seqkit stats -a` output as {name: string as printed},
    columns found by header name.  Raises ValueError when a needed column is
    missing, rather than reading a neighbouring one.

    The default (human) output is parsed, not -T.  The workbook has always held
    the default strings (thousands separators, "93.61" for Q20(%)), and the
    tab-separated form is not the same numbers formatted differently: seqkit
    2.13 prints Q20(%) and Q30(%) there rounded to whole numbers ("94", "90")
    while its default output keeps two decimals.  Parsing the default output by
    header name keeps every stored string as it was and survives a column being
    added, as N50_num was in seqkit 2.4.  Values never contain spaces (vSNP3
    refuses paths with whitespace), so splitting on runs of spaces is safe.
    '''
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError(f"No output from seqkit stats for {fastq_file}")
    header = lines[0].split()
    row = lines[-1].split()
    if len(row) != len(header):
        raise ValueError(f"Unexpected seqkit stats output format for {fastq_file}: "
                         f"{len(header)} header columns, {len(row)} values")
    by_header = dict(zip(header, row))
    missing = [column for column in SEQKIT_COLUMNS if column not in by_header]
    if missing:
        raise ValueError(f"seqkit stats output for {fastq_file} lacks column(s) "
                         f"{', '.join(missing)}; seqkit >= 2.1 with -a is required")
    return {name: by_header[column] for column, name in SEQKIT_COLUMNS.items()}


class FASTQ_Stats(Setup):
    """Class to calculate and store FASTQ statistics using seqkit"""

    def __init__(self, SAMPLE_NAME=None, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)

    def _run_seqkit_stats(self, fastq_file):
        """Run seqkit stats command and return parsed results"""
        # Named for the sample so two runs in one directory cannot read each
        # other's table.
        temp_file = f'{self.sample_name}_temp_fastq_seqkit_stats.txt'

        try:
            # Default output, read by header name (see parse_seqkit_stats_table
            # for why not -T).
            result = subprocess.run(
                ["seqkit", "stats", "-a", "-o", temp_file, fastq_file],
                stderr=subprocess.DEVNULL,
                check=True
            )

            # Parse output
            with open(temp_file, 'r', encoding='utf-8', errors='replace') as fopen:
                return parse_seqkit_stats_table(fopen.read(), fastq_file)

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
        file_patterns = ['*.log', '*out']
        files_to_move = []
        
        for pattern in file_patterns:
            files_to_move.extend(sorted(glob.glob(pattern)))
        
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