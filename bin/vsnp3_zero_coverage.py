#!/usr/bin/env python

__version__ = "3.35"

import os
import re
import shutil
import locale
import argparse
import textwrap
import pandas as pd
from Bio import SeqIO
import subprocess
import sys

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

class Zero_Coverage(Setup):
    ''' 
    Class to identify zero coverage regions in BAM alignments and 
    incorporate these into VCF files
    '''
    def __init__(self, FASTA=None, bam=None, vcf=None, debug=False):

        def count_ac1_positions(vcf_file):
            """Count positions with AC=1 in VCF file"""
            ac1_count = 0
            
            try:
                with open(vcf_file, 'r') as f:
                    for line in f:
                        # Skip header lines
                        if line.startswith('#'):
                            continue
                        
                        # Split the line into fields
                        fields = line.strip().split('\t')
                        
                        # Check if we have enough fields and INFO field contains AC=1
                        if len(fields) >= 8:
                            info = fields[7]
                            if 'AC=1' in info.split(';'):
                                ac1_count += 1
            except (IOError, OSError) as e:
                print(f"Error reading VCF file {vcf_file}: {e}")
                return 0
            
            return ac1_count

        # Validate input files exist
        if not os.path.exists(vcf):
            raise FileNotFoundError(f"VCF file not found: {vcf}")
        if not os.path.exists(bam):
            raise FileNotFoundError(f"BAM file not found: {bam}")
        if not os.path.exists(FASTA):
            raise FileNotFoundError(f"FASTA file not found: {FASTA}")

        self.ac1_count = count_ac1_positions(vcf)
        print(f"Number of positions with AC=1: {self.ac1_count:,}")

        Setup.__init__(self, FASTA=FASTA, debug=debug)
        self.print_run_time('Zero Coverage')
        self.sample_name = re.sub('[_.].*', '', bam)
        self.reference_name = re.sub('[_.].*', '', os.path.basename(FASTA))
        coverage_dict = {}
        
        # Use subprocess with proper error handling
        try:
            result = subprocess.run(['samtools', 'depth', bam], 
                                  capture_output=True, text=True, check=True)
            for line in result.stdout.splitlines():
                parts = line.split('\t')
                if len(parts) >= 3:
                    chrom, position, depth = parts[0], parts[1], parts[2]
                    coverage_dict[chrom + "-" + position] = depth
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools depth: {e}")
            raise
        except FileNotFoundError:
            print("Error: samtools not found. Please ensure samtools is installed and in your PATH.")
            raise
            
        coverage_df = pd.DataFrame.from_dict(coverage_dict, orient='index', columns=["depth"])
        zero_dict = {}
        reference_length = 0
        
        # Parse FASTA with error handling
        try:
            for record in SeqIO.parse(FASTA, "fasta"):
                chrom = record.id
                total_len = len(record.seq)
                reference_length = reference_length + len(record.seq)
                for pos in list(range(1, total_len + 1)):
                    zero_dict[str(chrom) + "-" + str(pos)] = 0
        except Exception as e:
            print(f"Error parsing FASTA file {FASTA}: {e}")
            raise
                
        zero_df = pd.DataFrame.from_dict(zero_dict, orient='index', columns=["depth"])
        
        # Improved pandas merge operation with error handling
        try:
            coverage_df = zero_df.merge(coverage_df, left_index=True, right_index=True, 
                                      how='outer', suffixes=('_x', '_y'))
            coverage_df = coverage_df.drop(columns=['depth_x'])
            coverage_df = coverage_df.rename(columns={'depth_y': 'depth'})
            coverage_df = coverage_df.fillna(0)
            coverage_df['depth'] = coverage_df['depth'].astype(int)
        except Exception as e:
            print(f"Error processing coverage data: {e}")
            raise
        
        total_length = len(coverage_df)
        ave_coverage = coverage_df['depth'].mean()
        zero_df = coverage_df[coverage_df['depth'] == 0]
        total_zero_coverage = len(zero_df)
        
        # Handle division by zero
        if reference_length > 0:
            percent_ref_with_zero_coverage = total_zero_coverage / reference_length * 100
        else:
            percent_ref_with_zero_coverage = 0
            
        print(f'\tPositions with no coverage: {total_zero_coverage:,}, ' +
              f'{bcolors.BLUE}{percent_ref_with_zero_coverage:,.6f}%{bcolors.ENDC} of reference\n')
        
        total_coverage = total_length - total_zero_coverage
        
        # Handle division by zero
        if total_length > 0:
            genome_coverage = float(total_coverage / total_length)
        else:
            genome_coverage = 0.0
        
        # Use context manager for file operations with better error handling
        try:
            vcf_df = pd.read_csv(vcf, sep='\t', header=None, 
                              names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"], 
                              comment='#')
        except pd.errors.EmptyDataError:
            # Handle empty VCF file
            vcf_df = pd.DataFrame(columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"])
        except Exception as e:
            print(f"Error reading VCF file {vcf}: {e}")
            raise
            
        # Use boolean indexing more efficiently with error handling
        try:
            # Handle potential missing columns or data type issues
            if not vcf_df.empty:
                good_snp_mask = ((vcf_df['ALT'].astype(str).str.len() == 1) & 
                               (vcf_df['REF'].astype(str).str.len() == 1) & 
                               (pd.to_numeric(vcf_df['QUAL'], errors='coerce') > 150))
                good_snp_count = good_snp_mask.sum()
            else:
                good_snp_count = 0
        except Exception as e:
            print(f"Error calculating good SNP count: {e}")
            good_snp_count = 0
            
        # Handle division by zero
        if reference_length > 0:
            percent_ref_with_good_snp_count = good_snp_count / reference_length * 100
        else:
            percent_ref_with_good_snp_count = 0
            
        zero_coverage_vcf = f'{self.sample_name}_zc.vcf'
        
        if total_zero_coverage > 0:
            try:
                # Write header
                with open('v_header.csv', 'w') as header_out:
                    with open(vcf, 'r') as fff:
                        for line in fff:
                            if re.search('^#', line):
                                print(line.strip(), file=header_out)
                
                # Create a new dataframe for zero coverage positions
                zero_positions_list = []
                for idx in zero_df.index:
                    try:
                        chrom, pos = idx.rsplit('-', 1)
                        zero_positions_list.append({
                            'CHROM': chrom,
                            'POS': int(pos),
                            'ID': '.',
                            'REF': 'N',
                            'ALT': '.',
                            'QUAL': '.',
                            'FILTER': '.',
                            'INFO': '.',
                            'FORMAT': 'GT',
                            'Sample': './.'
                        })
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing position {idx}: {e}")
                        continue
                
                # Create a dataframe from the list of dictionaries
                if zero_positions_list:
                    zero_positions_df = pd.DataFrame(zero_positions_list)
                    
                    # Combine the original VCF with zero coverage positions
                    combined_df = pd.concat([vcf_df, zero_positions_df], ignore_index=True)
                    combined_df = combined_df.sort_values(['CHROM', 'POS'])
                    
                    # Write the combined data to the output file
                    combined_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
                else:
                    # If no zero positions, just write the original VCF data
                    vcf_df.to_csv('v_annotated_body.csv', sep='\t', header=False, index=False)
                
                # Use context manager for safer file operations
                with open(zero_coverage_vcf, "wb") as outfile:
                    for cf in ['v_header.csv', 'v_annotated_body.csv']:
                        if os.path.exists(cf):
                            with open(cf, "rb") as infile:
                                outfile.write(infile.read())
                
                # Clean up temporary files
                for temp_file in ['v_header.csv', 'v_annotated_body.csv']:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                        
            except Exception as e:
                print(f"Error creating zero coverage VCF: {e}")
                # Fallback to copying original VCF
                shutil.copyfile(vcf, zero_coverage_vcf)
        else:
            shutil.copyfile(vcf, zero_coverage_vcf)
            
        # Store all attributes
        self.bam = bam
        self.reference_length = reference_length
        self.zero_coverage_vcf = zero_coverage_vcf
        self.good_snp_count = good_snp_count
        self.percent_ref_with_good_snp_count = percent_ref_with_good_snp_count
        self.ave_coverage = ave_coverage
        self.genome_coverage = genome_coverage
        self.total_zero_coverage = total_zero_coverage
        self.percent_ref_with_zero_coverage = percent_ref_with_zero_coverage

    def latex(self, tex):
        """Generate LaTeX output with proper escaping"""
        blast_banner = Banner(f'Coverage against {self.reference_name}')
        
        # Use regular string concatenation to avoid f-string backslash issues
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print(r'\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{ l | l | l | l | l | l | l | l }', file=tex)
        
        # Build the header string without f-string backslashes
        header_line = ('BAM File & Reference Length & Genome with Coverage & Average Depth & '
                      'No Coverage Bases & Ambiguous SNPs & Quality SNPs ' + r'\\')
        print(header_line, file=tex)
        print(r'\hline', file=tex)
        
        # Escape underscores in BAM filename for LaTeX
        bam_escaped = self.bam.replace('_', r'\_')
        
        # Build the data line without f-string backslashes
        data_line = (f'{bam_escaped} & {self.reference_length:,} & {(self.genome_coverage*100):,.2f}' +
                    r'\% & ' + f'{self.ave_coverage:,.1f}X & {self.total_zero_coverage:,} & ' +
                    f'{self.ac1_count:,} & {self.good_snp_count:,} ' + r'\\')
        print(data_line, file=tex)
        
        print(r'\hline', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\vspace{0.1 mm}', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\\', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        """Generate Excel statistics"""
        try:
            excel_dict['BAM File'] = f'{self.bam}'
            excel_dict['Reference'] = f'{self.reference_name}'
            excel_dict['Reference Length'] = f'{self.reference_length:,}'
            excel_dict['Genome with Coverage'] = f'{(self.genome_coverage*100):,.1f}%'
            excel_dict['Average Depth'] = f'{self.ave_coverage:,.1f}'
            excel_dict['No Coverage Bases'] = f'{self.total_zero_coverage:,}'
            excel_dict['Percent Ref with Zero Coverage'] = f'{self.percent_ref_with_zero_coverage:,.6f}%'
            excel_dict['Ambiguous SNPs'] = f'{self.ac1_count:,}'
            excel_dict['Quality SNPs'] = f'{self.good_snp_count:,}'
        except Exception as e:
            print(f"Error generating Excel statistics: {e}")

if __name__ == "__main__":  # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Add zero coverage positions to VCF file

    Usage:
    vsnp3_zero_coverage.py -r *fasta -b *_nodup.bam -c *vcf

    If multiple FASTAs were used the concatenated FASTA to build alignment is needed

    '''), epilog='''---------------------------------------------------------''')
    parser.add_argument('-r', '--reference', action='store', dest='FASTA', required=True, default=None, help="Reference used to build alignment")
    parser.add_argument('-b', '--bam', action='store', dest='bam', required=False, default=None, help='bam file used to make VCF')
    parser.add_argument('-c', '--vcf', action='store', dest='vcf', required=False, default=None, help='VCF file')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    # Validate required arguments
    if not args.bam:
        print("Error: BAM file is required (-b/--bam)")
        sys.exit(1)
    if not args.vcf:
        print("Error: VCF file is required (-c/--vcf)")
        sys.exit(1)
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    try:
        zero_coverage = Zero_Coverage(FASTA=args.FASTA, bam=args.bam, vcf=args.vcf, debug=args.debug)

        # LaTeX report
        latex_report = Latex_Report(zero_coverage.sample_name)
        zero_coverage.latex(latex_report.tex)
        latex_report.latex_ending()

        # Excel Stats
        excel_stats = Excel_Stats(zero_coverage.sample_name)
        zero_coverage.excel(excel_stats.excel_dict)
        excel_stats.post_excel()
        
    except Exception as e:
        print(f"Error running zero coverage analysis: {e}")
        sys.exit(1)

    # temp_dir = './temp'
    # if not os.path.exists(temp_dir):
    #     os.makedirs(temp_dir)
    # files_grab = []
    # for files in ('*.aux', '*.log', '*tex', '*png', '*out', '*_all.bam', '*.bai', '*_filtered_hapall.vcf', '*_mapfix_hapall.vcf', '*_unfiltered_hapall.vcf', '*.sam', '*.amb', '*.ann', '*.bwt', '*.pac', '*.fasta.sa', '*_sorted.bam', '*.dict', 'chrom_ranges.txt', 'dup_metrics.csv', '*.fai'):
    #     files_grab.extend(glob.glob(files))
    # for each in files_grab:
    #     shutil.move(each, temp_dir)

    # if args.debug is False:
    #     shutil.rmtree(temp_dir)

# Created March 2021 by Tod Stuber