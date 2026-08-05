#!/usr/bin/env python

from vsnp3_version import __version__

import os
import re
import shutil
import locale
import argparse
import textwrap
import numpy as np
import pandas as pd
from Bio import SeqIO
import subprocess
import sys

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
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
        # Per-contig depth as numpy arrays, streamed.
        #
        # This replaced a version that held four whole-genome copies at once: the
        # ~99 MB `samtools depth` stdout as one string, its 4.35 M-element
        # splitlines() list, a {chrom-pos: depth} dict of strings, and a second
        # {chrom-pos: 0} dict for every reference base -- then outer-merged two
        # DataFrames on a 4.35 M string index.  Measured at 7.5 s and 1.93 GB peak
        # for M. bovis.  At ten concurrent samples on one node that is ~19 GB, and
        # swapping turns a modest CPU cost into unbounded wall-clock, which is the
        # failure people actually hit on a laptop or under WSL.
        #
        # int32 per base is ~17 MB for the same genome.  `samtools depth -a` emits
        # every position rather than only covered ones, so the all-zeros dict and
        # the merge are unnecessary.
        try:
            contig_lengths = {record.id: len(record.seq)
                              for record in SeqIO.parse(FASTA, "fasta")}
        except Exception as e:
            print(f"Error parsing FASTA file {FASTA}: {e}")
            raise
        reference_length = sum(contig_lengths.values())

        # Sized from the FASTA, not from what samtools reports: `depth -a` covers
        # only contigs present in the BAM header, so a reference contig with no
        # reads at all must start as zeros here to be counted as uncovered.
        depth_arrays = {chrom: np.zeros(length, dtype=np.int32)
                        for chrom, length in contig_lengths.items()}

        unknown_contigs = set()
        try:
            proc = subprocess.Popen(['samtools', 'depth', '-a', bam],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    text=True)
            for line in proc.stdout:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 3:
                    continue
                array = depth_arrays.get(parts[0])
                if array is None:
                    unknown_contigs.add(parts[0])
                    continue
                index = int(parts[1]) - 1
                if 0 <= index < array.size:
                    array[index] = int(parts[2])
            proc.stdout.close()
            stderr = proc.stderr.read()
            proc.stderr.close()
            if proc.wait() != 0:
                raise subprocess.CalledProcessError(proc.returncode,
                                                    ['samtools', 'depth', '-a', bam],
                                                    stderr=stderr)
        except subprocess.CalledProcessError as e:
            print(f"Error running samtools depth: {e}")
            raise
        except FileNotFoundError:
            print("Error: samtools not found. Please ensure samtools is installed and in your PATH.")
            raise

        if unknown_contigs:
            # A BAM aligned against a different reference than the FASTA given here.
            # Previously these positions were silently folded into the denominator by
            # the outer merge, quietly changing every coverage percentage.
            print(f'Warning: {len(unknown_contigs)} contig(s) in {os.path.basename(bam)} '
                  f'are not in {os.path.basename(FASTA)} and were ignored: '
                  f'{", ".join(sorted(unknown_contigs)[:5])}')

        total_length = reference_length
        # dtype=np.int64 on the sum: int32 would overflow above ~2.1e9 total depth,
        # which a high-coverage bacterial genome can reach.
        total_depth = sum(int(a.sum(dtype=np.int64)) for a in depth_arrays.values())
        ave_coverage = (total_depth / total_length) if total_length else 0.0

        # (chrom, position) for every uncovered base, 1-based.  Replaces indexing a
        # DataFrame of every reference position and filtering it.
        zero_positions = [(chrom, int(offset) + 1)
                          for chrom, array in depth_arrays.items()
                          for offset in np.flatnonzero(array == 0)]
        total_zero_coverage = len(zero_positions)

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
                # zero_positions already holds (chrom, position) pairs, so there is no
                # 'chrom-pos' string to split back apart.  The old form did
                # idx.rsplit('-', 1), which was one hyphen in a contig name away from
                # silently mis-parsing.
                zero_positions_list = [{
                    'CHROM': chrom,
                    'POS': position,
                    'ID': '.',
                    'REF': 'N',
                    'ALT': '.',
                    'QUAL': '.',
                    'FILTER': '.',
                    'INFO': '.',
                    'FORMAT': 'GT',
                    'Sample': './.'
                } for chrom, position in zero_positions]
                
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
                        
            except (OSError, KeyError, ValueError, TypeError,
                    pd.errors.EmptyDataError) as e:
                # This used to fall back to shutil.copyfile(vcf, zero_coverage_vcf).
                # _zc.vcf is the artifact that records the difference between "no
                # coverage here" and "matches the reference here"; an unmodified
                # copy asserts a reference call at every uncovered position, and
                # nothing downstream can tell.  With total_zero_coverage > 0 we know
                # those positions exist, so writing the file anyway is not an option.
                if os.path.exists(zero_coverage_vcf):
                    os.remove(zero_coverage_vcf)
                for temp_file in ['v_header.csv', 'v_annotated_body.csv']:
                    if os.path.exists(temp_file):
                        os.remove(temp_file)
                raise RuntimeError(
                    f'could not insert {total_zero_coverage} zero-coverage positions '
                    f'into {os.path.basename(vcf)}: {type(e).__name__}: {e}. '
                    f'No {zero_coverage_vcf} was written -- an unmodified copy would '
                    f'report those positions as matching the reference.') from e
        else:
            # No zero-coverage positions to insert, so the VCF is already correct.
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
    #     files_grab.extend(sorted(glob.glob(files)))
    # for each in files_grab:
    #     shutil.move(each, temp_dir)

    # if args.debug is False:
    #     shutil.rmtree(temp_dir)

# Created March 2021 by Tod Stuber