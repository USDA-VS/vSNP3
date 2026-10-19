#!/usr/bin/env python

from vsnp3_version import __version__
from vsnp3_version import QUAL_THRESHOLD, MQ_THRESHOLD, AC_HOMOZYGOUS, AC_HETEROZYGOUS
from vsnp3_version import mixed_signature, MIXED_AC1_RATIO

import os
import io
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


def step2_snp_counts(vcf_df):
    '''
    (single-base AC=2 calls, those of them step 2 takes as SNP positions).

    The step 2 test at its default thresholds, as vsnp3_step2.py reads a
    record and Group._passing selects it: single-base REF and ALT, AC=2, QUAL
    above QUAL_THRESHOLD and MQ at least MQ_THRESHOLD, with QUAL and MQ
    truncated to integers and freebayes' MQM read as MQ.  A record missing
    any of these values does not pass.  tests/test_step2_inputs.py checks
    this agrees with Group._passing.
    '''
    if vcf_df.empty:
        return 0, 0
    info = vcf_df['INFO'].astype(str).str.replace(';MQM=', ';MQ=', regex=False)
    ac = pd.to_numeric(info.str.extract(r'(?:^|;)AC=([^;]*)', expand=False), errors='coerce')
    mq = np.trunc(pd.to_numeric(info.str.extract(r'(?:^|;)MQ=([^;]*)', expand=False),
                                errors='coerce'))
    qual = np.trunc(pd.to_numeric(vcf_df['QUAL'], errors='coerce'))
    snp = ((vcf_df['REF'].astype(str).str.len() == 1) & (vcf_df['ALT'].astype(str).str.len() == 1)
           & (ac == AC_HOMOZYGOUS))
    passing = snp & (qual > QUAL_THRESHOLD) & (mq >= MQ_THRESHOLD)
    return int(snp.sum()), int(passing.sum())


def confident_mixed_count(vcf_df):
    '''
    Single-base AC=1 calls with QUAL >= QUAL_THRESHOLD and MQ >= MQ_THRESHOLD:
    the confident mixed calls the mixed-sample signature counts
    (vsnp3_version.mixed_signature).  Read the same way as step2_snp_counts.
    '''
    if vcf_df.empty:
        return 0
    info = vcf_df['INFO'].astype(str).str.replace(';MQM=', ';MQ=', regex=False)
    ac = pd.to_numeric(info.str.extract(r'(?:^|;)AC=([^;]*)', expand=False), errors='coerce')
    mq = np.trunc(pd.to_numeric(info.str.extract(r'(?:^|;)MQ=([^;]*)', expand=False),
                                errors='coerce'))
    qual = np.trunc(pd.to_numeric(vcf_df['QUAL'], errors='coerce'))
    mixed = ((vcf_df['REF'].astype(str).str.len() == 1) & (vcf_df['ALT'].astype(str).str.len() == 1)
             & (ac == AC_HETEROZYGOUS) & (qual >= QUAL_THRESHOLD) & (mq >= MQ_THRESHOLD))
    return int(mixed.sum())


def mixed_sample_notice(step2_snp_count, confident_ac1):
    '''
    The step 1 report notice for a sample that has SNPs passing the step 2
    thresholds but whose confident mixed calls outnumber them by more than
    MIXED_AC1_RATIO to 1: the signature of a mixture of strains.  None when the signature is absent or the counts
    are unknown.  A sample with no passing SNP gets no_snp_notice instead.
    '''
    if not step2_snp_count or confident_ac1 is None:
        return None
    if not mixed_signature(step2_snp_count, confident_ac1):
        return None
    return (f'Possible mixed sample: {confident_ac1:,} confident mixed call(s) (AC=1, '
            f'QUAL >= {QUAL_THRESHOLD}, MQ >= {MQ_THRESHOLD}) against {step2_snp_count:,} '
            f'SNP(s) passing the step 2 thresholds. A single strain rarely has more than '
            f'{MIXED_AC1_RATIO:g} mixed calls per SNP; two strains together do. The sample is '
            f'highlighted in the Groupings table of the step 2 summary report.')


def sample_name_for(bam, sample_name=None):
    '''
    The sample name the _zc.vcf is written under: the caller's when it gives
    one, otherwise the BAM file name up to its first '_' or '.'.

    The name used to be re-derived from the BAM argument every time, so
    -n/--SAMPLE_NAME with an underscore or a dot in it (B26-0290_S3) named the
    BAM, the report and the statistics B26-0290_S3 but the VCF B26-0290, and two
    forced names sharing a prefix collided in step 2.  From the command line
    the whole path was used, so vsnp3_zero_coverage.py -b out/x_nodup.bam wrote
    out/x_zc.vcf... except that re.sub left 'out/x', a path, as the name.
    '''
    if sample_name:
        return sample_name
    return re.sub('[_.].*', '', os.path.basename(bam))


VCF_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample"]


def read_vcf_lines(vcf):
    '''
    (header lines, record lines) of a VCF, header lines being those that start
    with '#', each with its line ending removed.  Non-UTF-8 bytes are replaced,
    not fatal: freebayes writes the reference path into the header.
    '''
    header, body = [], []
    with open(vcf, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('#'):
                header.append(line)
            elif line.strip():
                body.append(line)
    return header, body


def vcf_body_frame(body_lines):
    '''
    The records as a DataFrame with VCF_COLUMNS.  The header is skipped by line
    here, not with read_csv(comment='#'), which cut every record short at a '#'
    anywhere in the line (a sample or contig name, an INFO value).
    '''
    if not body_lines:
        return pd.DataFrame(columns=VCF_COLUMNS)
    return pd.read_csv(io.StringIO('\n'.join(body_lines) + '\n'), sep='\t',
                       header=None, names=VCF_COLUMNS)


def no_snp_notice(step2_snp_count, ac1_count, ac2_snp_count):
    '''
    The step 1 report notice for a sample none of whose calls step 2 will
    take as a SNP position, or None when some will (or the count is unknown).
    '''
    if step2_snp_count != 0:
        return None
    ac1 = ac1_count or 0
    ac2 = ac2_snp_count or 0
    rule = (f'the step 2 thresholds (single-base, AC=2, QUAL > {QUAL_THRESHOLD}, '
            f'MQ >= {MQ_THRESHOLD})')
    if ac1 > 0:
        return (f'No SNP passed {rule}, but {ac1:,} position(s) '
                f'were called mixed (AC=1). This usually means a mixed or '
                f'contaminated sample. In step 2 this sample carries no SNPs of '
                f'its own and is flagged in the summary report.')
    if ac2 > 0:
        return (f'No SNP passed {rule}: all {ac2:,} single-base AC=2 call(s) '
                f'have lower QUAL or MQ. In step 2 this sample carries no SNPs '
                f'of its own.')
    return (f'No SNP passed {rule} and no position was called '
            f'mixed: this sample matches the reference wherever it has coverage.')


class Zero_Coverage(Setup):
    ''' 
    Class to identify zero coverage regions in BAM alignments and 
    incorporate these into VCF files
    '''
    def __init__(self, FASTA=None, bam=None, vcf=None, debug=False, sample_name=None):

        def count_ac1_positions(vcf_file):
            """Count positions with AC=1 in VCF file"""
            ac1_count = 0

            try:
                with open(vcf_file, 'r', encoding='utf-8', errors='replace') as f:
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
        self.sample_name = sample_name_for(bam, sample_name)
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
            with open(FASTA, encoding='utf-8', errors='replace') as fasta_handle:
                contig_lengths = {record.id: len(record.seq)
                                  for record in SeqIO.parse(fasta_handle, "fasta")}
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
        
        # One read of the VCF: the header lines are kept for the _zc.vcf, the
        # records go into the frame.
        vcf_header_lines = []
        try:
            vcf_header_lines, vcf_body_lines = read_vcf_lines(vcf)
            vcf_df = vcf_body_frame(vcf_body_lines)
        except pd.errors.EmptyDataError:
            # Handle empty VCF file
            vcf_df = pd.DataFrame(columns=VCF_COLUMNS)
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

        # "Quality SNPs" above counts any AC, so a mixed sample whose AC=1 calls
        # are of high quality has plenty of them and none that step 2 will use.
        # The step 1 no-SNP notice is decided on the step 2 definition instead.
        try:
            ac2_snp_count, step2_count = step2_snp_counts(vcf_df)
            print(f"\tSNPs passing the step 2 thresholds: {step2_count:,} "
                  f"(of {ac2_snp_count:,} single-base AC=2 calls)")
        except Exception as e:                                   # noqa: BLE001
            print(f"Error counting SNPs by the step 2 thresholds: {e}")
            ac2_snp_count = step2_count = None
        try:
            confident_ac1 = confident_mixed_count(vcf_df)
            print(f"\tConfident mixed calls (AC=1, QUAL >= {QUAL_THRESHOLD}, "
                  f"MQ >= {MQ_THRESHOLD}): {confident_ac1:,}")
        except Exception as e:                                   # noqa: BLE001
            print(f"Error counting confident mixed calls: {e}")
            confident_ac1 = None

        # Handle division by zero
        if reference_length > 0:
            percent_ref_with_good_snp_count = good_snp_count / reference_length * 100
        else:
            percent_ref_with_good_snp_count = 0
            
        zero_coverage_vcf = f'{self.sample_name}_zc.vcf'
        # Named for the sample: two samples run side by side in one directory
        # used to share v_header.csv and v_annotated_body.csv, so a _zc.vcf
        # could carry the other sample's header.
        header_csv = f'{self.sample_name}_v_header.csv'
        body_csv = f'{self.sample_name}_v_annotated_body.csv'

        if total_zero_coverage > 0:
            try:
                # Write header
                with open(header_csv, 'w', encoding='utf-8', errors='replace') as header_out:
                    for line in vcf_header_lines:
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
                    combined_df.to_csv(body_csv, sep='\t', header=False, index=False)
                else:
                    # If no zero positions, just write the original VCF data
                    vcf_df.to_csv(body_csv, sep='\t', header=False, index=False)

                # Use context manager for safer file operations
                with open(zero_coverage_vcf, "wb") as outfile:
                    for cf in [header_csv, body_csv]:
                        if os.path.exists(cf):
                            with open(cf, "rb") as infile:
                                outfile.write(infile.read())

                # Clean up temporary files
                for temp_file in [header_csv, body_csv]:
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
                for temp_file in [header_csv, body_csv]:
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
        self.ac2_snp_count = ac2_snp_count
        self.step2_snp_count = step2_count
        self.confident_ac1_count = confident_ac1
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
    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, default=None, help='Optional: name the _zc.vcf for this sample; default is the BAM file name up to its first "_" or "."')
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
        zero_coverage = Zero_Coverage(FASTA=args.FASTA, bam=args.bam, vcf=args.vcf, debug=args.debug,
                                      sample_name=args.SAMPLE_NAME)

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