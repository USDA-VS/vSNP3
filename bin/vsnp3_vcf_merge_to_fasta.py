#!/usr/bin/env python

from vsnp3_version import __version__

import os
import re
import io
import sys
import argparse
import textwrap
import tempfile
import subprocess
import shutil
from contextlib import contextmanager

try:
    import numpy as np
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError as e:
    print("Error: Required dependencies not found. Please install:")
    print("pip install numpy pandas biopython")
    sys.exit(1)

# Handle collections import for different Python versions
try:
    from collections.abc import Mapping
except ImportError:
    from collections import Mapping

try:
    from vsnp3_file_setup import Setup
except ImportError:
    print("Warning: vsnp3_file_setup module not found. Using minimal Setup class.")
    class Setup:
        def __init__(self, **kwargs):
            pass


class Merge_VCF(Setup):
    """ 
    Merges VCF file SNPs into a FASTA reference sequence
    """

    def __init__(self, vcf=None, fasta=None, qual_threshold=300, mq_threshold=56, ac_threshold=1, indels=False):
        """
        Initialize the Merge_VCF class
        
        Parameters:
        -----------
        vcf : str
            Path to the VCF file
        fasta : str
            Path to the FASTA file
        qual_threshold : int
            Quality threshold for filtering SNPs
        mq_threshold : int
            Mapping quality threshold for filtering SNPs
        ac_threshold : int
            Allele count threshold for filtering SNPs
        indels : bool
            Whether to include indels in the output
        """
        try:
            Setup.__init__(self, FASTA=fasta)
        except:
            pass
            
        self.fasta = fasta
        self.vcf = vcf
        self.qual_threshold = int(qual_threshold) if qual_threshold is not None else 300
        self.mq_threshold = int(mq_threshold) if mq_threshold is not None else 56
        self.ac_threshold = int(ac_threshold) if ac_threshold is not None else 1
        self.indels = bool(indels)
        
        # Validate input files
        if not os.path.exists(self.fasta):
            raise FileNotFoundError("FASTA file not found: {}".format(self.fasta))
        if not os.path.exists(self.vcf):
            raise FileNotFoundError("VCF file not found: {}".format(self.vcf))
            
        # Generate safe names using string formatting instead of f-strings
        self.fasta_name = re.sub(r'[.].*', '', os.path.basename(self.fasta))
        vcf_basename = re.sub(r'[.].*', '', os.path.basename(self.vcf))
        self.vcf_name = vcf_basename.rstrip('_zc')
        self.output_name = "{}_vcf_merged_{}".format(self.vcf_name, self.fasta_name)

    @contextmanager
    def temporary_file(self, suffix='.vcf'):
        """Context manager for safe temporary file handling"""
        temp_fd, temp_path = tempfile.mkstemp(suffix=suffix)
        try:
            os.close(temp_fd)
            yield temp_path
        finally:
            try:
                os.unlink(temp_path)
            except OSError:
                pass

    def filter_vcf_safely(self, input_vcf, output_vcf):
        """
        Safely filter VCF file to remove entries with multiple REF or ALT values
        Uses subprocess instead of os.system for better security
        """
        try:
            # First awk command: filter out rows where column 4 (REF) contains comma
            cmd1 = ['awk', 'BEGIN { OFS=FS="\t" } $4 !~ /,/', input_vcf]
            
            # Second awk command: filter out rows where column 5 (ALT) contains comma
            cmd2 = ['awk', 'BEGIN { OFS=FS="\t" } $5 !~ /,/']
            
            # Execute pipeline
            with open(output_vcf, 'w') as output_file:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=output_file)
                proc1.stdout.close()
                proc2.wait()
                
            if proc2.returncode != 0:
                raise RuntimeError("VCF filtering failed")
                
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            # Fallback to Python-based filtering if awk is not available
            print("Warning: awk not available, using Python fallback")
            self.filter_vcf_python(input_vcf, output_vcf)
    
    def filter_vcf_python(self, input_vcf, output_vcf):
        """Python fallback for VCF filtering when awk is not available"""
        try:
            with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
                for line in infile:
                    if line.startswith('#'):
                        outfile.write(line)
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 5:
                        ref_field = fields[3]
                        alt_field = fields[4]
                        
                        # Skip lines with multiple REF or ALT values
                        if ',' not in ref_field and ',' not in alt_field:
                            outfile.write(line)
        except Exception as e:
            raise RuntimeError("Python VCF filtering failed: {}".format(str(e)))

    def read_vcf(self, path):
        """
        Read a VCF file and return a pandas DataFrame
        
        Parameters:
        -----------
        path : str
            Path to the VCF file
            
        Returns:
        --------
        pandas.DataFrame
            DataFrame containing the VCF data
        """
        try:
            with open(path, 'r') as f:
                lines = [l for l in f if not l.startswith('##')]
            
            if not lines:
                raise ValueError("No valid VCF data found")
            
            df = pd.read_csv(
                io.StringIO(''.join(lines)),
                dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
                    'QUAL': str, 'FILTER': str, 'INFO': str},
                sep='\t'
            ).rename(columns={'#CHROM': 'CHROM'})
            
            # Convert columns to appropriate types with error handling
            df['POS'] = pd.to_numeric(df['POS'], errors='coerce').fillna(0).astype(int)
            df['QUAL'] = pd.to_numeric(df['QUAL'], errors='coerce').fillna(0).astype(int)
            
            # Helper function to parse INFO field
            def parse_info(info_str):
                if pd.isna(info_str) or not isinstance(info_str, str):
                    return {}
                info_dict = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        info_dict[key] = value
                return info_dict
            
            # Extract INFO fields safely
            df['INFO_DICT'] = df['INFO'].apply(parse_info)
            df['AC'] = df['INFO_DICT'].apply(lambda x: x.get('AC', '0') if isinstance(x, dict) else '0')
            df['DP'] = df['INFO_DICT'].apply(lambda x: x.get('DP', '0') if isinstance(x, dict) else '0')
            df['MQ'] = df['INFO_DICT'].apply(lambda x: x.get('MQ', '0') if isinstance(x, dict) else '0')
            
            # Convert to numeric with error handling
            df['AC'] = pd.to_numeric(df['AC'], errors='coerce').fillna(0).astype(int)
            df['DP'] = pd.to_numeric(df['DP'], errors='coerce').fillna(0).astype(int)
            df['MQ'] = pd.to_numeric(df['MQ'], errors='coerce').fillna(0).astype(int)
            
            # Drop unnecessary columns safely
            cols_to_drop = ['INFO', 'INFO_DICT']
            optional_cols = ['ID', 'FILTER', 'FORMAT']
            
            for col in optional_cols:
                if col in df.columns:
                    cols_to_drop.append(col)
                    
            df = df.drop(columns=cols_to_drop)

            return df
            
        except Exception as e:
            raise RuntimeError("Failed to read VCF file: {}".format(str(e)))

    def run(self):
        """
        Run the VCF to FASTA merging process
        """
        try:
            print("Starting VCF to FASTA merge process...")
            
            # Use context manager for temporary file
            with self.temporary_file() as temp_vcf:
                print("Filtering VCF file...")
                self.filter_vcf_safely(self.vcf, temp_vcf)
                
                print("Reading and processing VCF data...")
                df = self.read_vcf(temp_vcf)
                
                if df.empty:
                    raise ValueError("No valid VCF entries found after filtering")
                
                print("Applying quality filters...")
                # Apply filters based on parameters
                initial_count = len(df)
                
                if self.indels:
                    df = df[(df['QUAL'] > self.qual_threshold) & 
                           (df['MQ'] >= self.mq_threshold) & 
                           (df['AC'] >= self.ac_threshold)]
                else:
                    df['REF'] = df['REF'].astype('str')
                    df['ALT'] = df['ALT'].astype('str')
                    df = df[(df['QUAL'] > self.qual_threshold) & 
                           (df['MQ'] >= self.mq_threshold) & 
                           (df['AC'] >= self.ac_threshold) & 
                           (df['REF'].str.len() == 1) & 
                           (df['ALT'].str.len() == 1)]
                
                filtered_count = len(df)
                print("Variants: {} initial, {} after filtering".format(initial_count, filtered_count))
                
                if df.empty:
                    print("Warning: No variants passed quality filters")
                    return
                
                # Keep only position and alt allele
                df = df[['POS', 'ALT']].copy()
                
                print("Reading reference sequence...")
                # Read reference sequence
                try:
                    sequence = SeqIO.read(self.fasta, "fasta")
                except Exception as e:
                    raise RuntimeError("Failed to read FASTA file: {}".format(str(e)))
                
                print("Processing {} variants...".format(len(df)))
                
                # Create dictionary of positions to alt alleles
                dfdict = dict(zip(df['POS'].tolist(), df['ALT'].tolist()))
                
                # Convert to 0-indexed positions
                vcf_dict_zero_index = {max(0, key-1): value for key, value in dfdict.items()}
                
                # Convert sequence to list and then to dictionary
                sequence_list = list(str(sequence.seq))
                seq_dict = dict(enumerate(sequence_list))
                
                # Merge dictionaries, with VCF variants overriding reference
                merge_dicts = seq_dict.copy()
                merge_dicts.update(vcf_dict_zero_index)
                
                # Convert back to sequence
                updated_seq = "".join(str(merge_dicts[i]) for i in sorted(merge_dicts.keys()))

                # Create new sequence record
                record = SeqRecord(
                    Seq(updated_seq),
                    id=self.output_name,
                    name=self.output_name,
                    description=""
                )

                # Write to file
                output_file = "{}.fasta".format(self.output_name)
                print("Writing output to: {}".format(output_file))
                
                with open(output_file, "w") as output_handle:
                    SeqIO.write(record, output_handle, "fasta")
                
                print("VCF to FASTA merge completed successfully!")
                print("Output file: {}".format(output_file))
                print("Final sequence length: {}".format(len(updated_seq)))
                
        except Exception as e:
            print("Error during VCF to FASTA merge: {}".format(str(e)))
            raise


def main():
    """Main function to handle command line execution"""
    parser = argparse.ArgumentParser(
        prog='PROG', 
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=textwrap.dedent('''\
        ---------------------------------------------------------
        Merge VCF SNPs into a FASTA reference sequence
        This script takes a VCF file and the FASTA file that was used
        to generate it, and creates a new FASTA file with the variants
        from the VCF incorporated into the reference sequence.
        '''), 
        epilog='''---------------------------------------------------------'''
    )

    parser.add_argument(
        '-v', '--vcf', 
        action='store', 
        dest='vcf', 
        required=True, 
        help='Required: VCF file will be merged into FASTA file it was made from'
    )
    parser.add_argument(
        '-f', '--fasta', 
        action='store', 
        dest='fasta', 
        required=True,
        help='Required: FASTA used to build the VCF file'
    )
    parser.add_argument(
        '-x', '--qual_threshold', 
        action='store', 
        dest='qual_threshold', 
        default=300, 
        type=int,
        required=False, 
        help='Optional: Minimum QUAL threshold for using a SNP (default: 300)'
    )
    parser.add_argument(
        '-y', '--mq_threshold', 
        action='store', 
        dest='mq_threshold', 
        default=56, 
        type=int,
        required=False, 
        help='Optional: Minimum MQ threshold for using a SNP (default: 56)'
    )
    parser.add_argument(
        '-z', '--ac', 
        action='store', 
        dest='ac_threshold', 
        default=1, 
        type=int,
        required=False, 
        help='Optional: Minimum AC threshold for using a SNP (default: 1)'
    )
    parser.add_argument(
        '-i', '--indels', 
        action='store_true', 
        dest='indels', 
        default=False, 
        help='Optional: Include indels'
    )
    parser.add_argument(
        '--version', 
        action='version', 
        version='{}: version {}'.format(os.path.basename(__file__), __version__)
    )
    
    try:
        args = parser.parse_args()
        
        print('\n{} SET ARGUMENTS:'.format(os.path.basename(__file__)))
        print(args)
        print("\n")

        merge_vcf = Merge_VCF(
            vcf=args.vcf, 
            fasta=args.fasta, 
            qual_threshold=args.qual_threshold, 
            mq_threshold=args.mq_threshold, 
            ac_threshold=args.ac_threshold, 
            indels=args.indels
        )
        merge_vcf.run()
        
    except Exception as e:
        print("Error: {}".format(str(e)))
        sys.exit(1)


if __name__ == "__main__":
    main()


# Created 2022 by Tod Stuber