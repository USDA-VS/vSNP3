#!/usr/bin/env python3

__version__ = "3.33"

import gzip
import os
import argparse
import random
import sys
import logging

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
except ImportError as e:
    print("Error: Required BioPython module not found. Please install with: pip install biopython")
    sys.exit(1)

# Ambiguity codes
ambiguity_codes = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'T', 'C', 'G']
}

class Fasta_to_Paired_Fastq:

    def __init__(self, fasta_file, coverage, read_length):
        self.fasta_file = fasta_file
        self.coverage = coverage
        self.read_length = read_length
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        
        # Validate inputs
        if not os.path.exists(fasta_file):
            raise FileNotFoundError("Input FASTA file not found: {}".format(fasta_file))
        if coverage <= 0:
            raise ValueError("Coverage must be positive, got: {}".format(coverage))
        if read_length <= 0:
            raise ValueError("Read length must be positive, got: {}".format(read_length))
            
        self.process_fasta(fasta_file, coverage, read_length)

    def fake_quality_scores(self, length):
        """Generate a fake quality score string with broader range (Phred 20-40)."""
        # Use ASCII 53-73 for Phred scores 20-40 (!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI)
        return ''.join(chr(random.randint(53, 73)) for _ in range(length))

    def calculate_total_genome_length(self, fasta_file):
        """Calculate the total genome length from the input FASTA file."""
        total_length = 0
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                total_length += len(record.seq)
        except Exception as e:
            raise RuntimeError("Error reading FASTA file {}: {}".format(fasta_file, str(e)))
        
        if total_length == 0:
            raise ValueError("FASTA file appears to be empty or contains no valid sequences")
        return total_length

    def replace_ambiguities(self, sequence):
        """Replace ambiguity codes in the sequence with random nucleotides."""
        return ''.join(random.choice(ambiguity_codes[base]) if base in ambiguity_codes else base for base in sequence)

    def generate_paired_reads(self, sequence, num_reads, read_length):
        """Generate paired-end reads from the given sequence."""
        seq_len = len(sequence)
        reads = []
        
        # Validate sequence length
        if seq_len < read_length:
            logging.warning("Sequence length ({}) is shorter than read length ({}). Padding will be applied.".format(seq_len, read_length))
        
        # Calculate insert size (distance between read starts)
        insert_size = min(read_length + 50, max(seq_len // 2, read_length))  # Ensure minimum viable insert size
        
        for i in range(num_reads):
            # Ensure we don't go beyond sequence boundaries
            max_start = max(0, seq_len - insert_size)
            
            if max_start <= 0:
                # If sequence is too short, just take reads from the beginning
                start_pos = 0
            else:
                start_pos = random.randint(0, max_start)
            
            # R1: forward read from start_pos
            end_pos_r1 = min(start_pos + read_length, seq_len)
            seq_r1 = self.replace_ambiguities(sequence[start_pos:end_pos_r1].upper())
            
            # R2: reverse read from end of insert, reverse complemented
            r2_start = max(0, start_pos + insert_size - read_length)
            r2_end = min(r2_start + read_length, seq_len)
            seq_r2 = self.replace_ambiguities(sequence[r2_start:r2_end].upper())
            
            # Pad reads if they are shorter than read_length
            if len(seq_r1) < read_length:
                seq_r1 = seq_r1 + 'N' * (read_length - len(seq_r1))
            if len(seq_r2) < read_length:
                seq_r2 = seq_r2 + 'N' * (read_length - len(seq_r2))

            # Reverse complement R2 (standard paired-end orientation)
            try:
                seq_r2 = str(Seq(seq_r2).reverse_complement())
            except Exception as e:
                logging.error("Error creating reverse complement for read {}: {}".format(i+1, str(e)))
                # Fallback: simple reverse complement
                complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
                seq_r2 = ''.join(complement.get(base, base) for base in seq_r2[::-1])
            
            # Calculate coordinates safely
            x_coord = start_pos + 1
            y_coord = start_pos + insert_size
            
            reads.append((seq_r1, seq_r2, x_coord, y_coord, i + 1))

        return reads

    def process_fasta(self, fasta_file, coverage, read_length):
        """Process the input FASTA file and generate paired-end FASTQ files."""
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
        total_length = self.calculate_total_genome_length(fasta_file)
        
        # Use float division to ensure precision, then convert to int
        num_reads = int(float(total_length * coverage) / float(read_length))
        
        if num_reads <= 0:
            raise ValueError("Calculated number of reads is zero or negative. Check coverage and read length parameters.")
        
        logging.info("Total genome length: {} bp".format(total_length))
        logging.info("Target coverage: {}X".format(coverage))
        logging.info("Generating {} reads of length {} bp".format(num_reads, read_length))
        
        fastq_r1_file = "{}_R1.fastq.gz".format(file_prefix)
        fastq_r2_file = "{}_R2.fastq.gz".format(file_prefix)
        
        # Fixed barcode without invalid characters
        barcode = "AAACCCGGGTTT"

        try:
            with gzip.open(fastq_r1_file, 'wt') as r1, gzip.open(fastq_r2_file, 'wt') as r2:
                read_counter = 1
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequence = str(record.seq)
                    if len(sequence) == 0:
                        logging.warning("Skipping empty sequence record: {}".format(record.id))
                        continue
                        
                    reads = self.generate_paired_reads(sequence, num_reads, read_length)
                    
                    for seq_r1, seq_r2, x_coord, y_coord, unique_id in reads:
                        # Standard Illumina-like headers without problematic characters
                        # Using string formatting instead of f-strings for broader Python compatibility
                        header_r1 = "@SIM:{}:FC:1:1:{}:{} 1:N:0:{}".format(read_counter, x_coord, y_coord, barcode)
                        header_r2 = "@SIM:{}:FC:1:1:{}:{} 2:N:0:{}".format(read_counter, x_coord, y_coord, barcode)
                        
                        # Validate sequences contain only valid nucleotides for final output
                        valid_bases = set('ATCGN')
                        seq_r1_clean = ''.join(base if base in valid_bases else 'N' for base in seq_r1)
                        seq_r2_clean = ''.join(base if base in valid_bases else 'N' for base in seq_r2)
                        
                        record_r1 = "{}\n{}\n+\n{}\n".format(header_r1, seq_r1_clean, self.fake_quality_scores(len(seq_r1_clean)))
                        record_r2 = "{}\n{}\n+\n{}\n".format(header_r2, seq_r2_clean, self.fake_quality_scores(len(seq_r2_clean)))
                        
                        r1.write(record_r1)
                        r2.write(record_r2)
                        read_counter += 1
        
        except IOError as e:
            raise RuntimeError("Error writing FASTQ files: {}".format(str(e)))
        except Exception as e:
            raise RuntimeError("Unexpected error during FASTQ generation: {}".format(str(e)))

        self.fastq_r1_file = fastq_r1_file
        self.fastq_r2_file = fastq_r2_file
        logging.info("Successfully created paired-end FASTQ files")

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Convert a FASTA file into paired-end FASTQ files with fake quality scores and specific coverage. If FASTA is reporting ambiguity codes, they will be replaced with representative mix of nucleotides.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-c", "--coverage", type=int, default=100, help="Desired coverage (default: 100X)")
    parser.add_argument("-l", "--read_length", type=int, default=150, help="Read length (default: 150)")
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.basename(__file__), __version__))
    
    args = parser.parse_args()

    try:
        fasta_to_paired_fastq = Fasta_to_Paired_Fastq(args.input, args.coverage, args.read_length)
        print("Paired-end FASTQ files created: \n\t{}\n\t{}".format(
            fasta_to_paired_fastq.fastq_r1_file, 
            fasta_to_paired_fastq.fastq_r2_file
        ))
    except Exception as e:
        logging.error("Script failed: {}".format(str(e)))
        sys.exit(1)