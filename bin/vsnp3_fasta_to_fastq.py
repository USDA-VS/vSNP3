#!/usr/bin/env python3

__version__ = "3.22"

import gzip
import os
import argparse
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import seq1
import logging
from datetime import datetime

# Ambiguity codes
ambiguity_codes = {
    'R': ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'], 
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'T', 'C', 'G']
}
class Fasta_to_Paired_Fastq():

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


    def fake_quality_scores(self, length):
        """Generate a fake quality score string of given length with varying quality scores."""
        return ''.join(chr(random.randint(58, 72)) for _ in range(length))

    def calculate_total_genome_length(self, fasta_file):
        """Calculate the total genome length from the input FASTA file."""
        total_length = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_length += len(record.seq)
        return total_length

    def generate_paired_reads(self, sequence, num_reads, read_length):

        """Generate paired-end reads from the given sequence."""
        seq_len = len(sequence)
        reads = []
        for _ in range(num_reads):
            start_pos = random.randint(0, seq_len - read_length)
            end_pos = start_pos + read_length

            # Create read pair
            seq_r1 = sequence[start_pos:start_pos + read_length // 2]
            seq_r2 = sequence[end_pos - read_length // 2:end_pos]

            # Pad reads if they are shorter than the read length
            if len(seq_r1) < read_length // 2:
                seq_r1 = seq_r1 + 'N' * (read_length // 2 - len(seq_r1))
            if len(seq_r2) < read_length // 2:
                seq_r2 = seq_r2 + 'N' * (read_length // 2 - len(seq_r2))

            reads.append((seq_r1, seq_r2))
        
        return reads

    def is_valid_sequence(self, sequence):
        global ambiguity_codes
        """Check if the sequence contains only valid nucleotide characters."""
        valid_nucleotides = set('ATCG') | set(ambiguity_codes.keys())
        invalid_chars = set(sequence.upper()) - valid_nucleotides
        if invalid_chars:
            logging.warning(f"Invalid characters found in sequence: {invalid_chars}. Treating them as 'N'.")
            return False
        return True

    def handle_ambiguity(self, sequence):
        global ambiguity_codes
        """Handle ambiguity in sequence by randomly choosing one of the possible nucleotides."""
        new_sequence = []
        for base in sequence:
            if base in ambiguity_codes:
                new_sequence.append(random.choice(ambiguity_codes[base]))
            else:
                new_sequence.append(base)
        return ''.join(new_sequence)

    def reverse_complement(self, sequence):
        """Generate the reverse complement of a sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement[base] for base in reversed(sequence))

    def generate_barcode(self, length=6):
        """Generate a random barcode sequence of a given length."""
        return ''.join(random.choice('ATCG') for _ in range(length))

    def __init__(self, fasta_file, coverage, read_length):
        """Convert a FASTA file into paired-end FASTQ files with fake quality scores and specific coverage."""
        print(f'Converting FASTA to FASTQ')
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

        # Determine output file names based on input file name
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        fastq_r1_file = base_name + "_from_FASTA_R1.fastq.gz"
        fastq_r2_file = base_name + "_from_FASTA_R2.fastq.gz"

        # Calculate total genome length
        total_genome_length = self.calculate_total_genome_length(fasta_file)

        # Calculate the number of reads needed for the desired coverage
        num_reads = (coverage * total_genome_length) // read_length

        # Open output FASTQ files for writing
        with gzip.open(fastq_r1_file, "wt") as r1, gzip.open(fastq_r2_file, "wt") as r2:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = str(record.seq).upper()
                if not self.is_valid_sequence(sequence):
                    # Replace invalid characters with 'N'
                    sequence = ''.join([base if base in {'A', 'T', 'C', 'G', 'N'} else 'N' for base in sequence])

                reads = self.generate_paired_reads(sequence, num_reads, read_length)

                for i, (seq_r1, seq_r2) in enumerate(reads):
                    # Handle ambiguities in each read pair
                    seq_r1 = self.handle_ambiguity(seq_r1)
                    seq_r2 = self.handle_ambiguity(seq_r2)
                    seq_r2 = self.reverse_complement(seq_r2)  # Reverse complement for R2

                    # Create varying quality scores between Phred 20 and 35
                    qual_r1 = self.fake_quality_scores(len(seq_r1))
                    qual_r2 = self.fake_quality_scores(len(seq_r2))

                    # Generate a random barcode sequence
                    barcode = self.generate_barcode()

                    # Create Illumina-like headers with barcode
                    header_base = f"{base_name}:{i+1}:N:0:1:{barcode}"
                    record_r1 = f"@{header_base}/1\n{seq_r1}\n+\n{qual_r1}\n"
                    record_r2 = f"@{header_base}/2\n{seq_r2}\n+\n{qual_r2}\n"

                    # Write to FASTQ files
                    r1.write(record_r1)
                    r2.write(record_r2)

                    self.fastq_r1_file = fastq_r1_file
                    self.fastq_r2_file = fastq_r2_file

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Convert a FASTA file into paired-end FASTQ files with fake quality scores and specific coverage.  If FASTA is reporting ambiguity codes, they will be replaced with resepentative mix of nucleotides.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-c", "--coverage", type=int, default=100, help="Desired coverage (default: 100X)")
    parser.add_argument("-l", "--read_length", type=int, default=300, help="Read length (default: 300)")
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    fasta_to_paired_fastq = Fasta_to_Paired_Fastq(args.input, args.coverage, args.read_length)
    print(f"Paired-end FASTQ files created: \n\t{fasta_to_paired_fastq.fastq_r1_file}\n\t{fasta_to_paired_fastq.fastq_r2_file}")