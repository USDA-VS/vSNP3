#!/usr/bin/env python3

__version__ = "3.26"

import gzip
import os
import argparse
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
from datetime import datetime

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
        self.process_fasta(fasta_file, coverage, read_length)

    def fake_quality_scores(self, length):
        """Generate a fake quality score string of given length with varying quality scores."""
        return ''.join(chr(random.randint(58, 69)) for _ in range(length))

    def calculate_total_genome_length(self, fasta_file):
        """Calculate the total genome length from the input FASTA file."""
        total_length = 0
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_length += len(record.seq)
        return total_length

    def replace_ambiguities(self, sequence):
        """Replace ambiguity codes in the sequence with random nucleotides."""
        return ''.join(random.choice(ambiguity_codes[base]) if base in ambiguity_codes else base for base in sequence)

    def generate_paired_reads(self, sequence, num_reads, read_length):
        """Generate paired-end reads from the given sequence."""
        seq_len = len(sequence)
        reads = []
        for i in range(num_reads):
            adjusted_read_length = min(read_length, seq_len)
            start_pos = random.randint(0, seq_len - adjusted_read_length)
            end_pos = start_pos + adjusted_read_length

            # Create read pair
            seq_r1 = self.replace_ambiguities(sequence[start_pos:start_pos + adjusted_read_length // 2].upper())
            seq_r2 = self.replace_ambiguities(sequence[end_pos - adjusted_read_length // 2:end_pos].upper())

            # Pad reads if they are shorter than the adjusted read length
            if len(seq_r1) < adjusted_read_length // 2:
                seq_r1 = seq_r1 + 'N' * (adjusted_read_length // 2 - len(seq_r1))
            if len(seq_r2) < adjusted_read_length // 2:
                seq_r2 = seq_r2 + 'N' * (adjusted_read_length // 2 - len(seq_r2))

            seq_r2 = str(Seq(seq_r2).reverse_complement())  # Reverse complement for R2
            x_coord = start_pos + 1
            y_coord = end_pos
            reads.append((seq_r1, seq_r2, x_coord, y_coord, i + 1))

        return reads

    def process_fasta(self, fasta_file, coverage, read_length):
        """Process the input FASTA file and generate paired-end FASTQ files."""
        file_prefix = os.path.basename(fasta_file).split('.')[0]
        total_length = self.calculate_total_genome_length(fasta_file)
        num_reads = int((total_length * coverage) / read_length)
        
        fastq_r1_file = f"{file_prefix}_R1.fastq.gz"
        fastq_r2_file = f"{file_prefix}_R2.fastq.gz"
        constant_overlap_seq = "TTCAAGTATG+CGATACCATC"

        with gzip.open(fastq_r1_file, 'wt') as r1, gzip.open(fastq_r2_file, 'wt') as r2:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = str(record.seq)
                reads = self.generate_paired_reads(sequence, num_reads, read_length)
                for seq_r1, seq_r2, x_coord, y_coord, unique_id in reads:
                    header_r1 = f"@{file_prefix}:8:FASTQFROMFASTA:1:{x_coord}:{y_coord}:{unique_id} 1:N:0:{constant_overlap_seq}"
                    header_r2 = f"@{file_prefix}:8:FASTQFROMFASTA:1:{x_coord}:{y_coord}:{unique_id} 2:N:0:{constant_overlap_seq}"
                    
                    record_r1 = f"{header_r1}\n{seq_r1}\n+\n{self.fake_quality_scores(len(seq_r1))}\n"
                    record_r2 = f"{header_r2}\n{seq_r2}\n+\n{self.fake_quality_scores(len(seq_r2))}\n"
                    
                    r1.write(record_r1)
                    r2.write(record_r2)

        self.fastq_r1_file = fastq_r1_file
        self.fastq_r2_file = fastq_r2_file

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description="Convert a FASTA file into paired-end FASTQ files with fake quality scores and specific coverage. If FASTA is reporting ambiguity codes, they will be replaced with representative mix of nucleotides.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-c", "--coverage", type=int, default=100, help="Desired coverage (default: 100X)")
    parser.add_argument("-l", "--read_length", type=int, default=300, help="Read length (default: 300)")
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    fasta_to_paired_fastq = Fasta_to_Paired_Fastq(args.input, args.coverage, args.read_length)
    print(f"Paired-end FASTQ files created: \n\t{fasta_to_paired_fastq.fastq_r1_file}\n\t{fasta_to_paired_fastq.fastq_r2_file}")
