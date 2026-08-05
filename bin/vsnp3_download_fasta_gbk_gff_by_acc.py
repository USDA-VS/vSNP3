#!/usr/bin/env python

from vsnp3_version import __version__

import os
import argparse
import textwrap
from Bio import SeqIO
from Bio import Entrez


class Downloader:

    def __init__(self, accession):
        self.entrezDbName = 'nucleotide'
        # NCBI requires a real contact address so they can reach you before
        # blocking an address that is overloading E-utilities.  This was
        # 'mickey_mouse@gmail.com', which violates their usage policy and puts
        # every user sharing an institutional IP at risk of being throttled.
        self.email = os.environ.get('NCBI_EMAIL', '').strip()
        if not self.email:
            raise SystemExit(
                'NCBI requires a contact email address for E-utilities requests.\n'
                '  export NCBI_EMAIL="you@example.org"\n'
                'Optionally also set NCBI_API_KEY for a higher rate limit:\n'
                '  https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/')
        Entrez.email = self.email
        Entrez.tool = 'vsnp3'
        api_key = os.environ.get('NCBI_API_KEY', '').strip()
        if api_key:
            Entrez.api_key = api_key
        self.accession = accession

    def gbk(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gbk")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gbwithparts')
        writeFile = self.accession + ".gbk"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

    def gff(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} gff3")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='gff3')
        writeFile = self.accession + ".gff"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

    def fasta(self):
        Entrez.email = self.email
        print(f"Downloading {self.accession} FASTA")
        entryData = Entrez.efetch(db=self.entrezDbName, id=self.accession, retmode="text", rettype='fasta')
        writeFile = self.accession + ".fasta"
        local_file=open(writeFile,"w")
        local_file.write(entryData.read())
        entryData.close()
        local_file.close()

        handle = open(writeFile, "r")
        for record in SeqIO.parse(handle, "fasta"):
            print(f"{record.description}")
            print(f'Sequence length: {len(record.seq):,}\n')
        handle.close()
        return record.description
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------

    Usage:
        fasta_gbk_gff_by_acc.py -a NC_002945 -f
        fasta_gbk_gff_by_acc.py -a NC_006932 -fg
        fasta_gbk_gff_by_acc.py -a NC_006933 -fbg
        fasta_gbk_gff_by_acc.py -a CP023243 -fbg
        fasta_gbk_gff_by_acc.py -a NZ_CP023243 -fbg
        **bad request: fasta_gbk_gff_by_acc.py -a NZ_AQME0000000 -fbg # must be complete chromosome

    Search genomes: https://www.ncbi.nlm.nih.gov/genome

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-a', '--accession', action='store', dest='accession', help='NCBI chromosome number')
    parser.add_argument('-f', '--fasta', action='store_true', dest='fasta', help='get FASTA file')
    parser.add_argument('-b', '--gbk', action='store_true', dest='gbk', help='get gbk file')
    parser.add_argument('-g', '--gff', action='store_true', dest='gff', help='get gff file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    download = Downloader(args.accession)
    if args.fasta:
        download.fasta()
    if args.gbk:
        download.gbk()
    if args.gff:
        download.gff()
    
