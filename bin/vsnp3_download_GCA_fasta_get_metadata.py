#!/usr/bin/env python

__version__ = "3.30"

import os
import sys
import re
import shutil
import glob
import time
from pathlib import Path
from zipfile import ZipFile
import pandas as pd
import argparse
import textwrap
from datetime import datetime


class Downloader():
    ''' 
    Class to download genome data from NCBI using the datasets tool
    '''

    def __init__(self, gca_accession=None, annotation_file=False, rename=False, debug=False):
        '''
        Initialize the downloader with the given GCA accession number
        
        Parameters:
        -----------
        gca_accession : str
            The GCA accession number
        annotation_file : bool
            Whether to download annotation files (gff3, gbff)
        rename : bool
            Whether to rename files with metadata
        debug : bool
            Whether to enable debug mode
        '''
        cwd = os.getcwd()
        gca_accession = os.path.basename(gca_accession)
        self.gca_accession = gca_accession
        self.sample_name = gca_accession
     
        if annotation_file:
            file_download = f'datasets download genome accession {gca_accession} --filename {self.sample_name}.zip --include genome,gff3,gbff'
        else:
            file_download = f'datasets download genome accession {gca_accession} --filename {self.sample_name}.zip --include genome'
        os.system(file_download)
        count = 1
        while True and count < 5:
            try:
                with ZipFile(f'{self.sample_name}.zip', 'r') as zf:
                    zf.extractall()
            except Exception:  # More specific exception handling
                count += 1
                print(f'Download try {count} in progress')
                for secs in range(30, 0, -1):
                    time.sleep(1)
                    print('Pausing:', secs, end='\r')
                os.system(file_download)
                continue
            else:
                break

        # Use list comprehension instead of for loop for better readability
        jsonl_files = list(Path('.').rglob('*_report.jsonl'))
        if jsonl_files:
            jsonl = f'{cwd}/{jsonl_files[0].as_posix()}'
            df = pd.read_json(jsonl, typ='series')
            acc = df['accession']
            
            # Extract species information
            try:
                species = df['checkmInfo']['checkmMarkerSet'].split()[1]
                if species == 'tuberculosis':
                    species = "tbc"
            except (KeyError, IndexError):
                species = "species-not-listed"
                
            # Extract strain information
            try:
                strain = df['organism']['organismName']
            except KeyError:
                try:
                    strain = df['organism']['infraspecificNames']['strain']
                except KeyError:
                    try:
                        strain = df['organism']['infraspecificNames']['isolate']
                    except KeyError:
                        strain = "strain-not-listed"
                        
            # Clean up species and strain names
            species = re.sub(r'[ ?\./;:]', '_', species)
            strain = re.sub(r'[ ?\./;:]', '_', strain)
            
            print(f'{acc}\t{acc}_{species}_{strain}')
            metadata_string = f'{acc}_{species}_{strain}'
            print(f'metadata string: {metadata_string}')

            self.acc = acc
            self.species = species
            self.strain = strain
            self.metadata_string = metadata_string

            # Process FASTA file
            fasta_files = list(Path('.').rglob('*.fna'))
            if fasta_files:
                fasta = f'{cwd}/{fasta_files[0].as_posix()}'
                base_only = os.path.basename(fasta)
                base_name = re.sub('.fna', '', base_only)
                base_ext = re.sub('.fna', '.fasta', base_only)
                updated_file = f'{cwd}/{base_ext}'
                
                if rename:
                    os.rename(fasta, f'{metadata_string}.fasta')
                else:
                    os.rename(fasta, updated_file)
                    
                fasta = updated_file if not rename else f'{metadata_string}.fasta'
                print(f'fasta: {fasta}')
                self.fasta = fasta

                # Process annotation files if requested
                if annotation_file or rename:
                    gff_files = list(Path('.').rglob('*genomic.gff'))
                    if gff_files and annotation_file:
                        gff = f'{cwd}/{gff_files[0].as_posix()}'
                        if rename:
                            os.rename(gff, f'{cwd}/{metadata_string}.gff')
                            gff = f'{cwd}/{metadata_string}.gff'
                        else:
                            os.rename(gff, f'{cwd}/{base_name}.gff')
                            gff = f'{cwd}/{base_name}.gff'

                    gbk_files = list(Path('.').rglob('*genomic.gbff'))
                    if gbk_files and annotation_file:
                        gbk = f'{cwd}/{gbk_files[0].as_posix()}'
                        if rename:
                            os.rename(gbk, f'{cwd}/{metadata_string}.gbk')
                            gbk = f'{cwd}/{metadata_string}.gbk'
                        else:
                            os.rename(gbk, f'{cwd}/{base_name}.gbk')
                            gbk = f'{cwd}/{base_name}.gbk'

                        print(f'gff: {gff}')
                        print(f'gbk: {gbk}')
                        self.gff = gff
                        self.gbk = gbk

    def excel(self, excel_dict):
        """Update the Excel dictionary with sample metadata"""
        excel_dict['sample'] = f'{self.gca_accession}'
        excel_dict['metadata'] = f'{self.acc}_{self.species}_{self.strain}'


class Excel_Stats:
    """Class to handle Excel statistics"""

    def __init__(self, sample_name):
        """Initialize the Excel stats with the given sample name"""
        self.sample_name = sample_name
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.excel_filename = f'{sample_name}_{date_stamp}_stats.xlsx'
        excel_dict = {}
        excel_dict['sample'] = sample_name
        excel_dict['date'] = date_stamp
        self.excel_dict = excel_dict 

    def post_excel(self):
        """Create an Excel file with the collected statistics"""
        df = pd.DataFrame.from_dict(self.excel_dict, orient='index').T
        df = df.set_index('sample')
        df.to_excel(self.excel_filename)


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Uses NCBI's datasets v14 tool to download GCF and GCA genomes. 
    This tool is considered alpha and will undergo changes.  As Feb 2023 it seems to work well for bacteria genomes, but seemed unpredictable for downloads not consistent with the GCA format.

    GCA for GenBank assemblies and GCF for RefSeq assemblies.

    https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/

    GCA_000195835 # this is the type of number to feed script.
    NC_002945.4 # will NOT work, see vsnp3_download_fasta_gbk_gff_by_acc.py

    Usage:  module load datasets-14.12.0-gcc-9.2.0-46pggay
            fasta_GCA_get_metadata.py -a GCA_000195835 -gr

    Metadata will be downloaded along with the genome.  Metadata will be placed into an Excel file to be manipulated and used later.  Files can be immediately renamed with metadata info using --rename option.

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-a', '--accession', action='store', dest='accession', required=True, help='NCBI chromosome number -- GCA NUMBER, FASTA downloaded')
    parser.add_argument('-g', '--gfiles', action='store_true', dest='gfiles', help='download gff3 and gbk files')
    parser.add_argument('-r', '--rename', action='store_true', dest='rename', help='rename file with metadata')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    downloader = Downloader(gca_accession=args.accession, annotation_file=args.gfiles, rename=args.rename)

    # Excel Stats
    excel_stats = Excel_Stats(downloader.sample_name)
    downloader.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    # Clean up temporary files
    if os.path.exists(f'{args.accession}.zip'):
        os.remove(f'{args.accession}.zip')
    if os.path.exists('README.md'):
        os.remove('README.md')
    if os.path.exists('ncbi_dataset'):
        shutil.rmtree('ncbi_dataset')

    # Move all generated files to a dedicated directory
    files_grab = []
    for files in ('*.fasta', '*.gff', '*gbk', '*xlsx'):
        files_grab.extend(glob.glob(files))
    acc_dir = args.accession
    if not os.path.exists(acc_dir):
        os.makedirs(acc_dir)
    for each in files_grab:
        shutil.move(each, acc_dir)

# Created 2023 by Tod Stuber