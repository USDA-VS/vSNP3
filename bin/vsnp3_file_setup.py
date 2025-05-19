#!/usr/bin/env python

__version__ = "3.29"

import os
import shutil
import re
import sys
import locale
import pandas as pd
import multiprocessing
from datetime import datetime
import svgwrite
from cairosvg import svg2png
import yaml

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE='\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

class Setup:
    ''' 
    Standarize setup
    '''

    def __init__(self, SAMPLE_NAME=None, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, reference=None, gbk=None, debug=False):
        self.cwd = os.getcwd()
        self.FASTA = FASTA
        self.FASTQ_R1 = FASTQ_R1
        self.FASTQ_R2 = FASTQ_R2
        self.reference = reference
        self.gbk = gbk
        self.debug = debug
        self.remove_copied_gbk = False
        if FASTQ_R2:
            self.paired = True
        else:
            self.paired = False
        if FASTA:
            try: #IF FASTA provided as path variable copy local
                shutil.copy(FASTA, self.cwd)
            except shutil.SameFileError:
                pass
            FASTA = os.path.basename(FASTA)
            self.FASTA = f'{self.cwd}/{FASTA}'
            self.sample_name = re.sub('[_.].*', '', FASTA)
            self.fasta_name = re.sub('[.].*', '', FASTA) #explict FASTA name if needed
        if gbk:
            try:
                updated_gbk_list=[]
                for each in gbk:
                    shutil.copy(each, self.cwd)
                    updated_gbk_list.append(f'{self.cwd}/{os.path.basename(each)}')
                self.gbk = updated_gbk_list
                self.remove_copied_gbk = True
            except shutil.SameFileError:
                self.gbk = gbk
                self.remove_copied_gbk = False
        if reference:
            self.reference_name = re.sub('[.].*', '', reference)
            self.reference = reference
            try: #IF FASTA provided as path variable copy local
                shutil.copy(reference, self.cwd)
            except shutil.SameFileError:
                pass
            reference = os.path.basename(reference)
            self.reference = f'{self.cwd}/{reference}'
        elif FASTA: # if no explict reference make the FASTA a "reference"
            self.reference_name = re.sub('[.].*', '', FASTA)
            self.reference = self.FASTA
        self.FASTQ_list=[]
        self.FASTQ_dict={}
        try:
            try: #IF FASTQ provided as path variable copy local
                shutil.copy(FASTQ_R1, self.cwd)
            except shutil.SameFileError:
                pass
            FASTQ_R1 = os.path.basename(FASTQ_R1)
            self.FASTQ_R1 = f'{self.cwd}/{FASTQ_R1}'
            self.sample_name = re.sub('[_.].*', '', FASTQ_R1) # default to FASTQ_R1 as sample name if also FASTA
            self.fastq_name = re.sub('[_.].*', '', FASTQ_R1) #explict FASTQ name if needed
            self.FASTQ_list.append(self.FASTQ_R1)
            self.FASTQ_dict = {'FASTQ_R1': self.FASTQ_R1}
        except TypeError:
            self.FASTQ_R1 = None
        if FASTQ_R2:
            try: #IF FASTQ provided as path variable copy local
                shutil.copy(FASTQ_R2, self.cwd)
            except shutil.SameFileError:
                pass
            FASTQ_R2 = os.path.basename(FASTQ_R2)
            self.FASTQ_R2 = f'{self.cwd}/{FASTQ_R2}'
            self.FASTQ_list.append(self.FASTQ_R2)
            self.FASTQ_dict['FASTQ_R2'] = self.FASTQ_R2
        else:
            self.FASTQ_R2 = None

        self.startTime = datetime.now()
        self.cpus = multiprocessing.cpu_count() - 2
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        if SAMPLE_NAME:
            self.sample_name = SAMPLE_NAME

    def print_time(self,):
        '''
        description
        '''
        self.run_time = datetime.now() - self.startTime
        print (f'\n\nruntime: {datetime.now() - self.startTime}\n')

    def print_run_time(self, tool):
        print(f'{bcolors.RED}\n{tool} {bcolors.ENDC}')
        now = datetime.now()
        print (f'{bcolors.WHITE}{now.strftime("%Y-%m-%d %H:%M:%S")}{bcolors.ENDC}')

class Banner:
    ''' 
    '''

    def __init__(self, title, hexcode="56, 68, 117"):
        '''
        Banner Bar
        '''
        width = 2600
        height = 90
        svgimg = svgwrite.Drawing(size=(width, height))
        svgimg.add(svgimg.rect([0, 0], [width, height], rx=None, ry=None, fill=f'rgb({hexcode})', stroke="black"))
        #text to banner
        svgimg.add(svgimg.text(title, insert=(30, 60), fill='white', font_size='50px', font_weight='bold'))
        svgimg.saveas('banner.svg')
        with open('banner.svg', 'r') as content_file:
            content = content_file.read()
        title = title.replace(" ", "_")
        svg2png(bytestring=content, write_to=f'{title}-banner.png')
        os.remove('banner.svg')
        self.banner = f'{os.getcwd()}/{title}-banner.png'
    
class Latex_Report:
    '''
    LaTeX report generator with configurable logo from YAML file
    '''
    def __init__(self, sample_name, config_file=None, report_description=None):
        # If no config file is provided, look for one in the dependencies directory
        if config_file is None:
            # Get the directory where the script is located (bin folder)
            script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
            # Go up one level to parent directory
            parent_dir = os.path.dirname(script_dir)
            # Look for dependencies folder at the same level as bin
            dependencies_dir = os.path.join(parent_dir, 'dependencies')
            # Default config filename in dependencies folder
            default_config_path = os.path.join(dependencies_dir, 'latex_config.yaml')
            if os.path.exists(default_config_path):
                config_file = default_config_path
        
        # Load configuration from YAML if provided
        self.config = self._load_config(config_file)
        
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        tex_file = f'{sample_name}_{date_stamp}_report.tex'
        self.tex_file = tex_file
        tex = open(tex_file, 'w')
        print(r'\documentclass{article}', file=tex)
        print(r'\usepackage[margin=1in]{geometry}', file=tex)
        print(r'\usepackage{adjustbox}', file=tex)
        print(r'\usepackage{graphicx}', file=tex)
        print(r'\usepackage{fancyhdr}', file=tex)
        print(r'\usepackage{hyperref}', file=tex)
        print(r'\pagestyle{fancy}', file=tex)
        print(r'\usepackage[scaled]{helvet}', file=tex)
        print(r'\renewcommand\familydefault{\sfdefault}', file=tex)
        print(r'\usepackage[T1]{fontenc}', file=tex)
        print(r'\usepackage{xcolor}', file=tex)
        
        # Add logo only if specified in config and file exists
        logo_path = self.config.get('logo_path')
        if logo_path and os.path.exists(logo_path):
            print(r'\lhead{\begin{picture}(0,0) \put(0,0){\includegraphics[width=40mm]{' + logo_path + r'}} \end{picture}}', file=tex)
        
        current_date = datetime.now().strftime('%B %d, %Y')
        print(r'\rhead{\begin{date} \textbf{\large ' + current_date + r'} \end{date}}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth} \\', file=tex)
        if report_description:
            print(r'\centerline{\textbf{\large{' + report_description + r'}}}', file=tex)
        print(r'\end{adjustbox} \\', file=tex)
        print(r'\vspace{0.5 mm} \\', file=tex)
        print(r'\textbf{\large {\fontfamily{\sfdefault}\selectfont Sample: ' + sample_name + r'}}', file=tex)
        print(r'\vspace{0.5 mm}', file=tex)
        print(r'\begin{document}', file=tex)
        self.tex = tex
    
    def _load_config(self, config_file):
        """Load configuration from YAML file or use defaults"""
        default_config = {
            'logo_path': None  # Default is no logo
        }
        
        if config_file and os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    user_config = yaml.safe_load(f)
                    if user_config:
                        default_config.update(user_config)
            except Exception as e:
                print(f"Error loading config file: {e}")
        
        return default_config
    
    def latex_ending(self):
        print(r'\end{document}', file=self.tex)
        self.tex.close()
        os.system(f'pdflatex -interaction=nonstopmode {self.tex_file} > /dev/null 2>&1')
        os.system(f'pdflatex -interaction=nonstopmode {self.tex_file} > /dev/null 2>&1')

class Excel_Stats:

    def __init__(self, sample_name):
        self.sample_name = sample_name
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.excel_filename = f'{sample_name}_{date_stamp}_stats.xlsx'
        excel_dict = {}
        excel_dict['sample'] = sample_name
        excel_dict['date'] = date_stamp
        self.excel_dict = excel_dict 

    def post_excel(self,):
        df = pd.DataFrame.from_dict(self.excel_dict, orient='index').T
        df = df.set_index('sample')
        df.to_excel(self.excel_filename)

# Created 2021 by Tod Stuber