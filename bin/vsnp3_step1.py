#!/usr/bin/env python

import os
import sys

from vsnp3_version import __version__

# --version is answered before anything heavy is imported.  pandas, Biopython
# and dask below take seconds to import from a slow network file system, and
# the vSNP GUI's provenance writer probes `vsnp3_step1.py --version` under a
# 10 s timeout, recording an empty version when the imports were slow.  The
# text is the one argparse prints for -v below, so nothing that reads it sees a
# difference.
if __name__ == "__main__" and any(arg in ('-v', '--version') for arg in sys.argv[1:]):
    print(f'{os.path.basename(__file__)}: version {__version__}')
    sys.exit(0)

import subprocess
import shutil
import glob
import re
import locale
import argparse
import textwrap
import importlib.metadata
import json
import platform

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import HtmlReport
from vsnp3_file_setup import Excel_Stats
from vsnp3_file_setup import fastq_usability_verdict, reference_usability_verdict
from vsnp3_input_validator import validate_file_inputs

from vsnp3_fastq_stats_seqkit import FASTQ_Stats
from vsnp3_best_reference_sourmash import Best_Reference
from vsnp3_reference_options import Ref_Options
from vsnp3_alignment_vcf import Alignment, concatenated_reference_name
from vsnp3_spoligotype import Spoligo
from vsnp3_version import NANOPORE_AVG_READ_LEN_CUTOFF
from vsnp3_zero_coverage import no_snp_notice, mixed_sample_notice
from vsnp3_group_reporter import GroupLookupError, GroupReporter, NoDefiningSnpWorkbook
from vsnp3_fasta_to_fastq import Fasta_to_Paired_Fastq

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class vSNP3_Step1(Setup):
    ''' 
    vSNP3 Step1 class for processing FASTQ files and generating VCF files
    '''
    def __init__(self, SAMPLE_NAME=None, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, gbk=None, reference_type=None, nanopore=False, assemble_unmap=None, spoligo=None, debug=False):
        '''
        Initialize vSNP3_Step1 class
        '''
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1)
        self.assemble_unmap = assemble_unmap
        self.spoligo = spoligo
        self.report = HtmlReport(self.sample_name)
        self.excel_stats = Excel_Stats(self.sample_name)
        self.nanopore = nanopore
        
        if FASTA: # IF -f REFERENCE FASTA PROVIDED USE IT, ie .fasta
            concatenated_FASTA = self.concat_fasta(FASTA) # -f option for FASTA will take wildcard for multiple FASTAs
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=gbk, debug=debug)
            self.reference_type = None
            with open(self.FASTA, encoding='utf-8', errors='replace') as f:
                self.top_header_description = f.readline()
                
        elif reference_type: # -t directory name
            reference_options = Ref_Options(reference_type)
            concatenated_FASTA = self.concat_fasta(reference_options.fasta)
            self.excel_stats.excel_dict["Reference"] = f'{reference_type} Forced'
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=reference_options.gbk, debug=debug)
            self.reference_type = reference_type
            with open(self.FASTA, encoding='utf-8', errors='replace') as f:
                self.top_header_description = f.readline()
                
        else: # IF NO REFERENCE PROVIDED SEEK A "BEST REFERENCE", for tb, brucella, paraTB or SARS-CoV-2
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2)
            self.best_reference = Best_Reference(FASTQ_R1=FASTQ_R1)
            self.best_reference.run()
            self.top_header_description = self.best_reference.top_fasta_header
            
            if self.best_reference.reference_set:
                reference_options = Ref_Options(self.best_reference.reference_set)
            else:
                print("\n" + "="*80)
                print("🔥💥 REFERENCE FAILURE! 💥🔥")
                print("="*80)
                print("❌ NO REFERENCE TYPE AVAILABLE!")
                print("❌ NO FASTA PROVIDED!")
                print("\n🔧 SOLUTION:")
                print("   Run: vsnp3_path_adder.py -s")
                print("   Check: reference_options_paths.txt")
                print("\n💡 This will show you available reference options!")
                print("="*80 + "\n")
                print(f'#### {self.sample_name} PROVIDE A REFERENCE, REFERENCE NOT FOUND')
                self.excel_stats.excel_dict["Reference"] = 'PROVIDE A REFERENCE, REFERENCE NOT FOUND'
                print('Getting Stats before exiting')
                
                fastq_stats = FASTQ_Stats(FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
                fastq_stats.run()
                fastq_stats.excel(self.excel_stats.excel_dict)
                self.report.add_notice(
                    'No reference could be determined, so no alignment was performed. '
                    'Only read statistics are reported.')
                # The ranking is the evidence for why nothing matched, so it is
                # worth more on this path than on a successful one.
                self.report.add_table('Sourmash reference similarity',
                                      ('Similarity', 'ID'),
                                      self.best_reference.similarity_rows())
                self.write_report()
                self.excel_stats.post_excel()
                
                temp_dir = './temp'
                if not os.path.exists(temp_dir):
                    os.makedirs(temp_dir)
                files_grab = []
                for files in ('*.log', '*out', "*_seqkit_stats.txt"):
                    files_grab.extend(sorted(glob.glob(files)))
                for each in files_grab:
                    shutil.move(each, temp_dir)

                if not debug:  # Use the debug parameter directly
                    shutil.rmtree(temp_dir)
                # No reference means no alignment and no VCF.  Exiting 0 told every
                # wrapper and scheduler this sample had succeeded.
                sys.exit(1)

            # Checked before concat_fasta, which iterates it.  The check used to sit
            # after the call, so a missing FASTA raised "'NoneType' object is not
            # iterable" from inside concat_fasta and this message was unreachable.
            if not reference_options.fasta:
                print('Sourmash unable to find a suitable match. Provide a FASTA (-f) or directory name with dependencies (-n)')
                sys.exit(1)

            concatenated_FASTA = self.concat_fasta(reference_options.fasta)
            self.gbk = reference_options.gbk
            self.excel_stats.excel_dict["Reference"] = f'{reference_options.select_ref} by Best Reference'
            Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=concatenated_FASTA, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, gbk=reference_options.gbk, debug=debug)
            self.best_reference.excel(self.excel_stats.excel_dict)
            self.report.add_table('Sourmash reference similarity',
                                  ('Similarity', 'ID'),
                                  self.best_reference.similarity_rows())
            self.reference_type = self.best_reference.reference_set

    def write_report(self):
        '''
        Write the report as HTML, then as a PDF rendered from that HTML.

        Both paths that finish a run come through here, so neither can end up
        with only one of the two files.
        '''
        html_path = self.report.write(self.excel_stats.excel_dict)
        print(f'\tReport: {bcolors.WHITE}{html_path}{bcolors.ENDC}')
        pdf_path = self.report.to_pdf()
        if pdf_path:
            print(f'\tReport PDF: {bcolors.WHITE}{pdf_path}{bcolors.ENDC}')
        return html_path, pdf_path

    def concat_fasta(self, FASTAs):
        """
        Concatenate multiple FASTA files into a single file
        """
        self.excel_stats.excel_dict["FASTA/s"] = ", ".join(os.path.basename(each) for each in FASTAs)

        # Concatenate FASTA list, under a name that is not a file already in
        # the working directory (the user's own reference given with -r).
        concatenated_FASTA = concatenated_reference_name(FASTAs)
        concatenated_TEMP = f'{concatenated_FASTA}.tmp'

        with open(concatenated_TEMP, 'wb') as wfd:
            for each in FASTAs:
                with open(each, 'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        shutil.move(concatenated_TEMP, concatenated_FASTA)
        return concatenated_FASTA

    def run(self):
        '''
        Main execution method for vSNP3 Step1 processing
        '''
        fastq_stats = FASTQ_Stats(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
        fastq_stats.run()
        fastq_stats.excel(self.excel_stats.excel_dict)

        # Test for Mycobacterium
        MYCO = None
        try:
            if any(x in self.top_header_description for x in ['Mycobacterium bovis', 'Mycobacterium tuberculosis', 'Mycobacterium caprae', 'Mycobacterium orygis']):
                MYCO = True
        except AttributeError:
            pass
        self.MYCO = MYCO
        
        if self.spoligo:
            spoligo = Spoligo(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, debug=self.debug)
            spoligo.spoligo()
            spoligo.excel(self.excel_stats.excel_dict)

        # Platform detection.  This tested the MAXIMUM read length, so a single long
        # read in an Illumina library rerouted the whole sample to minimap2 and
        # bcftools, whose QUAL was then incomparable to the thresholds (now the
        # long-read path writes the evidence QUAL; see vsnp3_longread.py).
        #
        # Average length is the right statistic - Illumina 2x300 averages ~300, ONT
        # well over 1000 - and the presence of an R2 file settles it outright,
        # because paired-end reads mean Illumina whatever the lengths look like.
        if self.nanopore:
            nanopore = True
            platform_reason = 'requested with --nanopore'
        else:
            avg_len = float(str(fastq_stats.R1.avg_len).replace(',', ''))
            if self.FASTQ_R2:
                nanopore = False
                platform_reason = (f'paired-end reads present, average length '
                                   f'{avg_len:.0f}')
            elif avg_len > NANOPORE_AVG_READ_LEN_CUTOFF:
                nanopore = True
                platform_reason = (f'single-end with average read length '
                                   f'{avg_len:.0f} > {NANOPORE_AVG_READ_LEN_CUTOFF}')
            else:
                nanopore = False
                platform_reason = (f'single-end with average read length '
                                   f'{avg_len:.0f} <= {NANOPORE_AVG_READ_LEN_CUTOFF}')
        platform = 'nanopore/long-read' if nanopore else 'Illumina/short-read'
        # Recorded because it was previously invisible, which is why the max_len
        # behaviour went unnoticed.
        print(f'\tPlatform: {platform} ({platform_reason})')
        self.excel_stats.excel_dict['Platform'] = f'{platform} ({platform_reason})'


        alignment = Alignment(SAMPLE_NAME=self.sample_name, FASTQ_R1=self.FASTQ_R1, FASTQ_R2=self.FASTQ_R2, reference=self.reference, nanopore=nanopore, gbk=self.gbk, assemble_unmap=self.assemble_unmap, debug=self.debug)
        alignment.run()
        self.tool_output = list(getattr(alignment, 'tool_output', None) or [])

        # Group reporting.  Three outcomes, previously collapsed into two: a
        # failed lookup was reported as "group file not provided", which reads as
        # a configuration choice rather than an error.  A reference type with no
        # defining-SNP workbook is that configuration choice, and says so.
        # Everything else is caught: this step only reports, and the _zc.vcf is
        # already written, so no exception from it may end the run.
        if not self.reference_type:
            groups = "group file not provided"
        else:
            try:
                group_reporter = GroupReporter(alignment.zero_coverage_vcf_file_path, self.reference_type)
                groups = ", ".join(group_reporter.get_groups())
            except NoDefiningSnpWorkbook:
                groups = "group file not provided"
            except Exception as e:                                # noqa: BLE001
                groups = f"group lookup FAILED: {str(e) or type(e).__name__}"
                print(f'\n### {groups}\n')
        self.excel_stats.excel_dict['Groups'] = groups

        alignment.excel(self.excel_stats.excel_dict)

        # Usability verdicts.  The limits live in vsnp3_version.py, and the
        # verdict is reported with the measurement that set it, so that
        # "Questionable" is never a bare word.  Shown in the HTML report only:
        # vsnp3_file_setup.HTML_ONLY_KEYS keeps them out of the stats workbook.
        zc = alignment.zero_coverage
        fastq_usability, fastq_reason = fastq_usability_verdict(
            depth=zc.ave_coverage, duplicate_ratio=alignment.DUPLICATION_RATIO,
            q20_percent=fastq_stats.R1.passQ20, q30_percent=fastq_stats.R1.passQ30)
        colors = {'Poor': bcolors.RED, 'Questionable': bcolors.YELLOW}
        color_quality = colors.get(fastq_usability, bcolors.GREEN)
        print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}{fastq_usability}{bcolors.ENDC} FASTQ Usability ({fastq_reason})')
        self.excel_stats.excel_dict['FASTQ Usability'] = f'{fastq_usability} ({fastq_reason})'

        reference_usability, reference_reason = reference_usability_verdict(
            genome_coverage=zc.genome_coverage,
            zero_coverage_percent=zc.percent_ref_with_zero_coverage,
            good_snp_percent=zc.percent_ref_with_good_snp_count)
        reference_usability_test = reference_usability == 'Acceptable'
        color_quality = colors.get(reference_usability, bcolors.GREEN)
        print(f'{bcolors.WHITE}{self.fastq_name}{bcolors.ENDC} {color_quality}{reference_usability}{bcolors.ENDC} Reference Usability ({reference_reason})')
        self.excel_stats.excel_dict['Reference Usability'] = f'{reference_usability} ({reference_reason})'

        # A sample with no SNP passing the step 2 thresholds has an empty row
        # in step 2.  Say why now, while the reads are at hand: mixed (AC=1)
        # calls are the signature of a mixed or contaminated sample, and no
        # call at all is a sample identical to the reference.  Decided on the
        # step 2 definition (zc.step2_snp_count), not on "Quality SNPs", which
        # counts any AC: a mixed sample with high-quality AC=1 calls has many
        # quality SNPs, so this notice never appeared for it.
        notice = no_snp_notice(getattr(zc, 'step2_snp_count', None),
                               getattr(zc, 'ac1_count', 0), getattr(zc, 'ac2_snp_count', 0))
        # A mixture can also have a few passing SNPs: a site where one strain's
        # reads happen to dominate is called homozygous.  It is still flagged
        # when its confident mixed calls outnumber those SNPs several times over
        # (vsnp3_version.mixed_signature).
        if not notice:
            notice = mixed_sample_notice(getattr(zc, 'step2_snp_count', None),
                                         getattr(zc, 'confident_ac1_count', None))
        if notice:
            print(f'\n### {self.sample_name}: {notice}\n')
            self.report.add_notice(notice)
        
        # Store alignment results for reporting
        self.programs = alignment.programs
        self.alignment_vcf_run_summary = alignment.alignment_vcf_run_summary
        if str(groups).startswith('group lookup FAILED'):
            self.report.add_notice(groups)
        self.report.program_versions = "\n".join(
            list(self.programs or []) + list(self.alignment_vcf_run_summary or []))
        self.write_report()
        self.excel_stats.post_excel()

        # Clean up copied gbk files if they exist
        if hasattr(self, 'remove_copied_gbk') and self.remove_copied_gbk:
            for each in self.gbk:
                if os.path.exists(each):
                    os.remove(each)


def get_system_info():
    """
    Collect comprehensive system information for logging
    """
    system_info = []
    
    # Operating System details
    os_name = platform.system()
    os_version = platform.version()
    os_release = platform.release()
    
    if os_name == "Darwin":  # macOS
        # Check if running on Apple Silicon (ARM) or Intel
        arch = platform.machine()
        cpu_type = "ARM" if arch == "arm64" else "Intel"
        
        # Get macOS version details using sw_vers
        try:
            mac_version = subprocess.check_output(["sw_vers", "-productVersion"], 
                                            universal_newlines=True).strip()
            system_info.append(f"Operating System: macOS {mac_version} ({cpu_type})")
        except:
            system_info.append(f"Operating System: macOS {os_version} ({cpu_type})")
    
    elif os_name == "Linux":
        # Get CPU architecture
        arch = platform.machine()
        
        # Try to get Linux distribution info
        try:
            # Try reading from /etc/os-release first
            if os.path.exists('/etc/os-release'):
                with open('/etc/os-release', 'r', encoding='utf-8', errors='replace') as f:
                    lines = f.readlines()
                    distro_id = ""
                    distro_version = ""
                    for line in lines:
                        if line.startswith('ID='):
                            distro_id = line.split('=')[1].strip().strip('"')
                        elif line.startswith('VERSION_ID='):
                            distro_version = line.split('=')[1].strip().strip('"')
                    
                    if distro_id and distro_version:
                        system_info.append(f"Operating System: Linux {distro_id.capitalize()} {distro_version} ({arch})")
                    else:
                        system_info.append(f"Operating System: Linux {os_release} ({arch})")
            else:
                # Fallback to lsb_release if available
                lsb_output = subprocess.check_output(["lsb_release", "-ds"], 
                                                universal_newlines=True).strip()
                system_info.append(f"Operating System: {lsb_output} ({arch})")
        except:
            # Fallback to basic platform info
            system_info.append(f"Operating System: Linux {os_release} ({arch})")
            
        # Check if running on HPC (look for typical HPC indicators)
        try:
            # Check for common HPC environment variables or files
            if any(os.environ.get(var) for var in ['SLURM_JOB_ID', 'PBS_JOBID', 'SGE_TASK_ID']):
                hpc_type = "Unknown HPC"
                if os.environ.get('SLURM_JOB_ID'):
                    hpc_type = "SLURM"
                elif os.environ.get('PBS_JOBID'):
                    hpc_type = "PBS/Torque"
                elif os.environ.get('SGE_TASK_ID'):
                    hpc_type = "SGE"
                
                system_info.append(f"HPC Environment: {hpc_type}")
        except:
            pass
    
    elif os_name == "Windows":
        win_version = platform.win32_ver()[0]
        arch = platform.machine()
        system_info.append(f"Operating System: Windows {win_version} ({arch})")
    
    else:
        # Generic fallback
        system_info.append(f"Operating System: {os_name} {os_version} ({platform.machine()})")
    
    # CPU information
    try:
        import cpuinfo
        cpu_info = cpuinfo.get_cpu_info()
        system_info.append(f"CPU: {cpu_info['brand_raw']} ({cpu_info['count']} cores)")
    except:
        # Fallback CPU info if py-cpuinfo is not available
        try:
            if os_name == "Linux":
                with open('/proc/cpuinfo', 'r', encoding='utf-8', errors='replace') as f:
                    cpu_lines = f.readlines()
                model_name = "Unknown"
                cpu_count = 0
                
                for line in cpu_lines:
                    if line.startswith('model name'):
                        model_name = line.split(':')[1].strip()
                    if line.startswith('processor'):
                        cpu_count += 1
                
                system_info.append(f"CPU: {model_name} ({cpu_count} cores)")
            elif os_name == "Darwin":
                sysctl_output = subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"], 
                                                    universal_newlines=True).strip()
                cpu_count = subprocess.check_output(["sysctl", "-n", "hw.ncpu"], 
                                                universal_newlines=True).strip()
                system_info.append(f"CPU: {sysctl_output} ({cpu_count} cores)")
            else:
                system_info.append(f"CPU: {platform.processor()}")
        except:
            system_info.append(f"CPU: {platform.processor()}")
    
    return system_info


def get_program_versions():
    """
    Get version information for all required programs and dependencies
    """
    program_list = [
        'bcftools', 'biopython', 'bwa', 'minimap2',
        'dask', 'freebayes', 'humanize', 'numpy', 'pandas', 'openpyxl', 
        'xlsxwriter', 'parallel', 'pigz', 'regex', 'samtools', 'seqkit', 
        'sourmash', 'spades', 'py-cpuinfo', 'jinja2', 'weasyprint'
    ]

    # Dictionary for module name mapping (conda/pip name to import name)
    module_mapping = {
        'biopython': 'Bio',
        'numpy': 'numpy',
        'pandas': 'pandas',
        'openpyxl': 'openpyxl',
        'xlsxwriter': 'xlsxwriter',
        'regex': 're',
        'py-cpuinfo': 'cpuinfo',
        'dask': 'dask',
        'humanize': 'humanize',
        'jinja2': 'jinja2',
        'weasyprint': 'weasyprint',
    }

    # List to store program versions
    python_programs = []

    # One `conda list --json` for the whole environment, rather than one
    # `conda list <pkg>` per package.  Each of those calls took ~0.42 s and there
    # were 22 of them -- 9.3 s of pure overhead per sample, or roughly 74 minutes
    # across a 500-sample run -- to obtain information a single call already holds.
    conda_versions = {}
    try:
        conda_json = subprocess.check_output(
            ["conda", "list", "--json"], stderr=subprocess.DEVNULL,
            universal_newlines=True, timeout=120)
        for entry in json.loads(conda_json):
            name = entry.get('name')
            if name:
                conda_versions[name.lower()] = entry.get('version', 'nd')
    except (subprocess.SubprocessError, OSError, ValueError):
        # Not in a conda environment, or conda is unavailable.  The per-program
        # fallbacks below still report versions.
        pass

    python_source = "(conda)" if 'python' in conda_versions else "(system)"
    python_programs.append(f'Python, {sys.version.split()[0]} {python_source}')

    # Check all programs
    for program in program_list:
        version = "nd"  # Default to "nd" (no data)
        source = ""

        # First try conda
        if program.lower() in conda_versions:
            version = conda_versions[program.lower()]
            source = "(conda)"
            python_programs.append(f'{program}, {version} {source}')
            continue
        
        # If not found in conda, check Python modules
        if program in module_mapping:
            module_name = module_mapping[program]
            try:
                # Try importlib.metadata first
                try:
                    version = importlib.metadata.version(program)
                except (importlib.metadata.PackageNotFoundError, ModuleNotFoundError):
                    # Fall back to module's __version__ attribute
                    if module_name in sys.modules and hasattr(sys.modules[module_name], '__version__'):
                        version = sys.modules[module_name].__version__
                    elif module_name not in sys.modules:
                        # Try to import the module
                        try:
                            module = __import__(module_name)
                            if hasattr(module, '__version__'):
                                version = module.__version__
                        except ImportError:
                            pass
                # If we found the version, add it and continue
                if version != "nd":
                    source = "(system)"
                    python_programs.append(f'{program}, {version} {source}')
                    continue
            except Exception:
                pass
        
        # Finally, check command-line programs
        if shutil.which(program):
            try:
                # Different programs report versions differently
                if program in ['bcftools', 'samtools', 'bwa']:
                    cmd_output = subprocess.check_output([program, "--version"], 
                                                    stderr=subprocess.STDOUT, 
                                                    universal_newlines=True)
                    version = cmd_output.strip().split('\n')[0].split(' ')[1]
                elif program in ['minimap2', 'spades.py', 'seqkit']:
                    cmd_output = subprocess.check_output([program, "--version"], 
                                                    stderr=subprocess.STDOUT, 
                                                    universal_newlines=True)
                    version = cmd_output.strip()
                elif program == 'freebayes':
                    cmd_output = subprocess.check_output([program, "--version"], 
                                                    stderr=subprocess.STDOUT, 
                                                    universal_newlines=True)
                    version = cmd_output.strip().split(',')[0].split(' ')[-1]
                elif program in ['pigz', 'parallel']:
                    cmd_output = subprocess.check_output([program, "--version"], 
                                                    stderr=subprocess.STDOUT, 
                                                    universal_newlines=True)
                    version = cmd_output.strip().split('\n')[0]
                
                if version != "nd":
                    source = "(system)"
                    python_programs.append(f'{program}, {version} {source}')
            except Exception:
                pass

    return python_programs


def write_run_log(sample_name, args, alignment_vcf_run_summary, programs, tool_output=None):
    """
    Write comprehensive run log with system info and program versions.

    tool_output: [(title, [lines])] of external tool output kept by the
    alignment (a failed SPAdes run's last lines, freebayes-parallel's stderr),
    written at the end under "Tool output".  Not part of the HTML report.
    """
    # Get system information
    system_info = get_system_info()
    
    # Get program versions
    python_programs = get_program_versions()
    
    # Check if running in a conda environment
    conda_env = os.environ.get('CONDA_DEFAULT_ENV')
    
    try:
        log_filename = f'{sample_name}_run_log.txt'
        # Explicit encoding: the log carries tool command lines with the
        # reference and reads paths in them, and under LC_ALL=C on Linux the
        # default encoding is ASCII, so an accent in a path ended the run here.
        with open(log_filename, 'w', encoding='utf-8', errors='replace') as run_log:
            print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:', file=run_log)
            print(args, file=run_log)
            
            if hasattr(args, 'FASTAtoFASTQ') and args.FASTAtoFASTQ:
                print('Converted FASTA to FASTQ', file=run_log)
            
            print('\nCall Summary:', file=run_log)
            for each in alignment_vcf_run_summary:
                print(each, file=run_log)
            
            print('\nSystem Information:', file=run_log)
            for info in system_info:
                print(info, file=run_log)
            
            print('\nVersions:', file=run_log)
            
            # Display conda environment info if applicable
            if conda_env:
                print(f'Using conda environment: {conda_env}', file=run_log)
            
            print(f'vSNP3: {__version__}', file=run_log)
            
            for each in python_programs:
                print(each, file=run_log)
            
            # Skip programs we already reported
            programs_already_reported = [item.split(',')[0].lower() for item in python_programs]
            
            for each in programs:
                # The filter that used to live here suppressed
                # "Error getting program versions: cannot access local variable",
                # which was a real UnboundLocalError in vsnp3_alignment_vcf.py on
                # every run that did not use --assemble_unmap.  The cause is fixed,
                # so a message like that now means something and is not hidden.
                program_name = each.split(':')[0].strip() if ':' in each else each.split()[0].strip()
                if program_name.lower() not in programs_already_reported:
                    print(each, file=run_log)

            if tool_output:
                print('\nTool output:', file=run_log)
                for title, lines in tool_output:
                    print(f'\n{title}:', file=run_log)
                    for line in lines:
                        print(f'  {line}', file=run_log)

    except Exception as e:
        print(f"Error writing to log file: {str(e)}")
        # Fallback to a default filename
        try:
            with open("run_log.txt", 'w', encoding='utf-8', errors='replace') as run_log:
                print('\nSystem Information:', file=run_log)
                for info in system_info:
                    print(info, file=run_log)
                
                print('\nVersions:', file=run_log)
                if conda_env:
                    print(f'Using conda environment: {conda_env}', file=run_log)
                for each in python_programs:
                    print(each, file=run_log)
        except Exception as fallback_error:
            print(f"Error writing fallback log file: {str(fallback_error)}")


def cleanup_temp_files(debug=False):
    """
    Clean up temporary files generated during processing
    """
    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
        
    files_grab = []
    for files in ('*_report.out', '*.log', "*_seqkit_stats.txt"):
        files_grab.extend(sorted(glob.glob(files)))
    
    for each in files_grab:
        try:
            shutil.move(each, temp_dir)
        except (shutil.Error, OSError) as e:
            print(f"Warning: Could not move {each} to temp directory: {e}")

    if not debug:
        try:
            shutil.rmtree(temp_dir)
        except OSError as e:
            print(f"Warning: Could not remove temp directory: {e}")


if __name__ == "__main__": # execute if directly accessed by the interpreter
    parser = argparse.ArgumentParser(
        prog='vsnp3_step1.py', 
        formatter_class=argparse.RawDescriptionHelpFormatter, 
        description=textwrap.dedent(f'''\

    ---------------------------------------------------------
    Conda:
    conda install vsnp3={__version__} -c conda-forge -c bioconda
    ---------------------------------------------------------

    When running samples through step1 and 2 of vSNP, or when running a routine analysis, 
    set up dependencies using vsnp3_path_adder.py

    See vsnp3_path_adder.py -h for adding a reference type and for more information

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -t NC_045512_wuhan-hu-1

    For a few reference types (not all) sourmash can be used to find the best reference.  
    See vsnp3_best_reference_sourmash.py -h for list of references.

    If working or suspect one of these reference types the dependencies will automatically be found

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz

    When running a one-off dependencies can be provided explicitly: -r is the reference
    FASTA and -b its GenBank file.  (-f is different: it converts a FASTA into synthetic
    reads and analyses those.)

    Usage example:
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -r *fasta -b *gbk
    or
    vsnp3_step1.py -r1 ../myfastqs/*_R1*.fastq.gz -r2 ../myfastqs/*_R2*.fastq.gz -r ../myreference/*fasta -b ../myreference/*gbk

    At completion of running script it is recommended to store sample folders by reference.  
    For each reference there should be a "step1" and "step2" subfolder.  These sample folders 
    should go into step1 subfolders for the reference type used to generate SNPs stored in the VCF file.
    ---------------------------------------------------------

    '''), 
    epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=False, help='Provide R1 FASTQ gz file.  A single read file can also be supplied to this option')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: provide R2 FASTQ gz file')
    parser.add_argument('-f', '--FASTAtoFASTQ', action='store', dest='FASTAtoFASTQ', required=False, help='Input a FASTA file to be analysed as a sample: it is converted to synthetic paired FASTQ files and run.  This is not the reference; the reference FASTA is -r.')
    parser.add_argument('-r', '--FASTA', nargs='*', dest='FASTA', required=False, help='FASTA file to be used as reference.  Multiple can be specified with wildcard')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', required=False, default=None, help="Optional: Provide directory name with FASTA and GBK file/s")
    parser.add_argument('-p', '--nanopore', action='store_true', dest='nanopore', default=False, help='nanopore (long-read) reads: minimap2 alignment, bcftools calling and the vSNP3 long-read evidence QUAL, so the VCF is read by the same thresholds as an Illumina one.  Chosen automatically for single-end reads averaging more than 701 bases; this flag forces it')
    parser.add_argument('-o', '--output_dir', action='store', dest='output_dir', required=False, default=None, help="Optional: Provide a name.  This name will be a directory output files are written to.  Name can be a directory path, but doesn't have to be.")
    parser.add_argument('-assemble_unmap', '--assemble_unmap', action='store_true', dest='assemble_unmap', help='Optional: assemble the unmapped reads with SPAdes (not done by default).   See also vsnp3_assembly.py')
    parser.add_argument('-spoligo', '--spoligo', action='store_true', dest='spoligo', help='Optional: get spoligotype if TB complex.  See also vsnp3_spoligotype.py')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='keep spades output directory')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    if not args.FASTQ_R1 and not args.FASTAtoFASTQ:
        parser.error('No FASTQ or FASTA provided.  Use -r1 or -f option')

    # -f with -b and no reference is only ever the reference FASTA given under
    # the wrong flag (older help text and wrappers showed "-f *fasta -b *gbk").
    # Left alone, it converts the reference genome to reads and aligns them to
    # whichever reference sourmash picks, with no error.
    if args.FASTAtoFASTQ and args.gbk and not (args.FASTA or args.reference_type):
        parser.error('-f/--FASTAtoFASTQ converts a FASTA into reads to analyse as a sample; '
                     'it is not the reference.  With -b/--gbk, give the reference FASTA '
                     'with -r/--FASTA (or a reference type with -t).')

    # Handle FASTA to FASTQ conversion
    if args.FASTAtoFASTQ:
        # Validate FASTA file before conversion
        print(f"🔍 Validating input FASTA file: {args.FASTAtoFASTQ}")
        is_valid, errors = validate_file_inputs(fasta=[args.FASTAtoFASTQ], debug=args.debug)
        if not is_valid:
            print(f"\n❌ FASTA validation failed:")
            for error in errors:
                print(f"   {error}")
            print("\n💡 Please check your FASTA file and try again.")
            sys.exit(1)
        print("✅ FASTA validation passed")

        fasta_to_paired_fastq = Fasta_to_Paired_Fastq(args.FASTAtoFASTQ, coverage=100, read_length=300)
        args.FASTQ_R1 = fasta_to_paired_fastq.fastq_r1_file
        args.FASTQ_R2 = fasta_to_paired_fastq.fastq_r2_file
        print('Conversion to FASTQ completed')

    # COMPREHENSIVE INPUT VALIDATION
    print("\n🔍 VALIDATING INPUT FILES...")
    print("="*50)

    is_valid, errors = validate_file_inputs(
        fastq_r1=args.FASTQ_R1,
        fastq_r2=args.FASTQ_R2,
        fasta=args.FASTA,
        gbk_list=args.gbk,
        debug=args.debug
    )

    if not is_valid:
        print(f"\n❌ INPUT VALIDATION FAILED:")
        print("-" * 30)
        for error in errors:
            print(f"   {error}")

        print(f"\n💡 SOLUTIONS:")
        print("   • Check that all file paths are correct")
        print("   • Verify files exist and are readable")
        print("   • Ensure FASTQ files are properly formatted")
        print("   • Confirm FASTA files have valid sequences")
        print("="*50)
        sys.exit(1)

    print("✅ All input files validated successfully")
    print("="*50 + "\n")

    # Handle output directory setup
    if args.output_dir:
        output_dir = args.output_dir
        output_dir = os.path.expanduser(output_dir)
        output_dir = os.path.abspath(output_dir)
        args.output_dir = output_dir
        # Every path argument absolute before the chdir below.  A relative
        # -r ../ref/x.fasta or -b ../ref/x.gbk passed validation from the
        # original directory and was then looked up from output_dir, where
        # concat_fasta raised FileNotFoundError.
        for name in ('FASTQ_R1', 'FASTQ_R2', 'FASTAtoFASTQ'):
            value = getattr(args, name, None)
            if value:
                setattr(args, name, os.path.abspath(os.path.expanduser(value)))
        for name in ('FASTA', 'gbk'):
            value = getattr(args, name, None)
            if value:
                setattr(args, name, [os.path.abspath(os.path.expanduser(each)) for each in value])
        os.makedirs(output_dir, exist_ok=True)

        if args.FASTQ_R1:
            shutil.copy(args.FASTQ_R1, output_dir)
        if args.FASTQ_R2:
            shutil.copy(args.FASTQ_R2, output_dir)
        os.chdir(output_dir)

    # Main execution
    try:
        vsnp = vSNP3_Step1(
            SAMPLE_NAME=args.SAMPLE_NAME, 
            FASTQ_R1=args.FASTQ_R1, 
            FASTQ_R2=args.FASTQ_R2, 
            FASTA=args.FASTA, 
            gbk=args.gbk, 
            reference_type=args.reference_type, 
            nanopore=args.nanopore, 
            assemble_unmap=args.assemble_unmap, 
            spoligo=args.spoligo,
            debug=args.debug
        )
        vsnp.run()
        
        # Write comprehensive log file
        write_run_log(vsnp.sample_name, args, vsnp.alignment_vcf_run_summary, vsnp.programs,
                      tool_output=getattr(vsnp, 'tool_output', None))
        
    except Exception as e:
        print(f"Error during vSNP3 Step1 execution: {str(e)}")
        sys.exit(1)
    
    finally:
        # Clean up temporary files
        cleanup_temp_files(args.debug)

    print(f"\n{os.path.basename(__file__)} completed successfully!")

# Created 2021 by Tod Stuber