#!/usr/bin/env python

__version__ = "3.11"

import os
import sys
import shutil
import re
import pickle
import argparse
import textwrap
import zipfile
import glob
from datetime import datetime
import allel
from collections import defaultdict
from concurrent import futures
import multiprocessing
multiprocessing.set_start_method('spawn', True)

import warnings
warnings.filterwarnings('ignore')

from vsnp3_file_setup import Setup
from vsnp3_group_on_defining_snps import Group
from vsnp3_reference_options import Ref_Options
from vsnp3_remove_from_analysis import Remove_From_Analysis

global_date_stamp=None
global_working_dir='.'

class VCF_to_DF():
    ''' 
    '''

    def __init__(self, vcf_list=None, debug=False): #write_out=False, 
        '''
        Start at class call
        '''
        self.startTime = datetime.now()
        self.vcf_bad_list=[]
        cpu_count = int(multiprocessing.cpu_count() / 1.2)
        dataframes={}
        self.vcf_original_count = len(vcf_list)
        if debug:
            print(f'VCF file count {self.vcf_original_count}')
            for vcf in vcf_list:
                print(vcf)
                vcf, df, vcf_bad_list_temp = self.check_and_fix(vcf)
                try:
                    self.chrom = df['CHROM'].iloc[0]
                except TypeError:
                    pass
                if df is not None:
                    dataframes[os.path.basename(vcf)] = df
                self.vcf_bad_list = self.vcf_bad_list + vcf_bad_list_temp
        else:
            print(f'Fixing: Pool processing with {cpu_count} cpus...')
            with futures.ProcessPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## process works best for calling on multiple files
                for vcf, df, vcf_bad_list_temp in pool.map(self.check_and_fix, vcf_list):
                    try:
                        self.chrom = df['CHROM'].iloc[0]
                    except TypeError:
                        pass
                    if df is not None:
                        dataframes[os.path.basename(vcf)] = df
                    self.vcf_bad_list = self.vcf_bad_list + vcf_bad_list_temp
        self.dataframes = dataframes
        print (f'\n\nDictionary of dataframes to memory runtime: {datetime.now() - self.startTime}\n')
        # if write_out: # write out pickle file can be used in downstream applications
        with open('dictionary_of_dataframes.pickle', 'wb') as handle:
            pickle.dump(dataframes, handle, protocol=pickle.HIGHEST_PROTOCOL)
        if not debug:
            os.remove('dictionary_of_dataframes.pickle')

    def check_and_fix(self, vcf):
        vcf_bad_list_temp = []
        try:
            self.vcf_fix(vcf)
            df = allel.vcf_to_dataframe(vcf, fields=['variants/CHROM', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/ALT', 'variants/AC', 'variants/DP', 'variants/MQ'], alt_number=1)
            df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
        except RuntimeError:
            self.vcf_fix(vcf)
            try:
                df = allel.vcf_to_dataframe(vcf, fields=['variants/CHROM', 'variants/POS', 'variants/QUAL', 'variants/REF', 'variants/ALT', 'variants/AC', 'variants/DP', 'variants/MQ'], alt_number=1)
                df['abs_pos'] = df['CHROM'] + ':' + df['POS'].astype(str)
            except:
                vcf_bad_list_temp.append(vcf)
                os.remove(vcf)
                df = None
        try:
            df = df.drop_duplicates(subset=['abs_pos'])
        except AttributeError:
            # pass if df is empty, NoneType
            pass
        return vcf, df, vcf_bad_list_temp

    def vcf_fix(self, vcf):
        temp_file = vcf + ".temp"
        write_out = open(temp_file, 'w') #r+ used for reading and writing to the same file
        initial_file_time_stats = os.stat(vcf)
        with open(vcf, 'r') as file:
            for line in file:
                line = line.replace('\r\n', '\n')
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip()  #remove right white space
                    line = re.sub(r';MQM=', r';MQ=', line) #Allow Freebayes MQM to be read as MQ.  MQ is VCF standard
                    line = re.sub(r'ID=MQM,', r'ID=MQ,', line)
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub('[0-9]+\tGT\t.\/.$', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##') and line.endswith('"'):
                        line = re.sub('"$', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        os.rename(temp_file, vcf)
        os.utime(vcf, times=(initial_file_time_stats.st_mtime, initial_file_time_stats.st_mtime))


class HTML_Summary():

    def __init__(self, runtime=None, vcf_to_df=None, reference=None, groupings_dict=None, raxml_version=None, all_vcf_boolen=None, args=None, removed_samples=None):

        htmlfile = open(f'{global_working_dir}/vSNP_step2_summary-{global_date_stamp}.html', 'at')
        
        #MAKE HTML FILE:
        print("<html>\n<head><style> table { font-family: arial, sans-serif; border-collapse: collapse; width: 40%; } td, th { border: 1px solid #dddddd; padding: 4px; text-align: left; font-size: 11px; } </style></head>\n<body style=\"font-size:12px;\">", file=htmlfile)

        print(f"<h2>Script ran using <u>{reference} </u> variables<br>", file=htmlfile)

        if args.metadata:
            print(f"<h4>Metadata:  {args.metadata}<br>", file=htmlfile)
        else:
            print(f"No metadata for describing samples in trees and tables<br>", file=htmlfile)
        if args.defining_snps:
            print(f"Defining SNPs:  {args.defining_snps}<br>", file=htmlfile)
        else:
            print(f"No defining SNPs files for grouping and filtering<br>", file=htmlfile)
        if args.gbk:
            for each in args.gbk:
                print(f"gbk:  {each}<br>", file=htmlfile)
        else:
            print(f"No gbk for annotation<br>", file=htmlfile)

        print(f"SNP calling thresholds:  REF: QUAL <u><{args.n_threshold}</u>, N: QUAL <u>{args.n_threshold}-{args.qual_threshold}</u>, ALT: QUAL <u>>{args.qual_threshold}</u>, Ambigious: <u>AC=1</u>, MQ: <u>>{args.mq_threshold}</u></h4>", file=htmlfile)

        print(f"<h4>{vcf_to_df.vcf_original_count} VCF files initial count<br>", file=htmlfile)
        print(f"{len(vcf_to_df.dataframes)} VCF files in this run<br>", file=htmlfile)
        print(f"{len(vcf_to_df.vcf_bad_list)} VCF files in this run were corrupt and therefore removed</h4>", file=htmlfile)
        
        if all_vcf_boolen:
            print("\n<h4>All_VCFs is available</h4>", file=htmlfile)

        #TIME
        print(f"Total run time: {runtime}: </h4>", file=htmlfile)

        # ERROR LIST
        if len(vcf_to_df.vcf_bad_list) < 1:
            print("<h2>No corrupt files found</h2>", file=htmlfile)
        else:
            print("\n<h2>Corrupt files removed</h2>", file=htmlfile)
            for i in vcf_to_df.vcf_bad_list:
                print(f"{i} <br>", file=htmlfile)
            print("<br>", file=htmlfile)

        #GROUPING TABLE
        group_vcfs_dict = defaultdict(list) #invert the key, values
        for group, dataframes in groupings_dict.items():
            for vcf in dataframes.keys():
                group_vcfs_dict[vcf].append(group)
        group_vcfs_dict = dict(sorted(group_vcfs_dict.items())) #sorts on key(vcf name)

        print(f'<h2>Groupings with {len(group_vcfs_dict):,} listed:</h2>', file=htmlfile)
        print("<table>", file=htmlfile)
        print("<tr align=\"left\"><th>Sample Name</th><tr>", file=htmlfile)

        for key, value in group_vcfs_dict.items():
            print("<tr>", file=htmlfile)
            print(f"<td>{key}</td>", end='\t', file=htmlfile)
            for group in value:
                print(f"<td>{group}</td>", end='\t', file=htmlfile)
            print("</tr>", file=htmlfile)
        print("</table>", file=htmlfile)

        # Removed from analysis
        if removed_samples:
            if len(removed_samples) < 1:
                print("<h2>No samples purposely removed from initial VCF file dataset</h2>", file=htmlfile)
            else:
                print('\n<h2>VCF files in "remove_from_analysis.xlsx" and not included from dataset.</h2>', file=htmlfile)
                for each in removed_samples:
                    print(f"{os.path.basename(each)} <br>", file=htmlfile)
                print("<br>", file=htmlfile)

        try:
            print("\n<h2>Program versions:</h2>", file=htmlfile)
            print(f'vSNP3: {__version__}', file=htmlfile)
            print(f'Python: {sys.version} <br>', file=htmlfile)
            program_list = ['Python', 'Bio', 'allel', 'numpy', 'pandas', 'scipy',]
            for name, module in sorted(sys.modules.items()): 
                if hasattr(module, '__version__') and name in program_list: 
                    print(f'{name}, {module.__version__} <br>', file=htmlfile)
            raxml_version = re.sub('This is ', '', raxml_version)
            raxml_version = re.sub(' released by.*', '', raxml_version)
            print(f'{raxml_version} <br>', file=htmlfile)
        except:
            pass

        print("</body>\n</html>", file=htmlfile)
        htmlfile.close()


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Store VCF files from vSNP step1 to step 2 directory.  VCF files must be stored by reference type.  Make a VCF file directory database that will build over time as samples are ran in step 1

    For example...
    <path/to/files>
        referenceA_dir
            step1_dir
                sample1_dir
                    <alignemnt_files>
                sample2_dir
                    <alignemnt_files>
            step2_dir
                vcf_source_dir
                    sample1.vcf
                    sample2.vcf
                comparison1_dir
                comparison2_dir
            
    <path/to/dependencies> (added path using vsnp3_path_adder.py)
        referenceA_dir
            defining_snps_for_referenceA.xlsx
            metadata_for_referenceA.xlsx
            FASTA/s for referenceA
            GBK/s for referenceA

    When running samples through step 1 and 2 of vSNP, or when running a routine analysis, set up dependencies using vsnp3_path_adder.py.  See vsnp3_path_adder.py -h for adding a reference type and for more information

    Usage:

    vsnp3_step2.py -a -d -t ASFV_Georgia_2007

    vsnp3_step2.py -wd <path/to/vcf_directory> -abs_pos chrom1:123456 -group test_group -m <path/to/metadata.xlsx>

    vsnp3_step2.py -wd ../vcf_source

    vsnp3_step2.py -a --remove

    ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-wd', '--wd', action='store', dest='wd', required=False, default='.', help='Optional: path to VCF files. By default .vcf in current working directory are used')
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', default=None, required=False, help='Optional: A valid reference_type name will be automatically found, but a valid reference_type name can be supplied.  See vsnp3_path_adder.py -s')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple gbk files can be specified with wildcard')
    parser.add_argument('-s', '--defining_snps', action='store', dest='defining_snps', default=None, required=False, help='Optional: Defining SNPs with positions to filter.  See template_define_filter.xlsx in vsnp dependency folder.  Recommended having this file in reference type folder')
    parser.add_argument('-m', '--metadata', action='store', dest='metadata', default=None, required=False, help='Optional: Two column Excel file, Column One: full VCF file name, Column Two: Updated name.  Recommended having this file in reference type folder')
    parser.add_argument('-remove_by_name', '--remove_by_name', action='store', dest='remove_by_name', required=False, help='Optional: Excel file containing samples to remove from analysis Column 1: to match sample name minus extension. No header allowed.   Recommended having this file in reference type folder')
    parser.add_argument('-n', '--no_filters', action='store_true', dest='no_filters', default=False, help='Optional: turn off filters')
    parser.add_argument('-w', '--qual_threshold', action='store', dest='qual_threshold', default=150, required=False, help='Optional: Minimum QUAL threshold for calling a SNP')
    parser.add_argument('-x', '--n_threshold', action='store', dest='n_threshold', default=50, required=False, help='Optional: Minimum N threshold.  SNPs between this and qual_threshold are reported as N')
    parser.add_argument('-y', '--mq_threshold', action='store', dest='mq_threshold', default=56, required=False, help='Optional: At least one position per group must have this minimum MQ threshold to be called.')
    parser.add_argument('-f', '--fix_vcfs', action='store_true', dest='fix_vcfs', help='Optional: Just fix VCF files and exit')
    parser.add_argument('-k', '--keep_ind_vcfs', action='store_true', dest='keep_ind_vcfs', default=False, help='Optional: Keep VCF files in current working directory when VCF files in current working director are used, VCF files are always saved and zipped in "vcf_starting_files.zip".')
    parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', required=False, help='Optional: create table with all isolates')
    parser.add_argument('-i', '--find_new_filters', action='store_true', dest='find_new_filters', help='Optional: find new positions to apply to the filter file.  Positions must be manually added to filter file.  They are not added by running this command.  Only text files are output showing position detail. Curant before adding filters')
    parser.add_argument('-abs_pos', '--abs_pos', action='store', dest='abs_pos', required=False, help='Optional: Make a group on defining SNP.  Must be supplied with --group option.  Format as chrom in VCF, chrom:10000.')
    parser.add_argument('-group', '--group', action='store', dest='group', required=False, help='Optional: Name a group on defining SNP.  Must be supplied with --abs_pos option')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling.  A pickle file will be kept for troubleshooting to be used directly in vsnp3_group_on_defining_snps.py.  This saves processing time')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    setup = Setup(debug=args.debug)
    global_date_stamp = setup.date_stamp
    global_working_dir = setup.cwd
    cwd_test = False
    if args.wd == '.':
        cwd_test = True
        vcf_list = glob.glob('*vcf')
    else:
        vcf_list = glob.glob(f'{args.wd}/*vcf')

    def zipit(src, dst):
        zf = zipfile.ZipFile("%s.zip" % (dst), "w", zipfile.ZIP_DEFLATED)
        abs_src = os.path.abspath(src)
        for dirname, subdirs, files in os.walk(src):
            for filename in files:
                absname = os.path.abspath(os.path.join(dirname, filename))
                arcname = absname[len(abs_src) + 1:]
                zf.write(absname, arcname)
        zf.close()
        shutil.rmtree(src)

    starting_files = f'{setup.cwd}/vcf_starting_files'
    os.makedirs(starting_files)
    for each_vcf in vcf_list:
        shutil.copy(each_vcf, starting_files)

    vcf_to_df = VCF_to_DF(vcf_list=vcf_list, debug=args.debug) #write_out=args.write_out
    if args.fix_vcfs:
        sys.exit(0)

    print(f'\nvcf_bad_list')
    for each in vcf_to_df.vcf_bad_list:
        print(f'\t {each}')

    if args.reference_type:
        ro = Ref_Options(args.reference_type)
    else:
        ro = Ref_Options(vcf_to_df.chrom)

    if args.abs_pos and not args.group:
        print('\n### -abs_pos must be used with -group option\n')
        sys.exit()
    if args.group and not args.abs_pos:
        print('\n### -group must be used with -abs_pos option\n')
        sys.exit()

    if ro.metadata and not args.metadata:
        args.metadata = ro.metadata
    if ro.defining_snps and not args.defining_snps:
        args.defining_snps = ro.defining_snps
    if ro.gbk and not args.gbk:
        args.gbk = ro.gbk
    if ro.remove and not args.remove_by_name:
        args.remove_by_name = ro.remove
    
    print(f'Before sample filter: {len(vcf_to_df.dataframes)}')
    if args.remove_by_name:
        remove_from_analysis = Remove_From_Analysis(working_directory=global_working_dir, excel_remove=args.remove_by_name, extension="vcf")
        for each in remove_from_analysis.remove_list:
            vcf_to_df.dataframes.pop(os.path.basename(each), None)
        remove_list = remove_from_analysis.remove_list
    else:
        remove_list = None
    print(f'After sample filter: {len(vcf_to_df.dataframes)}')

    if args.defining_snps:
        shutil.copy(args.defining_snps, starting_files) #package with starting files for the record
    zipit(starting_files, starting_files) # zip starting files directory

    group = Group(cwd=global_working_dir, metadata=args.metadata, defining_snps=args.defining_snps, excel_remove=args.remove_by_name, gbk_list=args.gbk, dataframes=vcf_to_df.dataframes, all_vcf=args.all_vcf, find_new_filters=args.find_new_filters, no_filters=args.no_filters, qual_threshold=int(args.qual_threshold), n_threshold=int(args.n_threshold), mq_threshold=int(args.mq_threshold), abs_pos=args.abs_pos, group=args.group, debug=args.debug)
    vcf_to_df.vcf_bad_list = vcf_to_df.vcf_bad_list + group.vcf_bad_list

    #by default the VCF files used are not deleted.  They are only deleted when using --remove option AND the files were ran from the current working directory.  ie the --wd option was not used.:
    if not args.keep_ind_vcfs and cwd_test:
        for each_vcf in vcf_list:
            try:
                os.remove(each_vcf)
            except FileNotFoundError:
                # if file was previously removed such as it was empty
                pass

    setup.print_time()
    HTML_Summary(runtime=setup.run_time, vcf_to_df=vcf_to_df, reference=ro.select_ref, groupings_dict=group.groupings_dict, raxml_version=group.raxml_version, all_vcf_boolen=args.all_vcf, args=args, removed_samples=remove_list) 

# Created 2021 by Tod Stuber