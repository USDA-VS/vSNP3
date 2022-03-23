#!/usr/bin/env python

__version__ = "3.0"

import os
import sys
import re
import pickle
import argparse
import textwrap
import numpy as np
import pandas as pd
from collections import defaultdict
from collections import Counter
from concurrent import futures
import multiprocessing
import time
from datetime import datetime
from cpuinfo import get_cpu_info

import warnings
warnings.filterwarnings('ignore')

from vsnp3_reference_options import Ref_Options
from vsnp3_fasta_to_snps_table import Tree
from vsnp3_fasta_to_snps_table import Tables
from vsnp3_annotation import Annotation


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

class Group():
    ''' 
    '''
    def __init__(self, cwd=None, metadata=None, excel_remove=None, gbk_list=None, defining_snps=None, dataframes=None, pickle_file=None, abs_pos=None, group=None, all_vcf=None, find_new_filters=None, no_filters=True, qual_threshold=150, n_threshold=50, mq_threshold=56, debug=False):

        self.qual_threshold = qual_threshold
        self.n_threshold = n_threshold
        self.mq_threshold = mq_threshold
        self.find_new_filters = find_new_filters
        self.vcf_bad_list=[]
        filter_all_list=None
        defining_snps_dict = None
        self.debug = debug

        if cwd == None:
            self.cwd = os.getcwd()
        else:
            self.cwd = cwd

        cpu_count = int(multiprocessing.cpu_count() / 1.2)
        
        if abs_pos and group:
            print(f'Dropping {defining_snps} for single grouping')
            defining_snps = None
            no_filters = True
        
        self.beginTime = datetime.now()
        self.startTime = datetime.now()
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
        print("\nSorting defining SNPs into groups...")
        all_vcf_name = None #declare variable as None/False as default
        if abs_pos and group:
            abs_pos = abs_pos
            group = group
            defining_snps_dict={}
            defining_snps_dict[abs_pos] = group
        elif defining_snps:
            defining_snps_df = pd.read_excel(defining_snps)  
            defining_snps_dict = defining_snps_df.iloc[:, 1:].head(n=1).to_dict(orient='records')[0] # drop first "all" column, just keep abs_pos and group, and make into dictionary of key=abs_pos, item=group
            if not no_filters:
                filter_all_list = defining_snps_df.iloc[:, 0].to_list()[1:]
                filter_all_list = [x for x in filter_all_list if str(x) != 'nan']
                filter_all_list = self.list_expansion(filter_all_list)
            filter_snps_df = pd.read_excel(defining_snps, header=1)
            group_filter_snps_dict = filter_snps_df.iloc[:, 1:].to_dict(orient='list')
            group_filter_snps_dict = {k:[elem for elem in v if elem is not np.nan] for k,v in group_filter_snps_dict.items()}
            self.group_filter_snps_dict = group_filter_snps_dict
            #get all_vcf column name for labeling group
            all_vcf_column = filter_snps_df.iloc[:, 0:1].to_dict(orient='list')
            all_vcf_name = [*all_vcf_column][0]

        if metadata:
            metadata_test = True
            metadata_df = pd.read_excel(metadata, header=None, index_col=0, usecols=[0, 1], names=['file_name', 'metadata'])
            #for names to match change to string types in case sample name are numbers/int
            metadata_df = metadata_df.reset_index()
            metadata_df['file_name'] = metadata_df.file_name.astype(str)
            metadata_df['metadata'] = metadata_df.metadata.astype(str)
            #fix metadata tags
            metadata_df['metadata'] = metadata_df.metadata.replace({'/':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\.':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\*':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\?':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\(':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\)':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\[':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'\]':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({' ':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'{':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'}':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'-_':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'_-':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'--':'_'}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'_$':''}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({'-$':''}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({"\'": ""}, regex=True)
            metadata_df['metadata'] = metadata_df.metadata.replace({',':''}, regex=True)
            
        else:
            metadata_test = False

        print(f'{bcolors.RED}\nDefining SNPs: {bcolors.ENDC}{bcolors.WHITE}{defining_snps}{bcolors.ENDC}')
        print(f'{bcolors.RED}Metadata: {bcolors.ENDC}{bcolors.WHITE}{metadata}{bcolors.ENDC}')
        print(f'{bcolors.RED}Remove From Analysis: {bcolors.ENDC}{bcolors.WHITE}{excel_remove}{bcolors.ENDC}')
        print(f'{bcolors.RED}gbks: {bcolors.ENDC}', end="")
        if gbk_list:
            for each in gbk_list:
                print(f'\t{bcolors.WHITE}{each}{bcolors.ENDC}')
        else:
            print(f'\t{bcolors.WHITE}No gbk{bcolors.ENDC}')
            gbk_list = None

        print(f'\nSorting defining SNPs  Selection Time: {datetime.now() - self.startTime}\n')

        if pickle_file:
            with open(pickle_file, 'rb') as handle:
                dataframes = pickle.load(handle)
        else:
            pass # dataframe was passed

        self.startTime = datetime.now()
        dataframe_essentials={}
        abs_pos_nt_dict={}
        annotation_dict={}
        map_quality_dict = defaultdict(list)
        print("Getting dataframe essential positions...")
        #### Change names with metadata #####
        dataframes_names_updated={} # collect all positions for second interation
        for sample, single_df in dataframes.items(): #just do once
            #work with basename if vcfs called from path
            sample = os.path.basename(sample)
            # single_df.loc[single_df['ALT'].str.len() > 1, 'ALT'] = 'N' # keep indel positions Ns, ie. ALT indels to N, REF indels handled in make_groupings function
            if metadata_test: #update name with metadata if provided
                try:
                    sample = metadata_df.loc[metadata_df['file_name'] == sample, 'metadata'].iloc[0]
                    dataframes_names_updated[sample] = single_df
                except IndexError:
                    try: #try stripping off vcf
                        sample = re.sub('.vcf$', '', sample)
                        sample = metadata_df.loc[metadata_df['file_name'] == sample, 'metadata'].iloc[0]
                        dataframes_names_updated[sample] = single_df
                    except IndexError:
                        try: #try stripping off _zc
                            sample = re.sub('_zc$', '', sample)
                            sample = metadata_df.loc[metadata_df['file_name'] == sample, 'metadata'].iloc[0]
                            dataframes_names_updated[sample] = single_df
                        except IndexError:
                            dataframes_names_updated[sample] = single_df
            else:
                dataframes_names_updated[sample] = single_df
            if not no_filters and filter_all_list:
                single_df = single_df[~single_df['abs_pos'].isin(filter_all_list)]
            # First iteration.  Find good SNPs for each VCF.  There must be at least one good SNP to include position in table
            try:
                single_df = single_df[(single_df['QUAL'] > self.qual_threshold) & (single_df['AC'] == 2) & (single_df['REF'].str.len() == 1) & (single_df['ALT'].str.len() == 1) & (single_df['MQ'] >= self.mq_threshold)]
            except AttributeError:
                print(f'\n### Error with sample {sample}\nSee VCF file and rerun\n')
                sys.exit(1)

            mq_dictionary = single_df[['abs_pos', 'MQ']].set_index('abs_pos').to_dict()['MQ']
            for abs_pos, MQ in mq_dictionary.items():
                map_quality_dict[abs_pos].append(MQ)

            if single_df.empty:
                self.vcf_bad_list.append(f'{sample}')
            else:
                dataframe_essentials[sample] = single_df
            sample_dict = dict(zip(single_df.abs_pos, single_df.ALT))
            abs_pos_nt_dict = {**abs_pos_nt_dict, **sample_dict}
        self.dataframe_essentials = dataframe_essentials
        self.dataframes_names_updated = dataframes_names_updated
        samples_with_dataframes_set = set(dataframe_essentials.keys())
        dataframes={} #empty memory
        mq_averages={}
        for abs_pos, mq_list in map_quality_dict.items():
            mq_averages[abs_pos] = np.mean(mq_list)
        self.average_mq_df =  pd.DataFrame.from_dict(mq_averages, orient='index')
        self.average_mq_df = self.average_mq_df.reset_index()
        self.average_mq_df.columns = ['abs_pos', 'MQ']
        if gbk_list:
            annotation = Annotation(gbk_list=gbk_list)
            for abs_pos, snp_nt in abs_pos_nt_dict.items():
                annotation.run(abs_pos, snp_nt)
                if annotation.reference_base_code == 'n/a':
                    annotation_dict[abs_pos] = 'position not annotated'
                else:
                    annotation_dict[abs_pos] = f'{annotation.reference_base_code}->{annotation.snp_base_code}, {annotation.gene}:{annotation.ref_aa}{annotation.aa_residue_pos}{annotation.snp_aa}, {annotation.product}, {annotation.mutation_type}'
        self.annotation_df = pd.DataFrame(annotation_dict.items(), columns=['abs_pos', 'annotation'])
        print(f'\n\tGetting dataframe essentials  Selection Time: {datetime.now() - self.startTime}\n')

        if defining_snps_dict:
            self.startTime = datetime.now()
            groupings_dict = {} # gather groups (key) and sample names (values as list)
            defining_snps_list=[]
            for abs_pos, group in defining_snps_dict.items():
                defining_snps_list.append(abs_pos)
            # if  debug:
            for abs_pos in defining_snps_list:
                group_found, sample_dict = self.group_selection(abs_pos) #sample list is list of sample dataframes with defining snp
                if group_found:
                    group = defining_snps_dict[abs_pos]
                    if not no_filters: #don't apply filters when option called
                        sample_dict_group_filter={}
                        for sample, each_df in sample_dict.items():
                            each_df = each_df[~each_df['abs_pos'].isin(self.group_filter_snps_dict[group])] #by group remove positions to filter
                            sample_dict_group_filter[sample] = each_df
                        sample_dict = sample_dict_group_filter
                    groupings_dict[group] = sample_dict #defining_snps_dict[abs_pos] provides group
                    
            # else:
            #     print(f'Grouping: Pool processing with {cpu_count} cpus...')
            #     with futures.ThreadPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## thread call is more efficient
            #         for group_found, sample_dict in pool.map(self.group_selection, defining_snps_list):
            #             if group_found:
            #                 group = defining_snps_dict[abs_pos]
            #                 if not no_filters:
            #                     sample_dict_group_filter={}
            #                     for sample, each_df in sample_dict.items():
            #                         each_df = each_df[~each_df['abs_pos'].isin(self.group_filter_snps_dict[group])] #by group remove positions to filter
            #                         sample_dict_group_filter[sample] = each_df
            #                     sample_dict = sample_dict_group_filter
            #                 groupings_dict[group] = sample_dict #defining_snps_dict[abs_pos] provides group

            combined_lists=[]
            for group, samples in groupings_dict.items():
                combined_lists = combined_lists + list(samples.keys())
            samples_with_group_set = set(combined_lists)
            samples_without_group_set = samples_with_dataframes_set - samples_with_group_set
            print(f'\n\tGroup Selection Time: {datetime.now() - self.startTime}\n')
        else:
            samples_without_group_set = samples_with_dataframes_set
            groupings_dict = {}

        if all_vcf or not defining_snps_dict:
            if all_vcf_name:
                groupings_dict[all_vcf_name] = self.dataframe_essentials
            else:
                groupings_dict['all_vcf'] = self.dataframe_essentials

        print(f'All relevant positions by group')
        self.startTime = datetime.now()
        #Get all position in each group
        ambigious_lookup={}
        ambigious_lookup['AG'] = 'R'
        ambigious_lookup['CT'] = 'Y'
        ambigious_lookup['GC'] = 'S'
        ambigious_lookup['AT'] = 'W'
        ambigious_lookup['GT'] = 'K'
        ambigious_lookup['AC'] = 'M'
        ambigious_lookup['GA'] = 'R'
        ambigious_lookup['TC'] = 'Y'
        ambigious_lookup['CG'] = 'S'
        ambigious_lookup['TA'] = 'W'
        ambigious_lookup['TG'] = 'K'
        ambigious_lookup['CA'] = 'M'
        self.ambigious_lookup = ambigious_lookup

        finished_groupings_dict={} #normalized/sequence aligned
        groupings_dict_list=[]
        for group, sample_dict in groupings_dict.items():
            groupings_dict_list.append((group, sample_dict))

        # if debug:
        for group_sample_dict in groupings_dict_list:
            parsmony_sample_dict, group = self.make_groupings(group_sample_dict)
            finished_groupings_dict[group] = parsmony_sample_dict
        finished_groupings_dict = {i:j for i,j in finished_groupings_dict.items() if j != {}} #remove items if vaule is empty
        # else:
        #     with futures.ThreadPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## thread is diffently better but putting this in futures is slightly slower.
        #         for parsmony_sample_dict, group in pool.map(self.make_groupings, groupings_dict_list):
        #             finished_groupings_dict[group] = parsmony_sample_dict
        
        #Find positions that need to be filtered
        if self.find_new_filters:
            for group_dict_of_df in groupings_dict_list:
                if len(group_dict_of_df[1]) > 3:
                    for df in group_sample_dict:
                        postion_list = open(f'{group_dict_of_df[0]}_postion_list.txt', 'w')
                        postion_detail_list = open(f'{group_dict_of_df[0]}_postion_detail_list.txt', 'w')
                        if not no_filters:
                            print('New positions to filter found after current filter positions applied but before noninformative SNP are removed', file=postion_detail_list)
                        else:
                            print('No previous filtering applied', file=postion_detail_list)
                        print('dd.QUAL.mean() < 700 and dd.QUAL.max() < 1300 or dd.MQ.mean() < 56', file=postion_list)
                        print('dd.QUAL.mean() < 700 and dd.QUAL.max() < 1300 or dd.MQ.mean() < 56', file=postion_detail_list)
                        cc = pd.concat(group_dict_of_df[1].values(), ignore_index=True)
                        cc['abs_pos'] = cc['CHROM'] + ':' + cc['POS'].astype(str)
                        ll = set(cc['abs_pos'].to_list())
                        for vv in ll:
                            dd = cc[cc['abs_pos'] == vv]
                            if len(dd) > 3:
                                if dd.QUAL.mean() < 700 and dd.QUAL.max() < 1300 or dd.MQ.mean() < 56:
                                    print(vv, file=postion_list)
                                    print(f'{vv} Average QUAL: {dd.QUAL.mean():0.2f}, Max QUAL: {dd.QUAL.max():0.2f}, Average MQ: {dd.MQ.mean():0.2f}', file=postion_detail_list)
                        postion_list.close()
                        postion_detail_list.close()
        print(f'\n\tAll relevant positions by group {datetime.now() - self.startTime}\n')

        print(f'FASTAs out and RAxML trees')
        self.startTime = datetime.now()
        group_fasta_dict={}
        group_dataframe_dict={}
        remove_list=[]
        for group, sample_dict in finished_groupings_dict.items():
            os.mkdir(group)
            fasta = f'{group}/{group}-{self.st}.fasta'
            self.dict_to_fasta(sample_dict, fasta)
            num_lines=0
            with open(fasta) as opened_file:
                for line in opened_file:
                    num_lines+=1
                    string = line.split()[0]
                    if string.startswith('>'):
                        pass
                    else:
                        read_length = len(string)
            if num_lines < 7 or read_length <= 3 : # ie 4 or more FASTAs and sequence length > 3 required
                with open(f'{group}/TOO_FEW_SAMPLES_OR_SHORT_SEQUENCE_TO_BUILD_TREE', 'w') as message_out:
                    print(f'check sample numbers or sequence lengths', file=message_out)
                remove_list.append(group)
            else:
                group_fasta_dict[group] = fasta
                group_df = self.dict_to_dataframe(sample_dict, group)
                group_dataframe_dict[group] = group_df
        
        self.group_fasta_dict = group_fasta_dict
        self.group_dataframe_dict = group_dataframe_dict
        finished_groupings_list = [*finished_groupings_dict]
        working_group_list = [x for x in finished_groupings_list if x not in remove_list]

        if debug:
            for group in working_group_list:
                tree = self.raxml_table_build(group)
        else:
            with futures.ThreadPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## thread works best for raxml calls
                for tree in pool.map(self.raxml_table_build, working_group_list):
                    pass

        print(f'\n\tFASTAs out and RAxML trees {datetime.now() - self.startTime}\n')
        # print(f'\n\nTotal Time: {datetime.now() - self.beginTime}\n')

        #Add back those that where a group was not found
        for sample in samples_without_group_set:
            groupings_dict = {**groupings_dict, 'Group Not Found': {sample: None}}
        self.groupings_dict = groupings_dict # will be passed to html summary

    def group_selection(self, abs_pos):
        sample_dict={}
        group_found = False
        for sample, single_df in self.dataframe_essentials.items():
            if any(single_df['abs_pos'] == abs_pos):
                group_found = True
                sample_dict[sample] = self.dataframe_essentials[sample]
            elif "!" in abs_pos:
                abs_pos_inverted = re.sub('!', '', abs_pos)
                if not any(single_df['abs_pos'] == abs_pos_inverted):
                    group_found = True
                    sample_dict[sample] = self.dataframe_essentials[sample]
        return group_found, sample_dict

    def list_expansion(self, target_list):
        expanded_list=[]
        for list_entry in target_list:
            list_entry = str(list_entry)
            if "-" not in list_entry.split(":")[-1]:
                expanded_list.append(list_entry)
            elif "-" in list_entry:
                try:
                    chrom, sequence_range = list_entry.split(":")
                except ValueError as e:
                    raise type(e)(str(e) + f' \n#### error in Defining SNPs/Filter worksheet\n#### see value "{list_entry}"').with_traceback(sys.exc_info()[2])
                list_entry = sequence_range.split("-")
                for position in range(int(list_entry[0].replace(',', '')), int(list_entry[1].replace(',', '')) + 1):
                    expanded_list.append(chrom + ":" + str(position))
        return expanded_list

    def dict_to_fasta(self, sample_dict, fasta): #sample_dict = [file name]:snp dataframe
        with open(fasta, 'w') as write_out:
            for name, df in sample_dict.items():
                print(f'>{name}', file=write_out)
                try:
                    df = df.sort_values(by=['abs_pos']) # sorting ensures positions are aligned
                    # print(f'{name}: {len("".join(df["ALT"].to_list()))}') # use to troubleshoot if FASTAs are not aligning to the same length
                    print("".join(df['ALT'].to_list()), file=write_out)
                except KeyError:
                    pass

    def dict_to_dataframe(self, sample_dict, group): # dict_to_fasta will not have abs_pos this dataframe will retain abs_pos
        try:
            dict_of_dfs = { sample: df.set_index('abs_pos') for sample, df in sample_dict.items()} #set ALT column to sample name and merge
            group_df = pd.concat(dict_of_dfs, axis=1)
            group_df.columns = group_df.columns.droplevel(-1)
            return group_df.T
        except pd.errors.InvalidIndexError as e:
            if self.debug:
                print(f'\n\t#####\n\t##### {e}, Group {group} \n\t#####\n')

    def make_groupings(self, group_sample_dict):
        df_list=[]
        group, sample_dict = group_sample_dict # sample_dict is from dataframe_essentials, ie good SNPs.
        for sample, df in sample_dict.items():
            df_list.append(df[['abs_pos', 'REF']])  
        df = pd.concat(df_list, ignore_index=True)
        df_ref = df.drop_duplicates(subset='abs_pos', keep="first") # df_ref will be all abs_pos per group.  completion of first iteration
        norm_sample_dict={}
        position_list=[]
        position_list_parse_test = []
        for sample in sample_dict.keys(): # update grouping sample_dict with normalized dataframes
            # dataframes_names_updated contains all SNPs
            # dataframe_essentials contains good SNP positions
            sample_df = self.dataframes_names_updated[sample]
            sample_df = sample_df[sample_df['abs_pos'].isin(df_ref['abs_pos'])] # this will normalize positions
            #https://stackoverflow.com/questions/27673231/why-should-i-make-a-copy-of-a-data-frame-in-pandas
            sample_df_parse_test = sample_df.copy() # to use below            
            sample_df.loc[sample_df['ALT'].str.len() > 1, 'ALT'] = 'N'
            sample_df.loc[sample_df['ALT'] == 'N', 'AC'] = 2 # allow the above line to pass ambigious if needed        
            try: #change AC=1 to ambigious
                for index, row in sample_df.loc[sample_df['AC'] == 1].iterrows():
                    sample_df.at[index, 'ALT'] = self.ambigious_lookup[row['REF'] + row['ALT']]
            except (KeyError, TypeError) as e:
                if self.debug:
                    print(f'\n\t#####\n\t##### {e}, Sample: {sample}\n\t#####\n')
            #change alt to N if QUAL 50 - 150
            sample_df.loc[(sample_df['QUAL'] >= self.n_threshold) & (sample_df['QUAL'] < self.qual_threshold), 'ALT'] = 'N' # this will overwrite ambigious calls
            # < 50 will default to REF... change ALT to REF
            try:
                sample_df.loc[sample_df['REF'].str.len() > 1, 'REF'] = 'N' #if REF call is indel change to N to maintain equal sequence length for all samples
                sample_df.loc[sample_df['QUAL'] < self.n_threshold, 'ALT'] = sample_df['REF']
            except (ValueError) as e:
                if self.debug:
                    print(print(f'\n\t#####\n\t##### {e}, Sample: {sample}\n\t#####\n'))

            sample_df = sample_df[['abs_pos', 'ALT']] # no longer need other columns
            sample_df = sample_df.replace(np.nan, '-', regex=True) # change zero coverage to -
            df_merged = sample_df.merge(df_ref, left_on='abs_pos', right_on='abs_pos', how='outer') # finish normalizing if df doesn't include all position in group
            df_merged['ALT'] = np.where(df_merged['ALT'].notna(), df_merged['ALT'], df_merged['REF']) # merge REF from df_ref to ALT column
            df_merged['ALT'] = df_merged['ALT'].fillna('-')
            
            # df_merged[~sample_df.isin(df_ref)].dropna() #make sure only those position in df_ref are being used | not sure this is needed
            position_list.extend(list(df_merged['abs_pos'] + list(df_merged['ALT']))) #make position unique on ALT call for parsimony selection
            norm_sample_dict[sample] = df_merged[['abs_pos', 'ALT']]
            
            ####
            #Do not include low quality calls when determining if position is parsimonious.  Only include calls with QUAL >= qual_threshold [150]
            try:
                sample_df_parse_test.loc[sample_df_parse_test['REF'].str.len() > 1, 'REF'] = 'N' #if REF call is indel change to N to maintain equal sequence length for all samples
                sample_df_parse_test.loc[sample_df_parse_test['QUAL'] < self.qual_threshold, 'ALT'] = sample_df_parse_test['ALT'] # change n_threshold from above to qual threshold skipping the Ns when determining if position is parsimonious AND changed to 'ALT'.  So, if there are just a few low quality represented the SNP position will be seen as parisomonious uninformative and removed.
            except (ValueError) as e:
                if self.debug:
                    print(print(f'\n\t#####\n\t##### {e}, Sample: {sample}\n\t#####\n'))
            sample_df_parse_test = sample_df_parse_test[['abs_pos', 'ALT']] # no longer need other columns
            sample_df_parse_test = sample_df_parse_test.replace(np.nan, '-', regex=True) # change zero coverage to -
            df_merged_parse_test = sample_df_parse_test.merge(df_ref, left_on='abs_pos', right_on='abs_pos', how='outer') # finish normalizing if df doesn't include all position in group
            df_merged_parse_test['ALT'] = np.where(df_merged_parse_test['ALT'].notna(), df_merged_parse_test['ALT'], df_merged_parse_test['REF']) # merge REF from df_ref to ALT column
            df_merged_parse_test['ALT'] = df_merged_parse_test['ALT'].fillna('-')
            position_list_parse_test.extend(list(df_merged_parse_test['abs_pos'] + list(df_merged_parse_test['ALT']))) #make position unique on ALT call for parsimony selection
            ####

        #find parsimonies uninformative positions
        counter = Counter(position_list_parse_test)
        most_common = counter.most_common()
        most_common = dict(most_common)
        parsimony_positions=[]
        for pos, count in most_common.items():
            if count == len(norm_sample_dict):
                parsimony_positions.append(pos[:-1]) #drop the ALT
        parsmony_sample_dict={}
        for sample, df_norm in norm_sample_dict.items():
            df_norm = df_norm[~df_norm['abs_pos'].isin(parsimony_positions)]
            if df_norm.empty:
                parsmony_sample_dict[sample] = pd.DataFrame()
            else:
                df_norm= self.sort_df(df_norm)
                parsmony_sample_dict[sample] = df_norm

        #Add root
        df_root = df_ref[df_ref['abs_pos'].isin(df_norm['abs_pos'])]
        df_root = df_norm.merge(df_root, left_on='abs_pos', right_on='abs_pos') #df with columns abs_pos, ALT, REF, in correct position order
        df_root = df_root[['abs_pos', 'REF']]
        df_root = df_root.rename(columns={"REF": "ALT"}) #fake the column name as ALT so it is as samples when dict_to_fasta function is called on the dataframes
        if df_root.empty:
            pass
        else:
            df_root= self.sort_df(df_root)
            parsmony_sample_dict['root'] = df_root
        return parsmony_sample_dict, group
    
    def sort_df(self, df):
        df[['chrom', 'pos']] = df['abs_pos'].str.split(':', expand=True)
        df['pos'] = df['pos'].astype(int)
        df = df.sort_values('pos').reset_index(drop=True)
        df = df[['abs_pos', 'ALT']]
        return df

    def raxml_table_build(self, group):
        tree = Tree(fasta_alignments=self.group_fasta_dict[group], write_path=f'{self.cwd}/{group}', tree_name=group)
        tables = Tables(df_alignments=self.group_dataframe_dict[group], tree=tree.newick, gbk=self.annotation_df, mq=self.average_mq_df, write_path=f'{self.cwd}/{group}', table_name=group, debug=False)
        tables.build_tables()
        self.raxml_version = tree.raxml_version



if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Usage:
    vsnp3_group_on_defining_snps.py -p dictionary_of_dataframes.pickle -abs_pos NC_002945.4:1295549 -group test_group -m path/to/metadata.xlsx -b /path/to/*gbk

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-p', '--pickle_file', action='store', dest='pickle_file', required=False, help='Pickle file: dictionary of dataframes')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-m', '--metadata', action='store', dest='metadata', default=None, required=False, help='Explicit metadata file.  Two column Excel file --> Column One: VCF file name, Column Two: Updated name')
    parser.add_argument('-s', '--defining_snps', action='store', dest='defining_snps', required=False, help='Defining SNPs with positions to filter.  See template_define_filter.xlsx')
    parser.add_argument('-abs_pos', '--abs_pos', action='store', dest='abs_pos', required=False, help='Must be supplied with --group option.  Format as chrom in VCF, likely chrom:10000... NC_002945.4:2138896.  Run: `vsnp3_step2.py --wd ../original -da` to obtain pickle for entire set, isolate pickle file and run `vsnp3_group_on_defining_snps.py -p dictionary_of_dataframes.pickle -a NC_002945.4:1295549`')
    parser.add_argument('-group', '--group', action='store', dest='group', required=False, help='Must be supplied with --abs_pos option')
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', required=False, help='Reference type group/directory with dependencies')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    if args.reference_type:
        ro = Ref_Options(args.reference_type)
        if ro.metadata and not args.metadata:
            args.metadata = ro.metadata
        if ro.excel and not args.defining_snps:
            args.defining_snps = ro.excel
        if ro.gbk and not args.gbk:
            args.gbk = ro.gbk

    group = Group(pickle_file=args.pickle_file, metadata=args.metadata, gbk_list=args.gbk, defining_snps=args.defining_snps, abs_pos=args.abs_pos, group=args.group, debug=args.debug)

# Created 2021 by Tod Stuber