#!/usr/bin/env python

from vsnp3_version import __version__

import os
import sys
import bisect
import re
import unicodedata
import pickle
import locale
import argparse
import textwrap
import numpy as np
import pandas as pd
from collections import defaultdict
from collections import Counter
from concurrent import futures
import multiprocessing
from multiprocessing import Process, Queue
from queue import Empty
import time
import traceback
from datetime import datetime

import warnings
# Targeted, not blanket.  The previous filterwarnings('ignore') hid real defects:
# a SettingWithCopyWarning on a write-to-a-slice in make_groupings, and pandas
# deprecations for positional Series access that is removed in pandas 3.0.  Only
# the known-noisy messages are suppressed so anything new is visible.
for _msg in (r'.*invalid value encountered.*',
             r'.*divide by zero encountered.*',
             r'.*DataFrame is highly fragmented.*',
             r'.*Passing a BlockManager.*'):
    warnings.filterwarnings('ignore', message=_msg)
warnings.filterwarnings('ignore', category=DeprecationWarning, module='openpyxl')

from vsnp3_version import (AC_HETEROZYGOUS, AC_HOMOZYGOUS, MQ_THRESHOLD,
                             N_THRESHOLD, QUAL_THRESHOLD, mixed_signature)
from vsnp3_reference_options import Ref_Options
from vsnp3_file_setup import worker_count
from vsnp3_fasta_to_snps_table import Tree
from vsnp3_fasta_to_snps_table import Tables
from vsnp3_annotation import Annotation
from vsnp3_html_tree import html_tree

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


# Define wrapper at module level
def wrapper(queue, func, *args, **kwargs):
    # The outcome is always reported, success or failure.  Before, an exception
    # in func ended the child without putting anything on the queue, and the
    # parent's unconditional queue.get() then waited forever: an html_tree
    # crash on one group hung the entire step 2 run with no message.
    try:
        result = func(*args, **kwargs)
    except BaseException as e:                                  # noqa: BLE001
        queue.put(('error', f'{type(e).__name__}: {e}\n{traceback.format_exc()}'))
        return
    queue.put(('ok', result))

def run_with_timeout(func, timeout_seconds, *args, **kwargs):
    name = getattr(func, '__name__', str(func))
    queue = Queue()
    process = Process(
        target=wrapper,
        args=(queue, func) + args,
        kwargs=kwargs
    )
    process.start()
    deadline = time.monotonic() + timeout_seconds
    message = None
    while message is None:
        try:
            message = queue.get(timeout=1)
        except Empty:
            if not process.is_alive():
                # Exited: its message is on the queue, or it was killed before
                # it could report (an OOM kill, for instance).
                try:
                    message = queue.get(timeout=5)
                except Empty:
                    process.join()
                    raise RuntimeError(f'{name} exited with status {process.exitcode} '
                                       f'without returning a result')
            elif time.monotonic() > deadline:
                process.terminate()
                process.join()
                raise TimeoutError("Function took longer than {} seconds".format(timeout_seconds))
    process.join(timeout=30)
    status, payload = message
    if status == 'error':
        raise RuntimeError(f'{name} failed: {payload}')
    return payload

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

def clean_unicode(text):
    """Clean unicode characters and normalize text for safe file operations"""
    if pd.isna(text):
        return text
    # Replace accented characters with their ASCII equivalents, drop the rest.
    text = unicodedata.normalize('NFD', str(text))
    return text.encode('ascii', errors='ignore').decode('ascii')


class DefiningSnpsError(ValueError):
    '''A mistake in the defining-SNP worksheet that the run must not continue with.'''


# Characters a Newick taxon name cannot carry unquoted.  RAxML exits on them,
# so they failed the whole group.  Metadata names are cleaned in read_metadata;
# this is for names derived from the VCF file name.
NEWICK_UNSAFE = re.compile(r"[\s:,;()\[\]']")


class Group():
    ''' 
    '''
    def __init__(self, cwd=None, metadata=None, excel_remove=None, gbk_list=None, defining_snps=None, dataframes=None, pickle_file=None, abs_pos=None, group=None, all_vcf=None, find_new_filters=None, no_filters=True, qual_threshold=QUAL_THRESHOLD, n_threshold=N_THRESHOLD, mq_threshold=MQ_THRESHOLD, show_groups=False, hash_groups=None, html_tree=False, dp=False, filter_density=False, density_threshold=3, density_window=20, debug=False, sample_stats=None):

        self.qual_threshold = qual_threshold
        self.n_threshold = n_threshold
        self.mq_threshold = mq_threshold
        self.find_new_filters = find_new_filters
        self.vcf_bad_list=[]
        # Samples analysed with no SNP passing the thresholds: name -> detail.
        # These used to go into vcf_bad_list -- reported as "corrupt" by the
        # step 2 summary, which they are not -- and were dropped from every
        # group, table and tree with nothing else said.
        self.no_snp_samples = {}
        # Samples that have passing SNPs but whose confident mixed calls
        # outnumber them several times over (vsnp3_version.mixed_signature):
        # name -> counts.  The summary highlights them in the Groupings table.
        self.mixed_samples = {}
        # Cells taken from inside a sample's multi-base record (make_groupings):
        # sample name -> number of table cells, summed over groups.
        self.embedded_cells = defaultdict(int)
        # group -> positions that left the group's table because a call read from
        # a multi-base record made every sample in the group carry the same base.
        self.embedded_now_uninformative = defaultdict(set)
        self._embedded_cache = {}
        # sample name -> VCF basename, for every sample that was analysed
        self.sample_sources = {}
        # parse statistics per VCF basename, from vsnp3_step2.py
        self.sample_stats = sample_stats or {}
        self.metadata_problem = None
        self.metadata_blank_rows = 0        # worksheet rows ignored for a blank cell
        self.fallback_renames = {}          # VCF basename -> (name before, name after)
        self.undetermined_in_group = defaultdict(list)   # group -> samples left out: no determined base
        self.rewrite_errors = {}            # (group, sample) -> error in the call rewrite
        self.worksheet_notes = []           # defining-SNP worksheet entries ignored or merged
        self.all_vcf_name = None
        self.html_tree_failures = {}
        self.phase_times = []
        self._phase_t0 = time.monotonic()
        # Assigned per group below.  Initialized here because vsnp3_step2.py reads
        # raxml_version at the very end of the run, and when no group produced a tree
        # the attribute never existed -- so a complete run died with AttributeError
        # after all the work was done.
        self.raxml_version = None
        self.group_failures = {}
        filter_all_list=None
        defining_snps_dict = None
        self.show_groups = show_groups
        self.html_tree = html_tree
        self.dp = dp
        self.debug = debug
        self.filter_density = filter_density
        self.density_threshold = density_threshold
        self.density_window = density_window
        self.density_filtered_positions = []  # Track positions filtered by density
        metadata_test = False
        
        self.sample_coverage_dict = {}

        if cwd == None:
            self.cwd = os.getcwd()
        else:
            self.cwd = cwd

        cpu_count = worker_count()     # honours VSNP3_MAX_CPUS

        if abs_pos and group:
            print('Dropping {} for single grouping'.format(defining_snps))
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
            if hash_groups:
                # Read the Excel file
                defining_snps_df = pd.read_excel(defining_snps, header=None)

                # Check the first row for "#" signs and remove them
                if not defining_snps_df.empty:
                    first_row = defining_snps_df.iloc[0]
                    cleaned_first_row = first_row.apply(lambda x: str(x).replace("#", "") if pd.notna(x) else x)
                    defining_snps_df.iloc[0] = cleaned_first_row

                    # Set the cleaned first row as column names and remove it from the data
                    defining_snps_df.columns = defining_snps_df.iloc[0]
                    defining_snps_df = defining_snps_df.drop(defining_snps_df.index[0])

                    # Reset the index
                    defining_snps_df = defining_snps_df.reset_index(drop=True)

                    # Remove empty columns
                    defining_snps_df = defining_snps_df.dropna(axis=1, how='all')

                    # Remove columns where all values are empty strings
                    defining_snps_df = defining_snps_df.loc[:, (defining_snps_df != '').any()]

                # Now defining_snps_df contains the data with "#" signs removed from the column names and empty columns removed
            else:
                defining_snps_df = pd.read_excel(defining_snps)  

            try:
                defining_snps_dict = defining_snps_df.iloc[:, 1:].head(n=1).to_dict(orient='records')[0] # drop first "all" column, just keep abs_pos and group, and make into dictionary of key=abs_pos, item=group
            except IndexError:
                defining_snps_dict = defining_snps_df.iloc[:, 0:].head(n=1).to_dict(orient='records')[0]
            if not no_filters:
                filter_all_list = defining_snps_df.iloc[:, 0].to_list()[1:]
                filter_all_list = [x for x in filter_all_list if str(x) != 'nan']
                filter_all_list = self.list_expansion(filter_all_list)
            filter_snps_df = pd.read_excel(defining_snps, header=1)
            group_filter_snps_dict = filter_snps_df.iloc[:, 1:].to_dict(orient='list')
            group_filter_snps_dict = {k:[elem for elem in v if not pd.isna(elem)] for k,v in group_filter_snps_dict.items()}
            group_filter_snps_dict = self._check_defining_snps(
                defining_snps, defining_snps_df, defining_snps_dict, group_filter_snps_dict,
                filter_all_list if not no_filters else [])
            self.group_filter_snps_dict = group_filter_snps_dict
            #get all_vcf column name for labeling group
            all_vcf_column = filter_snps_df.iloc[:, 0:1].to_dict(orient='list')
            all_vcf_name = list(all_vcf_column.keys())[0]

        # Metadata is read whether or not a defining-SNP worksheet is in use.
        # It used to be read inside the `elif defining_snps` branch above, so a
        # reference type without a metadata file -- and any -abs_pos/-group run
        # -- left metadata_df unassigned and the name resolution below died with
        # UnboundLocalError.  A metadata file is optional: without one the VCF
        # file names are the sample names, and the run goes on.  (The phrase
        # "kapurlab: metadata is optional" is the marker vsnp_gui's
        # deploy/vsnp3-patches/apply.sh checks before inserting its own copy of
        # this fix into an installed vsnp3; leave it in place.)
        metadata_df = None
        metadata_test = False
        if metadata:
            try:
                metadata_df = self.read_metadata(metadata)
                metadata_test = True
                self.metadata_blank_rows = int(metadata_df.attrs.get('blank_rows', 0))
                if self.metadata_blank_rows:
                    print(f'\n### Metadata: {self.metadata_blank_rows} row(s) with a blank cell '
                          f'ignored; the VCF file name is used for those samples.\n')
            except Exception as e:                              # noqa: BLE001
                self.metadata_problem = (
                    f'{os.path.basename(str(metadata))} could not be read '
                    f'({type(e).__name__}: {e})')
                print(f'\n### Metadata not used: {self.metadata_problem}. '
                      f'VCF file names are used as sample names.\n')

        print('{}Defining SNPs: {}{}{}'.format(bcolors.RED, bcolors.ENDC, bcolors.WHITE, defining_snps))
        print('{}{}'.format(bcolors.ENDC, ''))
        print('{}Metadata: {}{}{}'.format(bcolors.RED, bcolors.ENDC, bcolors.WHITE, metadata))
        print('{}{}'.format(bcolors.ENDC, ''))
        print('{}Remove From Analysis: {}{}{}'.format(bcolors.RED, bcolors.ENDC, bcolors.WHITE, excel_remove))
        print('{}{}'.format(bcolors.ENDC, ''))
        print('{}gbks: {}'.format(bcolors.RED, bcolors.ENDC), end="")
        if gbk_list:
            for each in gbk_list:
                print('\t{}{}{}'.format(bcolors.WHITE, each, bcolors.ENDC))
        else:
            print('\t{}No gbk{}'.format(bcolors.WHITE, bcolors.ENDC))
            gbk_list = None

        # Log density filtering status
        if self.filter_density:
            print('{}Density filtering: {}{}ENABLED - filtering SNPs when >={} found within {} bp window{}'.format(
                bcolors.RED, bcolors.ENDC, bcolors.WHITE, self.density_threshold, self.density_window, bcolors.ENDC))
        else:
            print('{}Density filtering: {}{}DISABLED{}'.format(
                bcolors.RED, bcolors.ENDC, bcolors.WHITE, bcolors.ENDC))

        print('\nSorting defining SNPs  Selection Time: {}\n'.format(datetime.now() - self.startTime))
        self._phase('Sort defining SNPs')

        if pickle_file:
            with open(pickle_file, 'rb') as handle:
                dataframes = pickle.load(handle)
        else:
            pass # dataframe was passed

        self.startTime = datetime.now()
        dataframe_essentials={}
        # Every ALT observed at each position, across all samples.  This used to be
        # a plain dict overwritten per sample, so the single ALT that reached the
        # annotator was whichever VCF happened to be processed last -- which made
        # the published amino acid change depend on file order and differ between
        # runs on identical data.
        abs_pos_nt_counts = defaultdict(Counter)
        annotation_dict={}
        map_quality_dict = defaultdict(list)
        print("Getting dataframe essential positions...")
        #### Change names with metadata #####
        dataframes_names_updated={} # collect all positions for second interation

        # Names are resolved for every input first, in their own pass, so a
        # collision is settled before any sample is processed.  Done inline
        # before, the loser was whichever file the dict happened to overwrite
        # last, and both dataframes_names_updated and dataframe_essentials could
        # end up describing different files under one name.
        resolved = {}
        matched_key = {}
        name_sources = defaultdict(list)
        # A file whose name matches several rows of the worksheet.  The lookup
        # takes the first row, so the others are discarded -- and because the file
        # still ends up with exactly one name, no count changes and nothing else
        # would ever mention it.  Only reported when the rows disagree; identical
        # duplicate rows pick the same name either way.
        self.metadata_ambiguous = {}
        for source in dataframes:
            base = os.path.basename(source)
            name, key, hits = self.resolve_sample_name_detail(
                base, metadata_df, metadata_test)
            if key is None:
                # A file-name-derived name goes to RAxML as a Newick taxon; a
                # space, quote, bracket, colon, comma or semicolon in the file
                # name failed RAxML, and with it the whole group.
                safe = self.safe_fallback_name(name)
                if safe != name:
                    self.fallback_renames[base] = (name, safe)
                    name = safe
            resolved[source] = name
            matched_key[source] = key
            if hits > 1:
                targets = list(dict.fromkeys(
                    metadata_df.loc[metadata_df['file_name'] == key,
                                    'metadata'].tolist()))
                if len(targets) > 1:
                    self.metadata_ambiguous[base] = {'key': key, 'targets': targets}
            name_sources[name].append(source)

        if self.fallback_renames:
            print(f'\n### {len(self.fallback_renames)} VCF file name(s) contain characters a '
                  f'tree cannot carry; those samples are named with "_" in their place:')
            for base, (before, after) in sorted(self.fallback_renames.items()):
                print(f'      {base}: {before!r} -> {after!r}')
            print()

        if self.metadata_ambiguous:
            print(f'\n### {len(self.metadata_ambiguous)} VCF file(s) match more than '
                  f'one row of the metadata worksheet, and those rows disagree on the '
                  f'name. The first matching row is used:')
            for base, info in sorted(self.metadata_ambiguous.items()):
                print(f'      {base}: "{info["key"]}" appears more than once, '
                      f'giving {", ".join(repr(t) for t in info["targets"])} '
                      f'-- using {info["targets"][0]!r}')
            print()

        # Two input files reducing to one name used to be fatal.  It is a mistake
        # in the metadata worksheet, not in the data, so the run now completes on
        # one file per name and the HTML summary carries the note -- a name clash
        # in a 7,992-sample sheet should not cost the whole analysis.  The file
        # used is the last of the colliding names in sorted order, chosen
        # explicitly so the result does not depend on directory listing order.
        self.name_collisions = {}
        superseded = set()
        for name, sources in sorted(name_sources.items()):
            if len(sources) > 1:
                ordered = sorted(sources, key=lambda s: os.path.basename(s))
                used, dropped = ordered[-1], ordered[:-1]
                bases = [os.path.basename(s) for s in ordered]
                # State the cause rather than guessing at one.  There are three,
                # they call for different corrections, and only the first is a
                # metadata problem at all.
                if len(set(bases)) < len(bases):
                    cause = ('two input files have the same file name, from '
                             'different directories')
                elif all(matched_key[s] is not None for s in ordered):
                    keys = {matched_key[s] for s in ordered}
                    cause = (f'the metadata worksheet gives this name to '
                             f'{" and ".join(sorted(repr(k) for k in keys))}')
                elif all(matched_key[s] is None for s in ordered):
                    cause = ('these file names reduce to the same name once '
                             '".vcf" and "_zc" are removed; no metadata row '
                             'matched either')
                else:
                    named = [os.path.basename(s) for s in ordered
                             if matched_key[s] is not None]
                    cause = (f'the metadata worksheet maps {", ".join(named)} to '
                             f'this name, which the other file already reduces to')
                self.name_collisions[name] = {
                    'used': os.path.basename(used),
                    'dropped': [os.path.basename(d) for d in dropped],
                    'cause': cause}
                superseded.update(dropped)
        if self.name_collisions:
            print(f'\n### {len(self.name_collisions)} sample name(s) resolved from '
                  f'more than one VCF. One file per name is analysed; the others are '
                  f'left out until the names are corrected:')
            for name, info in sorted(self.name_collisions.items()):
                print(f'      {name}: analysing {info["used"]}, leaving out '
                      f'{", ".join(info["dropped"])}')
                print(f'          because {info["cause"]}')
            print()

        for source, single_df in dataframes.items(): #just do once
            if source in superseded:
                continue
            sample = resolved[source]
            self.sample_sources[sample] = os.path.basename(source)
            # single_df.loc[single_df['ALT'].str.len() > 1, 'ALT'] = 'N' # keep indel positions Ns, ie. ALT indels to N, REF indels handled in make_groupings function
            dataframes_names_updated[sample] = single_df
            all_records = single_df
            # First iteration.  Find good SNPs for each VCF.  A sample needs no
            # good SNP of its own to be in the table: it is normalized against
            # the group's positions in make_groupings and carries the reference
            # call (or its own mixed/N/no-coverage call) at each of them.
            try:
                # A record with no MQ or QUAL is excluded from position selection
                # rather than read as zero: as zero these silently failed the
                # thresholds, so the positions vanished from selection.  Not
                # reported -- a zero-coverage record carries neither value by
                # construction (REF=N, ALT=., QUAL=.), so in a _zc VCF this is
                # most of the file and says nothing about the sample.  Those
                # records still reach the table as '-' through make_groupings;
                # what is excluded here is only their use in choosing positions.
                #
                # One pass over the sample's records.  The filtered selection
                # is these rows minus the filter positions (a row mask, so the
                # order of filtering and selecting does not matter), and the
                # counts below reuse the same result: three passes per sample
                # over every record, zero-coverage ones included, took about
                # 36 minutes of a 10,000-sample run.
                passing_all = self._passing(all_records)
            except AttributeError:
                print('\n### Error with sample {}\nSee VCF file and rerun\n'.format(sample))
                sys.exit(1)
            if not no_filters and filter_all_list:
                single_df = passing_all[~passing_all['abs_pos'].isin(filter_all_list)]
            else:
                single_df = passing_all

            mq_dictionary = single_df[['abs_pos', 'MQ']].set_index('abs_pos').to_dict()['MQ']
            for abs_pos, MQ in mq_dictionary.items():
                map_quality_dict[abs_pos].append(MQ)

            if single_df.empty:
                # Kept, and said.  This sample used to be appended to
                # vcf_bad_list as "Dataframe Empty" and left out of every
                # group, so the summary counted it as corrupt while the
                # validation section counted zero failures, and the sample was
                # nowhere in the tables.  A sample identical to the reference
                # looks exactly like this; so does one whose calls are all AC=1.
                passing_before_filter = 0
                if not no_filters and filter_all_list:
                    passing_before_filter = len(passing_all)
                self.no_snp_samples[sample] = self.describe_no_snps(
                    all_records, os.path.basename(source), passing_before_filter)
            else:
                # A mixture with a few passing SNPs (sites where one strain's reads
                # happen to dominate) is not caught above.  Flag it when its
                # confident mixed calls outnumber its passing SNPs several times
                # over (vsnp3_version.mixed_signature).
                passing_all_count = len(passing_all)
                confident_ac1 = self._confident_mixed(all_records)
                if mixed_signature(passing_all_count, confident_ac1):
                    self.mixed_samples[sample] = {
                        'file': os.path.basename(source), 'passing': passing_all_count,
                        'confident_ac1': confident_ac1}
            dataframe_essentials[sample] = single_df
            for abs_pos, alt in zip(single_df.abs_pos, single_df.ALT):
                abs_pos_nt_counts[abs_pos][alt] += 1
        # Name collisions were settled before this loop, so every remaining input
        # has its own name and nothing can be overwritten here.  Asserted rather
        # than assumed: the count silently disagreeing is what the old fatal check
        # existed to catch, and it would mean a sample vanished.
        expected = len(dataframes) - len(superseded)
        if len(dataframes_names_updated) != expected:
            raise RuntimeError(
                f'{expected} VCF files should have produced {expected} sample '
                f'names, but {len(dataframes_names_updated)} were filed. A sample '
                f'has been overwritten, which would drop it from the analysis '
                f'without further notice. This is a defect in vSNP3, not in the '
                f'input; please report it.')
        self.dataframe_essentials = dataframe_essentials
        self.dataframes_names_updated = dataframes_names_updated
        samples_with_dataframes_set = set(dataframe_essentials.keys())
        dataframes={} #empty memory
        mq_averages={}
        for abs_pos, mq_list in map_quality_dict.items():
            mq_averages[abs_pos] = np.mean(mq_list)
        self.average_mq_df = pd.DataFrame.from_dict(mq_averages, orient='index')
        self.average_mq_df = self.average_mq_df.reset_index()
        self.average_mq_df.columns = ['abs_pos', 'MQ']
        
        # Apply density filtering if enabled
        if self.filter_density:
            self.apply_density_filter()
        
        if gbk_list:
            annotation = Annotation(gbk_list=gbk_list)
            for abs_pos in sorted(abs_pos_nt_counts):
                counts = abs_pos_nt_counts[abs_pos]
                # Most frequently observed ALT wins, ties broken alphabetically, so
                # the result does not depend on which VCF was read last.
                snp_nt = min(counts, key=lambda alt: (-counts[alt], alt))
                annotation.run(abs_pos, snp_nt)
                if not getattr(annotation, 'feature_found', False):
                    annotation_dict[abs_pos] = 'position not annotated'
                elif annotation.reference_base_code == 'n/a':
                    # In a feature but no codon applies: tRNA, rRNA, ncRNA,
                    # repeat_region, mobile_element, a pseudogene, or an indel.
                    # Report what is known instead of discarding the annotation.
                    annotation_dict[abs_pos] = '{}, {}, {}'.format(
                        annotation.gene, annotation.product, annotation.mutation_type)
                else:
                    annotation_dict[abs_pos] = '{}->{}, {}:{}{}{}, {}, {}, codon_pos_{}'.format(
                        annotation.reference_base_code, annotation.snp_base_code,
                        annotation.gene, annotation.ref_aa, annotation.aa_residue_pos,
                        annotation.snp_aa, annotation.product, annotation.mutation_type, annotation.aa_pos)
                    if getattr(annotation, 'cds_overlap', 'n/a') != 'n/a':
                        # Short marker only; the annotated VCF carries the full list.
                        annotation_dict[abs_pos] += ', CDS overlap'
                if len(counts) > 1:
                    # More than one ALT is genuinely present at this position, so the
                    # single amino acid change above describes only the majority.  Kept
                    # terse because this row is rendered rotated at a fixed height and
                    # already truncates; the annotated VCF carries every consequence.
                    others = ','.join(sorted(counts))
                    annotation_dict[abs_pos] += f', multi-allelic [{others}]'
            if annotation.chroms_without_records:
                missing = ', '.join(sorted(annotation.chroms_without_records))
                unannotated = sum(1 for v in annotation_dict.values()
                                  if v == 'position not annotated')
                print(f'\nNote: no GenBank record was supplied for {missing}. '
                      f'{unannotated} position(s) are reported as unannotated. '
                      f'Add a gbk covering these contigs to annotate them.')
        self.annotation_df = pd.DataFrame(annotation_dict.items(), columns=['abs_pos', 'annotation'])
        print('\n\tGetting dataframe essentials  Selection Time: {}\n'.format(datetime.now() - self.startTime))
        self._phase('Select SNP positions and annotate')

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
                        # Expanded once per group, not once per sample: the argument
                        # depends only on `group`.  The worst group in the shipped
                        # Mbovis filter file expands to 73,945 positions, ~7 ms, and
                        # this was paying that for every sample in the group.
                        expanded_filter_list = self.list_expansion(self.group_filter_snps_dict[group])
                        sample_dict_group_filter={}
                        for sample, each_df in sample_dict.items():
                            each_df = each_df[~each_df['abs_pos'].isin(expanded_filter_list)] #by group remove positions to filter
                            sample_dict_group_filter[sample] = each_df
                        sample_dict = sample_dict_group_filter
                    if group in groupings_dict:
                        # Two defining positions name the same group: the group
                        # is their union.  The second used to overwrite the
                        # first, so the first position's samples left the
                        # group's table with nothing said.
                        merged = dict(groupings_dict[group])
                        merged.update(sample_dict)
                        groupings_dict[group] = merged
                    else:
                        groupings_dict[group] = sample_dict #defining_snps_dict[abs_pos] provides group
                    
            # else:
            #     print('Grouping: Pool processing with {} cpus...'.format(cpu_count))
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
            print('\n\tGroup Selection Time: {}\n'.format(datetime.now() - self.startTime))
            self._phase('Group selection')
        else:
            samples_without_group_set = samples_with_dataframes_set
            groupings_dict = {}

        if all_vcf or not defining_snps_dict:
            self.all_vcf_name = all_vcf_name or 'all_vcf'
            if all_vcf_name:
                groupings_dict[all_vcf_name] = self.dataframe_essentials
            else:
                groupings_dict['all_vcf'] = self.dataframe_essentials

        self.groupings_dict = groupings_dict
        print('All relevant positions by group')
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

        # A group with nothing to align still gets its directory and the same
        # marker file that a too-small group gets below.  Before the fail-loud
        # change these groups reached dict_to_fasta, wrote a '>name' header with no
        # sequence, and were caught by the line-count check there -- so the marker
        # is what users have been seeing for them.  Dropping the group silently
        # would leave no trace on disk that it had been considered at all.
        for skipped_group, _ in groupings_dict_list:
            if skipped_group not in finished_groupings_dict:
                os.makedirs(skipped_group, exist_ok=True)
                with open('{}/TOO_FEW_SAMPLES_OR_SHORT_SEQUENCE_TO_BUILD_TREE'.format(
                        skipped_group), 'w') as message_out:
                    print('no parsimony-informative positions: nothing to align, '
                          'so no tree was built', file=message_out)
        # else:
        #     with futures.ThreadPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## thread is diffently better but putting this in futures is slightly slower.
        #         for parsmony_sample_dict, group in pool.map(self.make_groupings, groupings_dict_list):
        #             finished_groupings_dict[group] = parsmony_sample_dict
        
        #Find positions that need to be filtered
        if self.find_new_filters:
            # The body used to be wrapped in `for df in group_sample_dict:`, which
            # iterated the 2-tuple left over from an earlier loop and never used
            # `df` -- so everything below ran exactly twice per group, doing the
            # work twice and rewriting the same two files.  Same leaked-loop-variable
            # class as the df_norm defect fixed in the fail-loud stage.
            for group_dict_of_df in groupings_dict_list:
                if len(group_dict_of_df[1]) > 3:
                    with open('{}_postion_list.txt'.format(group_dict_of_df[0]), 'w') as postion_list, \
                         open('{}_postion_detail_list.txt'.format(group_dict_of_df[0]), 'w') as postion_detail_list:
                        if not no_filters:
                            print('New positions to filter found after current filter positions applied but before noninformative SNP are removed', file=postion_detail_list)
                        else:
                            print('No previous filtering applied', file=postion_detail_list)
                        # Hardcoded values.  Looking to get a list of positions likely poor.  Curating findings before adding positions to filter shold be done post this step
                        print('dd.QUAL.mean() < 700 and dd.QUAL.max() < 1300 or dd.MQ.mean() < 40', file=postion_list)
                        print('dd.QUAL.mean() < 700 and dd.QUAL.max() < 1300 or dd.MQ.mean() < 40', file=postion_detail_list)
                        cc = pd.concat(group_dict_of_df[1].values(), ignore_index=True)
                        cc['abs_pos'] = cc['CHROM'] + ':' + cc['POS'].astype(str)
                        # One groupby instead of a full-frame boolean mask per
                        # position.  The mask form was O(positions x rows): 40 s at
                        # 40,000 rows, and a 200,000-row group did not finish in four
                        # minutes.  Sorted so the output order is stable -- the old
                        # form iterated a set, so it varied with PYTHONHASHSEED.
                        stats = cc.groupby('abs_pos').agg(
                            n=('QUAL', 'size'), qual_mean=('QUAL', 'mean'),
                            qual_max=('QUAL', 'max'), mq_mean=('MQ', 'mean'))
                        flagged = stats[(stats['n'] > 3)
                                        & (((stats['qual_mean'] < 700) & (stats['qual_max'] < 1300))
                                           | (stats['mq_mean'] < 40))]
                        for vv, row in flagged.sort_index().iterrows():
                            print(vv, file=postion_list)
                            print('{} Average QUAL: {:.2f}, Max QUAL: {:.2f}, Average MQ: {:.2f}'.format(
                                vv, row['qual_mean'], row['qual_max'], row['mq_mean']), file=postion_detail_list)
        print('\n\tAll relevant positions by group {}\n'.format(datetime.now() - self.startTime))
        self._phase('Normalize groups, parsimony-informative positions')

        print('FASTAs out and RAxML trees')
        self.startTime = datetime.now()
        group_fasta_dict={}
        group_dataframe_dict={}
        remove_list=[]
        for group, sample_dict in finished_groupings_dict.items():
            os.makedirs(group, exist_ok=True)
            fasta = '{}/{}-{}.fasta'.format(group, group, self.st)
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
                with open('{}/TOO_FEW_SAMPLES_OR_SHORT_SEQUENCE_TO_BUILD_TREE'.format(group), 'w') as message_out:
                    print('check sample numbers or sequence lengths', file=message_out)
                remove_list.append(group)
            else:
                group_fasta_dict[group] = fasta
                group_df = self.dict_to_dataframe(sample_dict, group)
                group_dataframe_dict[group] = group_df
        
        self.group_fasta_dict = group_fasta_dict
        self.group_dataframe_dict = group_dataframe_dict
        finished_groupings_list = list(finished_groupings_dict.keys())
        working_group_list = [x for x in finished_groupings_list if x not in remove_list]
        # Largest alignment first.  The all-VCF group is the last one defined
        # and by far the biggest, so it used to be submitted last and the
        # phase ended with it running alone.  Each group's output does not
        # depend on when it is built; a stable sort keeps the old order among
        # groups of the same size.
        working_group_list.sort(key=lambda g: -group_dataframe_dict[g].size)
        self.calculate_average_coverage()

        # One group's failure must not take the others with it.  pool.map returns a
        # generator, so an exception raised in any group propagates out of the for
        # loop and every group not yet iterated is silently abandoned -- which made
        # the fail-loud changes elsewhere in this stage strictly worse without this.
        def build_one(group):
            try:
                return self.raxml_table_build(group)
            except Exception as e:                              # noqa: BLE001
                self.group_failures[group] = f'{type(e).__name__}: {e}'
                return None

        if debug:
            for group in working_group_list:
                build_one(group)
        else:
            with futures.ThreadPoolExecutor(max_workers=cpu_count) as pool: #ProcessPoolExecutor ThreadPoolExecutor ## thread works best for raxml calls
                for tree in pool.map(build_one, working_group_list):
                    pass

        if self.group_failures:
            print(f'\n### {len(self.group_failures)} of {len(working_group_list)} '
                  f'group(s) produced no table or tree:')
            for group, reason in sorted(self.group_failures.items()):
                print(f'###   {group}: {reason}')
            print('###   The remaining groups completed normally.\n')
        if self.html_tree_failures:
            print(f'\n### HTML tree not written for {len(self.html_tree_failures)} group(s); '
                  f'their tables and RAxML trees were:')
            for group, reason in sorted(self.html_tree_failures.items()):
                print(f'###   {group}: {reason}')
            print()

        print('\n\tFASTAs, RAxML and HTML trees {}\n'.format(datetime.now() - self.startTime))
        self._phase('FASTAs, RAxML trees, tables and HTML trees')
        # print('\n\nTotal Time: {}\n'.format(datetime.now() - self.beginTime))

        #Add back those that where a group was not found
        if 'Group Not Found' not in groupings_dict:
            groupings_dict['Group Not Found'] = {}

        for sample in samples_without_group_set:
            groupings_dict['Group Not Found'][sample] = pd.DataFrame()
        self.groupings_dict = groupings_dict # will be passed to html summary

    def apply_density_filter(self):
        """Apply density filtering to remove SNPs when threshold+ SNPs are found within specified window"""
        print("Applying density filtering (>={} SNPs within {} bp window)...".format(
            self.density_threshold, self.density_window))
        
        # Log file for density filtering
        density_log = open('{}/density_filtering_log.txt'.format(self.cwd), 'w', encoding='utf-8')
        print("Density filtering applied: SNPs filtered when >={} SNPs found within {} bp window".format(
            self.density_threshold, self.density_window), file=density_log)
        print("Filtering applied on: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')), file=density_log)
        print("=" * 60, file=density_log)
        
        positions_to_remove = set()
        total_positions_examined = 0
        regions_filtered = 0
        
        # Group positions by chromosome.  Vectorized: this was `for _, row in
        # df.iterrows()` per sample at ~12 us a row, which is ~29 s for 500 samples
        # across 5,000 positions -- more than the window scan below used to cost.
        chrom_positions = defaultdict(set)
        for sample_name, df in self.dataframe_essentials.items():
            if df.empty:
                continue
            split = df['abs_pos'].str.rsplit(':', n=1, expand=True)
            for chrom, group in split.groupby(split[0])[1]:
                chrom_positions[chrom].update(group.astype(int).tolist())

        # Process each chromosome separately
        for chrom, positions in chrom_positions.items():
            unique_positions = sorted(positions)
            total_positions_examined += len(unique_positions)

            print("Examining chromosome {} with {} positions".format(chrom, len(unique_positions)), file=density_log)

            # Two-pointer sweep over the sorted positions instead of rescanning all
            # of them for every window.  The old form was O(P^2): 3.96 s at 20,000
            # positions and 15.9 s at 40,000.  `bisect_right` finds the end of each
            # forward window directly, so the same windows are produced in the same
            # order and the log output is unchanged.
            end = 0
            for start, pos in enumerate(unique_positions):
                window_end = pos + self.density_window - 1  # -1 because we want inclusive range
                if end < start:
                    end = start
                end = bisect.bisect_right(unique_positions, window_end, start)
                snps_in_window = unique_positions[start:end]

                # If threshold+ SNPs in window, mark all for removal
                if len(snps_in_window) >= self.density_threshold:
                    regions_filtered += 1
                    print("Dense region found: positions {} (window {}-{})".format(
                        snps_in_window, pos, window_end), file=density_log)
                    for dense_pos in snps_in_window:
                        abs_pos = "{}:{}".format(chrom, dense_pos)
                        positions_to_remove.add(abs_pos)
                        self.density_filtered_positions.append(abs_pos)
        
        # Remove dense positions from all samples
        original_counts = {}
        filtered_counts = {}
        
        for sample_name, df in self.dataframe_essentials.items():
            original_counts[sample_name] = len(df)
            filtered_df = df[~df['abs_pos'].isin(positions_to_remove)]
            self.dataframe_essentials[sample_name] = filtered_df
            filtered_counts[sample_name] = len(filtered_df)
        
        # Log summary statistics
        print("=" * 60, file=density_log)
        print("SUMMARY:", file=density_log)
        print("Filtering criteria: >={} SNPs within {} bp window".format(
            self.density_threshold, self.density_window), file=density_log)
        print("Total positions examined: {}".format(total_positions_examined), file=density_log)
        print("Dense regions identified: {}".format(regions_filtered), file=density_log)
        print("Positions filtered: {}".format(len(positions_to_remove)), file=density_log)
        print("", file=density_log)
        print("Per-sample filtering results:", file=density_log)
        
        for sample_name in original_counts:
            removed = original_counts[sample_name] - filtered_counts[sample_name]
            print("{}: {} -> {} (removed {})".format(
                sample_name, original_counts[sample_name], filtered_counts[sample_name], removed), file=density_log)
        
        density_log.close()
        
        print("Density filtering complete: {} positions removed from {} dense regions".format(
            len(positions_to_remove), regions_filtered))
        print("Detailed log written to: {}/density_filtering_log.txt".format(self.cwd))

    def calculate_average_coverage(self):
        """Calculate average depth of coverage for each sample, excluding zero values."""
        for sample, single_df in self.dataframes_names_updated.items():
            try:
                # Convert to numeric and exclude zeros
                coverage_values = pd.to_numeric(single_df['DP'], errors='coerce')
                # Filter out zeros and NaN values before calculating mean
                valid_coverage = coverage_values[coverage_values > 0]
                mean_coverage = valid_coverage.mean()
                # If all values were zero/NaN, mean_coverage will be NaN
                self.sample_coverage_dict[sample] = int(round(mean_coverage)) if pd.notna(mean_coverage) else 0
            except (KeyError, AttributeError):
                # Handle cases where DP column might not exist
                self.sample_coverage_dict[sample] = 0

    def _position_sets(self):
        """
        {sample: set of abs_pos} for every sample, built once.

        group_selection() is called once per defining SNP and needs each sample's
        position set to test membership.  It used to rebuild all of them on every
        call, so the work scaled with defining-SNPs x samples x positions for a
        result that never changes.
        """
        cached = getattr(self, '_position_sets_cache', None)
        if cached is None:
            cached = {sample: set(df['abs_pos'].values)
                      for sample, df in self.dataframe_essentials.items()}
            self._position_sets_cache = cached
        return cached

    def group_selection(self, abs_pos):
        sample_dict = {}
        group_found = False
        abs_pos = abs_pos.split(", ")
        
        # Separate normal and inverted positions
        normal_positions = []
        inverted_positions = []
        
        for pos in abs_pos:
            if pos.endswith("!"):
                # Remove the "!" and add to inverted positions
                inverted_positions.append(pos[:-1])
            else:
                normal_positions.append(pos)
        
        for sample, single_df in self.dataframe_essentials.items():
            # Cached across calls: group_selection runs once per defining SNP (197 in
            # the shipped Mbovis file) and rebuilt every sample's position set each
            # time, though the sets never change.  Measured 2 s at 200 samples.
            file_positions = self._position_sets().get(sample)
            if file_positions is None:
                file_positions = set(single_df['abs_pos'].values)

            # Check if sample qualifies for the group
            qualifies = True
            
            # All normal positions must be present in the file
            if normal_positions:
                if not set(normal_positions).issubset(file_positions):
                    qualifies = False
            
            # All inverted positions must NOT be present in the file
            if inverted_positions and qualifies:
                for inv_pos in inverted_positions:
                    if inv_pos in file_positions:
                        qualifies = False
                        break
            
            # If sample qualifies, add it to the group
            if qualifies and (normal_positions or inverted_positions):
                group_found = True
                sample_dict[sample] = self.dataframe_essentials[sample]
        
        return group_found, sample_dict

    def _check_defining_snps(self, path, defining_snps_df, defining_snps_dict,
                             group_filter_snps_dict, filter_all_list):
        '''
        The worksheet mistakes that used to crash the run after the inputs
        were archived, or mis-group silently, caught once after reading.

        Raises DefiningSnpsError for a position header that is not
        chrom:position text and for a blank group name.  A group name used by
        more than one position is allowed (the group is the union; see where
        groupings_dict is filled) and the filter lists of its columns, which
        pandas reads back as "name", "name.1", ..., are merged here.  Filter
        entries that are not chrom:position text match no record and are
        noted, as they always were ignored.
        '''
        name = os.path.basename(str(path))
        problems = []
        positions = list(defining_snps_df.columns[1:])
        for number, label in enumerate(positions, start=2):
            if not isinstance(label, str) or ':' not in label:
                problems.append(f'column {number}: the position header {label!r} is not '
                                f'chrom:position text')
        groups = {}
        for label in positions:
            if not isinstance(label, str):
                continue
            grp = defining_snps_dict.get(label)
            if grp is None or (isinstance(grp, float) and grp != grp) or not str(grp).strip():
                problems.append(f'position {label}: the group name in row 2 is blank')
            else:
                groups.setdefault(str(grp), []).append(label)
        if problems:
            raise DefiningSnpsError(
                f'Defining SNP worksheet {name} cannot be used:\n  '
                + '\n  '.join(problems)
                + '\nRow 1 holds chrom:position headers, row 2 the group names, and the '
                  'rows below them the positions to filter. Correct the worksheet and rerun; '
                  'the input VCF files have not been removed.')
        duplicates = sorted(g for g, cols in groups.items() if len(cols) > 1)
        if duplicates:
            self.worksheet_notes.append(
                f'{len(duplicates)} group name(s) are given by more than one defining '
                f'position and were merged into one group each: '
                + ', '.join(f'{g} ({len(groups[g])} positions)' for g in duplicates))
            # pandas names repeated columns "G", "G.1", "G.2": fold them back.
            merged = {}
            for column, entries in group_filter_snps_dict.items():
                base = re.sub(r'\.\d+$', '', str(column))
                target = base if (base in duplicates and column not in groups) else column
                merged.setdefault(target, [])
                merged[target].extend(e for e in entries if e not in merged[target])
            group_filter_snps_dict = merged
        odd = [(str(g), e) for g, entries in group_filter_snps_dict.items() for e in entries
               if not (isinstance(e, str) and ':' in e)]
        odd += [('(all)', e) for e in filter_all_list if not (isinstance(e, str) and ':' in e)]
        if odd:
            shown = ', '.join(f'{g}: {e!r}' for g, e in odd[:8])
            self.worksheet_notes.append(
                f'{len(odd)} filter entr{"y" if len(odd) == 1 else "ies"} in {name} are not '
                f'chrom:position text and match no record, so they filter nothing: {shown}'
                + (' ...' if len(odd) > 8 else ''))
        for note in self.worksheet_notes:
            print(f'\n### {note}\n')
        return group_filter_snps_dict

    def _phase(self, name):
        now = time.monotonic()
        self.phase_times.append((name, now - self._phase_t0))
        self._phase_t0 = now

    @staticmethod
    def _single_base(df):
        '''
        Mask of records whose REF and ALT are both one base.  The reader
        (vsnp3_step2.dataframe_from_vcf_body) writes it as the single_base
        column once per file; the two str.len() passes are the fallback for a
        frame that did not come from the reader.
        '''
        if 'single_base' in df.columns:
            return df['single_base'].astype(bool)
        return (df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1)

    def _passing(self, df):
        '''The records that qualify as SNP positions: the 3.36 selection, unchanged.'''
        return df[df['QUAL'].notna() & (df['QUAL'] > self.qual_threshold)
                  & (df['AC'] == AC_HOMOZYGOUS)
                  & self._single_base(df)
                  & df['MQ'].notna() & (df['MQ'] >= self.mq_threshold)]

    def _confident_mixed(self, df):
        '''Single-base AC=1 records at QUAL >= qual_threshold and MQ >= mq_threshold.'''
        return int((df['QUAL'].notna() & (df['QUAL'] >= self.qual_threshold)
                    & (df['AC'] == AC_HETEROZYGOUS)
                    & self._single_base(df)
                    & df['MQ'].notna() & (df['MQ'] >= self.mq_threshold)).sum())

    @staticmethod
    def _allele_bases(ref, alt, cigar):
        '''
        {reference offset: base} for one ALT allele of a multi-base record, or
        None when the allele cannot be placed on the reference.

        The CIGAR freebayes writes for each allele says which reference base each
        ALT base aligns to (M, X, = consume both; D deletes reference bases, which
        map to None; I inserts ALT bases).  Without a CIGAR that accounts for every
        base of REF and ALT, a same-length allele is read base by base; any other
        allele cannot be placed.
        '''
        if isinstance(cigar, str) and cigar:
            ops = re.findall(r'(\d+)([MIDX=])', cigar)
            if ops and ''.join(n + op for n, op in ops) == cigar:
                bases, r, a = {}, 0, 0
                for n, op in ops:
                    n = int(n)
                    if op in 'MX=':
                        for i in range(n):
                            bases[r + i] = alt[a + i] if a + i < len(alt) else None
                        r += n
                        a += n
                    elif op == 'D':
                        for i in range(n):
                            bases[r + i] = None
                        r += n
                    else:
                        a += n
                if r == len(ref) and a == len(alt):
                    return bases
        if len(ref) == len(alt):
            return dict(enumerate(alt))
        return None

    def _embedded_calls(self, sample):
        '''
        One row per reference position inside the sample's multi-base records
        (MNP, complex or indel records whose REF spans more than one base), other
        than the record's first position, where the sample's base is not the
        reference: abs_pos, REF (that reference base), ALT (the sample's base, or
        N where it is deleted, cannot be placed, or the record's alleles
        disagree), and the record's QUAL, AC and MQ.

        freebayes reports a SNP inside such a record when a nearby read carries an
        indel, e.g. REF=CAC ALT=CCC at position p for an A>C change at p+1.  The
        SNP's own position then has no record, and make_groupings used to write
        the reference base there.  These rows give it the base the record reports,
        under the same QUAL, AC and MQ rules as any other call.
        '''
        cache = self.__dict__.setdefault('_embedded_cache', {})
        if sample in cache:
            return cache[sample]
        df = self.dataframes_names_updated[sample]
        multi = df[(df['REF'].str.len() > 1) & (df['ALT'] != '.')]
        rows = {}
        for rec in multi.to_dict('records'):
            ref, alleles = rec['REF'], str(rec['ALT']).split(',')
            cigar = rec.get('CIGAR')
            cigars = cigar.split(',') if isinstance(cigar, str) else []
            if len(cigars) != len(alleles):
                cigars = [None] * len(alleles)
            placed = [self._allele_bases(ref, a, c) for a, c in zip(alleles, cigars)]
            for offset in range(1, len(ref)):
                ref_base = ref[offset]
                seen = set()
                for bases in placed:
                    if bases is None:
                        seen.add('?')                     # allele cannot be placed
                    else:
                        b = bases.get(offset)
                        seen.add('-' if b is None else b)  # '-' a deleted base
                if seen == {ref_base}:
                    continue                              # every allele keeps the reference base
                base = seen.pop() if len(alleles) == 1 and len(seen) == 1 else 'N'
                if base not in ('A', 'C', 'G', 'T'):
                    base = 'N'
                key = f"{rec['CHROM']}:{int(rec['POS']) + offset}"
                if key in rows:
                    if rows[key]['ALT'] != base:
                        rows[key]['ALT'] = 'N'            # two records disagree here
                    continue
                rows[key] = {'abs_pos': key, 'REF': ref_base, 'ALT': base,
                             'QUAL': rec['QUAL'], 'AC': rec['AC'], 'MQ': rec['MQ']}
        out = pd.DataFrame(list(rows.values()), columns=['abs_pos', 'REF', 'ALT', 'QUAL', 'AC', 'MQ'])
        out['QUAL'] = out['QUAL'].astype('Int64')
        out['MQ'] = out['MQ'].astype('Int64')
        out['AC'] = out['AC'].astype(np.int64)
        cache[sample] = out
        return out

    def describe_no_snps(self, df, basename, passing_before_filter=0):
        '''
        Why a sample has no SNP passing the thresholds, as counts.

        Written for the summary report, so that a sample with an empty row can
        be told apart at a glance: identical to the reference, no coverage,
        every call mixed (AC=1), every call an indel, or every call below QUAL
        or MQ.
        '''
        detail = {'file': basename, 'records': int(len(df)), 'called': 0, 'ac1': 0,
                  'mixed_only': False, 'reason': ''}
        if len(df) == 0:
            detail['reason'] = ('no records in the VCF: identical to the reference '
                                'everywhere, with full coverage')
            return detail
        called_mask = df['QUAL'].notna()
        called = int(called_mask.sum())
        zero_cov = len(df) - called
        detail['called'] = called
        if called == 0:
            detail['reason'] = (f'{len(df):,} records, all zero-coverage, no variant '
                                f'called: identical to the reference where covered')
            return detail
        cdf = df[called_mask]
        ref_len = cdf['REF'].str.len()
        alt_len = cdf['ALT'].str.len()
        indel = int(((ref_len > 1) | (alt_len > 1)).sum())
        ac1 = int((cdf['AC'] == AC_HETEROZYGOUS).sum())
        snp_mask = (cdf['AC'] == AC_HOMOZYGOUS) & (ref_len == 1) & (alt_len == 1)
        snps = int(snp_mask.sum())
        low_qual = int((snp_mask & ~(cdf['QUAL'] > self.qual_threshold)).sum())
        low_mq = int((snp_mask & (cdf['MQ'].isna() | (cdf['MQ'] < self.mq_threshold))).sum())
        detail['ac1'] = ac1
        # Called records that are none of the above: multi-allelic (AC=1,1),
        # or AC missing from INFO.
        other = max(0, called - ac1 - indel - snps)
        parts = [f'{len(df):,} records: {zero_cov:,} zero-coverage, {called:,} called',
                 f'{ac1:,} mixed (AC=1)', f'{indel:,} indel']
        if snps:
            parts.append(f'{snps:,} single-base AC=2, of which {low_qual:,} have QUAL <= '
                         f'{self.qual_threshold} and {low_mq:,} have MQ < {self.mq_threshold}')
        else:
            parts.append('0 single-base AC=2')
        if other:
            parts.append(f'{other:,} other (multi-allelic or no AC)')
        reason = '; '.join(parts) + '; 0 passed the thresholds'
        if passing_before_filter:
            reason += (f' ({passing_before_filter:,} passed but every one is at a '
                       f'position in the filter list)')
        detail['reason'] = reason
        # The mixed-sample signature: no homozygous single-base SNP at all, and
        # mixed calls make up nearly every non-indel call.  Nearly, not all: a
        # real sample of this kind usually has a stray multi-allelic record.
        non_indel = max(1, called - indel)
        detail['mixed_only'] = snps == 0 and ac1 > 0 and ac1 >= 0.9 * non_indel
        return detail

    @staticmethod
    def read_metadata(metadata):
        '''
        The metadata worksheet as a two-column frame: file_name, metadata.

        Names are cleaned to what the file system and RAxML accept, as before.
        '''
        metadata_df = pd.read_excel(metadata, header=None, index_col=0, usecols=[0, 1], names=['file_name', 'metadata'])

        metadata_df = metadata_df.reset_index()
        # A blank cell is not a name.  astype(str) turned a blank name cell
        # into the string 'nan', so one blank renamed its sample to "nan" and
        # two blanks collided and dropped a sample under "duplicate sample
        # name: nan".  Rows with a blank cell in either column are left out
        # (the VCF file name then names the sample, as for any unmatched
        # file) and counted for the summary.
        blank = (metadata_df['file_name'].isna() | metadata_df['metadata'].isna()
                 | (metadata_df['file_name'].astype(str).str.strip() == '')
                 | (metadata_df['metadata'].astype(str).str.strip() == ''))
        blank_rows = int(blank.sum())
        metadata_df = metadata_df[~blank].copy()
        # Convert to string types in case sample names are numbers/int
        metadata_df['file_name'] = metadata_df['file_name'].astype(str)
        metadata_df['metadata'] = metadata_df['metadata'].astype(str)

        # Clean unicode characters first
        metadata_df['metadata'] = metadata_df['metadata'].apply(clean_unicode)

        # Define all illegal characters that need to be replaced
        illegal_chars = ['\t', '\r', ' ', ':', ',', ')', '(', ';', ']', '[', "'",
                        '/', '.', '*', '?', '{', '}']

        # Replace each illegal character with underscore
        for char in illegal_chars:
            metadata_df['metadata'] = metadata_df['metadata'].str.replace(char, '_', regex=False)

        # Clean up multiple underscores and trailing characters
        metadata_df['metadata'] = metadata_df['metadata'].str.replace('_+', '_', regex=True)  # Multiple underscores to single
        metadata_df['metadata'] = metadata_df['metadata'].str.replace(r'_$', '', regex=True)   # Trailing underscore
        metadata_df['metadata'] = metadata_df['metadata'].str.replace(r'-$', '', regex=True)   # Trailing dash
        metadata_df.attrs['blank_rows'] = blank_rows
        return metadata_df

    @staticmethod
    def safe_fallback_name(name):
        '''
        A file-name-derived sample name as a Newick taxon.  Characters RAxML
        cannot carry become '_'; a name that needs no change is returned as is.
        '''
        cleaned = NEWICK_UNSAFE.sub('_', str(name))
        if cleaned == name:
            return name
        cleaned = re.sub('_+', '_', cleaned).strip('_')
        return cleaned or name

    @staticmethod
    def resolve_sample_name_detail(basename, metadata_df, metadata_test):
        '''
        (name, metadata key matched or None, number of metadata rows it matched).

        The key and the row count are returned so a caller can say why two files
        ended up with one name, and can see when one file matched several rows of
        the worksheet -- a lookup takes the first row, so the rest are discarded
        with nothing said about it.

        Lifted from the nested try/except chain this replaces, including two
        details worth stating because they are easy to "tidy" into a behaviour
        change.  The substitutions are cumulative: each failed lookup leaves its
        strip in place, so the fallback is the basename with '.vcf', then '_zc',
        then '_zc_*' removed in that order.  And the first pattern is '.vcf$',
        where '.' is the regex any-character -- kept as it was, because widening
        it to r'\\.vcf$' would change which names strip.
        '''
        sample = basename
        if not metadata_test:
            return sample, None, 0
        for strip in (None, '.vcf$', '_zc$', '_zc_.*$'):
            if strip is not None:
                sample = re.sub(strip, '', sample)
            hits = metadata_df.loc[metadata_df['file_name'] == sample, 'metadata']
            if len(hits):
                return hits.iloc[0], sample, len(hits)
        return sample, None, 0

    @staticmethod
    def resolve_sample_name(basename, metadata_df, metadata_test):
        '''The name a VCF is filed under.  See resolve_sample_name_detail.'''
        return Group.resolve_sample_name_detail(
            basename, metadata_df, metadata_test)[0]

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
                    raise type(e)(str(e) + ' \n#### error in Defining SNPs/Filter worksheet\n#### see value "{}"'.format(list_entry)).with_traceback(sys.exc_info()[2])
                list_entry = sequence_range.split("-")
                for position in range(int(list_entry[0].replace(',', '')), int(list_entry[1].replace(',', '')) + 1):
                    expanded_list.append(chrom + ":" + str(position))
        return expanded_list

    def dict_to_fasta(self, sample_dict, fasta): #sample_dict = [file name]:snp dataframe
        '''
        Write the alignment, then assert it is actually an alignment.

        RAxML requires equal-length sequences.  Previously a KeyError here was
        swallowed after the '>name' header had already been written, leaving a
        header with no sequence line, and an unequal-length alignment surfaced only
        as 'check sample numbers' in a SEE_RAXML_INFO file that nothing reads.
        '''
        sequences = {}
        for name, df in sample_dict.items():
            try:
                df = df.sort_values(by=['abs_pos']) # sorting ensures positions are aligned
                sequences[name] = "".join(df['ALT'].to_list())
            except KeyError as e:
                raise RuntimeError(
                    f'could not build an alignment sequence for {name}: missing '
                    f'column {e}') from e

        lengths = {name: len(seq) for name, seq in sequences.items()}
        if len(set(lengths.values())) > 1:
            expected = Counter(lengths.values()).most_common(1)[0][0]
            odd = {n: l for n, l in lengths.items() if l != expected}
            raise RuntimeError(
                f'{os.path.basename(fasta)} is not an alignment: most sequences are '
                f'{expected} positions, but {odd} differ. RAxML requires equal '
                f'lengths, so no tree could be built from this.')
        if 'root' not in sequences:
            raise RuntimeError(
                f'{os.path.basename(fasta)} has no "root" record, and RAxML is called '
                f'with -o root, so it would fail. Samples present: '
                f'{sorted(sequences)}')

        with open(fasta, 'w') as write_out:
            for name, seq in sequences.items():
                print('>{}'.format(name), file=write_out)
                print(seq, file=write_out)

    def dict_to_dataframe(self, sample_dict, group): # dict_to_fasta will not have abs_pos this dataframe will retain abs_pos
        try:
            dict_of_dfs = { sample: df.set_index('abs_pos') for sample, df in sample_dict.items()} #set ALT column to sample name and merge
            group_df = pd.concat(dict_of_dfs, axis=1)
            group_df.columns = group_df.columns.droplevel(-1)
            return group_df.T
        except pd.errors.InvalidIndexError as e:
            if self.debug:
                print('\n\t#####\n\t##### {}, Group {} \n\t#####\n'.format(e, group))
            return pd.DataFrame()

    def make_groupings(self, group_sample_dict):
        df_list=[]
        group, sample_dict = group_sample_dict # sample_dict is from dataframe_essentials, ie good SNPs.
        for sample, df in sample_dict.items():
            if not df.empty:
                df_list.append(df[['abs_pos', 'REF']])
        if not df_list:
            # Every sample in the group is reference-identical, or has no SNP
            # passing the thresholds: there is no position to align.
            print(f'  {group}: no sample has a SNP passing the thresholds. Nothing to '
                  f'align, so no tree; group skipped and the rest of the run continues.')
            return {}, group
        df = pd.concat(df_list, ignore_index=True)
        df_ref = df.drop_duplicates(subset='abs_pos', keep="first") # df_ref will be all abs_pos per group.  completion of first iteration
        norm_sample_dict={}
        # The parsimony test, kept O(positions) in memory.  It used to extend
        # two lists with one abs_pos+ALT string per sample per position (one of
        # them never read) and count them: about 27 GB for an all-VCF group of
        # 10,000 samples x 20,000 positions.  A position is uninformative when
        # every sample carries the same single-base call there, which is "no
        # sample differs from the first" -- the same test, as a boolean mask.
        undetermined = []       # samples with no determined base at the group's positions
        embedded_added = {}   # sample -> positions whose call came from a multi-base record
        # Set in __init__; created here too for a Group built without it (tests, pickles).
        if not hasattr(self, 'embedded_cells'):
            self.embedded_cells = defaultdict(int)
            self.embedded_now_uninformative = defaultdict(set)

        def normalize_one(sample, quiet=False):
            '''
            One sample over the group's positions: (normalized calls as an
            abs_pos/ALT frame, the parsimony-test calls as a Series indexed by
            abs_pos, the positions whose call came from a multi-base record).
            quiet suppresses the messages and records for the rare second pass
            below, which repeats the same arithmetic.
            '''
            embedded_positions = set()
            # dataframes_names_updated contains all SNPs
            # dataframe_essentials contains good SNP positions
            sample_df = self.dataframes_names_updated[sample]
            # .copy() because the boolean mask returns a view under some pandas
            # versions and the block below writes into it via .loc and .at.  The
            # SettingWithCopyWarning that would have said so is suppressed module-wide.
            sample_df = sample_df[sample_df['abs_pos'].isin(df_ref['abs_pos'])].copy() # this will normalize positions
            #https://stackoverflow.com/questions/27673231/why-should-i-make-a-copy-of-a-data-frame-in-pandas
            # A group position with no record of its own may lie inside one of the
            # sample's multi-base records.  Its base is read from that record
            # (_embedded_calls) instead of falling through to the reference below.
            embedded = self._embedded_calls(sample)
            if not embedded.empty:
                add = embedded[embedded['abs_pos'].isin(df_ref['abs_pos'])
                               & ~embedded['abs_pos'].isin(sample_df['abs_pos'])]
                if not add.empty:
                    sample_df = pd.concat([sample_df, add], ignore_index=True)
                    embedded_positions = set(add['abs_pos'])
            sample_df_parse_test = sample_df.copy() # to use below
            sample_df.loc[sample_df['ALT'].str.len() > 1, 'ALT'] = 'N'
            sample_df.loc[sample_df['ALT'] == 'N', 'AC'] = AC_HOMOZYGOUS # allow the above line to pass ambigious if needed
            sample_df['ALT'] = sample_df['ALT'].replace('.', '-')
            # One try per row.  The try used to wrap the whole loop, so the first
            # REF/ALT pair missing from ambigious_lookup abandoned every remaining
            # AC=1 row in that sample -- and those rows kept their raw ALT, which
            # publishes a heterozygous site as a confident homozygous call.  An
            # unrecognised pair now becomes N, matching how the lines below already
            # handle uncertainty, and the count is reported rather than hidden
            # behind --debug.
            # Vectorized.  This was an iterrows loop with a try/except per row,
            # the slowest Python loop of the stage on groups with many mixed
            # calls.  Same result: a known REF/ALT pair becomes its IUPAC code,
            # anything else becomes N and is counted.
            unresolved = []
            het = sample_df['AC'] == AC_HETEROZYGOUS
            if het.any():
                het_rows = sample_df.loc[het]
                keys = het_rows['REF'].fillna('').astype(str) + het_rows['ALT'].fillna('').astype(str)
                mapped = keys.map(self.ambigious_lookup)
                missing = mapped.isna()
                if missing.any():
                    unresolved = [f"{row.abs_pos} {row.REF}/{row.ALT}"
                                  for row in het_rows.loc[missing].itertuples()]
                sample_df.loc[het, 'ALT'] = mapped.fillna('N')
            if unresolved and not quiet:
                print(f'  {sample}: {len(unresolved)} heterozygous position(s) had no '
                      f'IUPAC code for their REF/ALT pair and were called N '
                      f'(first: {unresolved[0]})')
            # A missing QUAL is a no-call, applied BEFORE the quality bands below.
            # Previously QUAL was filled with 0, which put these positions under
            # n_threshold and rewrote them to the reference base - turning "we do not
            # know" into "this sample matches the reference".
            sample_df.loc[sample_df['QUAL'].isna() & (sample_df['ALT'] != '-'), 'ALT'] = 'N'
            #change alt to N if QUAL 50 - 150
            sample_df.loc[(sample_df['QUAL'] >= self.n_threshold) & (sample_df['QUAL'] < self.qual_threshold) & (sample_df['ALT'] != '-'), 'ALT'] = 'N' # this will overwrite ambigious calls
            # < 50 will default to REF... change ALT to REF
            try:
                sample_df.loc[sample_df['REF'].str.len() > 1, 'REF'] = 'N' #if REF call is indel change to N to maintain equal sequence length for all samples
                # notna() so an unknown QUAL does not fall into the "below threshold,
                # therefore reference" branch.
                mask = sample_df['QUAL'].notna() & (sample_df['QUAL'] < self.n_threshold) & (sample_df['ALT'] != '-')
                sample_df.loc[mask, 'ALT'] = sample_df['REF']
            except (ValueError) as e:
                # Said and recorded, not hidden behind --debug: if this ever
                # fires, positions below the N threshold keep their ALT and
                # would be published as confident calls.
                if not quiet:
                    self.__dict__.setdefault('rewrite_errors', {})[(group, sample)] = f'{type(e).__name__}: {e}'
                    print(f'\n### {group} / {sample}: the call rewrite failed ({e}); low-quality '
                          f'positions of this sample may keep their ALT. See the summary.\n')

            sample_df = sample_df[['abs_pos', 'ALT']] # no longer need other columns
            # sample_df = sample_df.replace(np.nan, '-') # change zero coverage to -
            df_merged = sample_df.merge(df_ref, left_on='abs_pos', right_on='abs_pos', how='outer') # finish normalizing if df doesn't include all position in group
            df_merged['ALT'] = np.where(df_merged['ALT'].notna(), df_merged['ALT'], df_merged['REF']) # merge REF from df_ref to ALT column
            df_merged['ALT'] = df_merged['ALT'].fillna('-')

            ####
            #Do not include low quality calls when determining if position is parsimonious.  Only include calls with QUAL >= qual_threshold [150]
            sample_df_parse_test['ALT'] = sample_df_parse_test['ALT'].replace('.', '-')
            try:
                sample_df_parse_test.loc[sample_df_parse_test['REF'].str.len() > 1, 'REF'] = 'N' #if REF call is indel change to N to maintain equal sequence length for all samples
                sample_df_parse_test.loc[sample_df_parse_test['QUAL'] < self.qual_threshold, 'ALT'] = sample_df_parse_test['ALT'] # change n_threshold from above to qual threshold skipping the Ns when determining if position is parsimonious AND changed to 'ALT'.  So, if there are just a few low quality represented the SNP position will be seen as parisomonious uninformative and removed.
            except (ValueError) as e:
                if not quiet:
                    self.__dict__.setdefault('rewrite_errors', {})[(group, sample)] = f'{type(e).__name__}: {e}'
                    print(f'\n### {group} / {sample}: the parsimony-test rewrite failed ({e}). '
                          f'See the summary.\n')
            sample_df_parse_test = sample_df_parse_test[['abs_pos', 'ALT']] # no longer need other columns
            # sample_df_parse_test = sample_df_parse_test.replace(np.nan, '-') # change zero coverage to -
            df_merged_parse_test = sample_df_parse_test.merge(df_ref, left_on='abs_pos', right_on='abs_pos', how='outer') # finish normalizing if df doesn't include all position in group
            df_merged_parse_test['ALT'] = np.where(df_merged_parse_test['ALT'].notna(), df_merged_parse_test['ALT'], df_merged_parse_test['REF']) # merge REF from df_ref to ALT column
            df_merged_parse_test['ALT'] = df_merged_parse_test['ALT'].fillna('-')
            calls = df_merged_parse_test.set_index('abs_pos')['ALT']
            ####
            return df_merged[['abs_pos', 'ALT']], calls, embedded_positions

        def parsimony_uninformative(calls_of_samples):
            '''
            The positions where every sample carries the same call.  The old
            count-based test built one string per sample per position and
            dropped its last character to recover the position, so it only ever
            removed positions whose shared call was one character long; that
            condition is kept explicitly, so the result is unchanged.
            '''
            first_calls = differs = None
            for calls in calls_of_samples:
                if first_calls is None:
                    first_calls = calls
                    differs = pd.Series(False, index=first_calls.index)
                else:
                    differs = differs | (calls.reindex(first_calls.index) != first_calls)
            if first_calls is None:
                return []
            same = ~differs & (first_calls.astype(str).str.len() == 1)
            return list(first_calls.index[same])

        def first_pass():
            for sample in sample_dict.keys(): # update grouping sample_dict with normalized dataframes
                norm_df, calls, emb = normalize_one(sample)
                if norm_df['ALT'].isin(['-', 'N']).all():
                    # No determined base at any of the group's positions ('-' or
                    # N everywhere).  RAxML exits on such a sequence ("sequences
                    # that consist entirely of undetermined values") and the
                    # whole group lost its table and tree.  The sample is left
                    # out of this group before the parsimony test; it is named in
                    # the terminal and the summary and stays in every other
                    # group it belongs to.
                    undetermined.append(sample)
                    continue
                if emb:
                    embedded_added[sample] = emb
                norm_sample_dict[sample] = norm_df
                yield calls

        parsimony_positions = parsimony_uninformative(first_pass())
        retained = df_ref[~df_ref['abs_pos'].isin(parsimony_positions)]

        # A sample whose only determined bases sit at positions the test just
        # removed is still all '-'/N over the retained positions, and RAxML
        # would exit on it just the same.  Leave those out too and repeat the
        # test without them until nothing changes (rare: once per group with
        # such a sample, never otherwise).
        while len(retained) and norm_sample_dict:
            kept = set(retained['abs_pos'])
            late = [s for s, d in norm_sample_dict.items()
                    if d.loc[d['abs_pos'].isin(kept), 'ALT'].isin(['-', 'N']).all()]
            if not late:
                break
            for s in late:
                del norm_sample_dict[s]
                embedded_added.pop(s, None)
                undetermined.append(s)
            parsimony_positions = parsimony_uninformative(
                normalize_one(s, quiet=True)[1] for s in norm_sample_dict)
            retained = df_ref[~df_ref['abs_pos'].isin(parsimony_positions)]

        if undetermined:
            # setdefault: a Group built without __init__ (tests, pickles) has no such attribute.
            self.__dict__.setdefault('undetermined_in_group', defaultdict(list))[group].extend(undetermined)
            print(f'  {group}: {len(undetermined)} sample(s) have no determined base at any '
                  f'of the group\'s informative positions and are left out of this '
                  f'group\'s tree and table: {", ".join(undetermined)}')
            if not norm_sample_dict:
                print(f'  {group}: no sample can be placed. Nothing to align, so no tree; '
                      f'group skipped and the rest of the run continues.')
                return {}, group
        # Every sample was normalized against df_ref by an outer merge, so all of
        # them carry the same positions, and the same positions are removed from
        # each.  They therefore empty out together: one sample having nothing left
        # means the whole group has nothing left.  Reported as a property of the
        # group rather than of a sample, because naming a sample here reads as
        # "that VCF is malformed" when the sample is only the first one iterated.
        # Count the calls read from multi-base records that stay in this group's
        # table.  A position can also leave the table because of them: with the
        # sample's base read correctly, every sample in the group may now carry
        # the same base, and the position is then not informative.
        kept = set(retained['abs_pos'])
        uninformative = set(parsimony_positions)
        for sample, positions in embedded_added.items():
            self.embedded_cells[sample] += len(positions & kept)
            self.embedded_now_uninformative[group] |= positions & uninformative
        if retained.empty:
            # A one-sample group is the common case and deserves to be named as
            # such: with a single sample every position is trivially uninformative,
            # which is a statement about the group's size, not about its calls.
            if len(norm_sample_dict) == 1:
                reason = ('it holds one sample, so no position can distinguish '
                          'samples')
            else:
                reason = (f'every position carries the same call in all '
                          f'{len(norm_sample_dict)} of its samples')
            print(f'  {group}: no parsimony-informative positions -- {reason}. '
                  f'Nothing to align, so no tree; group skipped and the rest of '
                  f'the run continues.')
            return {}, group

        parsmony_sample_dict={}
        for sample, df_norm in norm_sample_dict.items():
            df_norm = df_norm[~df_norm['abs_pos'].isin(parsimony_positions)]
            if df_norm.empty:
                # Columns declared even though there are no rows.  A bare
                # DataFrame() has none at all, and dict_to_fasta then failed with
                # KeyError: 'abs_pos' -- which reads as a missing column, i.e. a
                # structural defect, rather than as an empty result.  With the
                # columns present an empty frame reaches the alignment-length check,
                # which says what is actually wrong.
                parsmony_sample_dict[sample] = pd.DataFrame(columns=['abs_pos', 'ALT'])
            else:
                df_norm= self.sort_df(df_norm)
                parsmony_sample_dict[sample] = df_norm

        # Add root.  The retained positions are stated explicitly here; this used to
        # read df_norm, the loop variable left behind by the loop above, so the root
        # was derived from whichever sample happened to be iterated last.  When that
        # sample's frame was empty the root came out empty too, which produced a
        # FASTA with no 'root' record, and `raxml -o root` then failed with no tree.
        df_root = retained[['abs_pos', 'REF']].rename(columns={"REF": "ALT"}) #fake the column name as ALT so it is as samples when dict_to_fasta function is called on the dataframes
        if not df_root.empty:
            df_root = self.sort_df(df_root)
            parsmony_sample_dict['root'] = df_root
        return parsmony_sample_dict, group
    
    def sort_df(self, df):
        # .copy() because callers pass a boolean-mask slice of a larger frame, and
        # the two assignments below write into it.  Under copy-on-write that wrote
        # to a temporary and was silently discarded; before it, it wrote through to
        # the parent frame.  Either way it emitted a SettingWithCopyWarning per
        # sample per group, which is what filled the run log with hundreds of lines.
        df = df.copy()
        df[['chrom', 'pos']] = df['abs_pos'].str.split(':', expand=True)
        df['pos'] = df['pos'].astype(int)
        df = df.sort_values('pos').reset_index(drop=True)
        df = df[['abs_pos', 'ALT']]
        return df

    def raxml_table_build(self, group):
        tree = Tree(fasta_alignments=self.group_fasta_dict[group], write_path='{}/{}'.format(self.cwd, group), tree_name=group)
        tables = Tables(
            df_alignments=self.group_dataframe_dict[group],
            tree=tree.newick,
            gbk=self.annotation_df,
            mq=self.average_mq_df,
            write_path='{}/{}'.format(self.cwd, group),
            groupings_dict=self.groupings_dict,
            show_groups=self.show_groups,
            table_name=group,
            debug=False,
            sample_coverage_dict=self.sample_coverage_dict if self.dp else {}
        )
        tables.build_tables()
        self.raxml_version = tree.raxml_version
        # Only run html_tree if the option is specified
        if self.html_tree:
            self.html_tree_build(group, tables)

    @staticmethod
    def _html_tree_timeout(tables):
        '''
        Seconds allowed for one group's HTML tree: VSNP3_HTML_TREE_TIMEOUT when
        set, otherwise 900, raised for large groups.  Measured at about 34 s for
        1,000 samples and 124 s for 2,000, roughly quadratic; twice that fit is
        allowed, so the all-VCF tree of a 10,000-sample run (about 50 minutes)
        is not cut off by a constant written for a few hundred samples.
        '''
        try:
            base = int(str(os.environ.get('VSNP3_HTML_TREE_TIMEOUT') or '').strip() or 900)
        except ValueError:
            base = 900
        df = getattr(tables, 'table_to_tree_df', None)
        samples = len(df) if df is not None else 0
        return max(base, int(6.2e-5 * samples * samples))

    def html_tree_build(self, group, tables):
        '''
        The HTML tree for one group, in a child process with a timeout.  A
        failure is recorded in html_tree_failures and the run goes on.
        '''
        try:
            if hasattr(tables, 'table_to_tree_path') and tables.table_to_tree_path:
                # The cascade1 frame goes with the path, so the child does
                # not read the workbook it was just written from back
                # through openpyxl (5.6 s at 400 samples, longer than the
                # tree itself now takes).
                # raise_errors: an exception in html_tree reaches the
                # wrapper with its message.  Without it html_tree prints
                # the error and returns 1, which came back here as a
                # normal result, so the group was never listed as failed.
                timeout = self._html_tree_timeout(tables)
                status = run_with_timeout(html_tree, timeout, tables.table_to_tree_path,
                                          df=getattr(tables, 'table_to_tree_df', None),
                                          raise_errors=True)
                if status not in (0, None):
                    raise RuntimeError(f'html_tree returned {status}; its error is '
                                       f'printed in the terminal output')
            else:
                # Said in the summary too: the tree was silently absent for
                # exactly the largest groups.
                reason = (getattr(tables, 'html_tree_skipped', None)
                          or 'table split into several workbooks; no HTML tree is built from a split table')
                self.html_tree_failures[group] = reason
                print(f'{group} HTML tree not built: {reason}')
        except TimeoutError:
            self.html_tree_failures[group] = (f'timed out after {timeout} s; set VSNP3_HTML_TREE_TIMEOUT '
                                              f'(seconds) to allow more')
            print(f'{group} HTML tree generation timed out after {timeout} s '
                  f'(VSNP3_HTML_TREE_TIMEOUT raises the limit)')
        except RuntimeError as e:
            # The tables and the tree for this group are already written;
            # only the HTML view is missing.  Before run_with_timeout
            # reported this, a crash inside html_tree hung the whole run.
            self.html_tree_failures[group] = str(e).splitlines()[0]
            print('{} HTML tree generation failed: {}'.format(group, e))

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
    parser.add_argument('-w', '--qual_threshold', action='store', dest='qual_threshold', default=QUAL_THRESHOLD, required=False, help='Optional: Minimum QUAL threshold for calling a SNP')
    parser.add_argument('-x', '--n_threshold', action='store', dest='n_threshold', default=N_THRESHOLD, required=False, help='Optional: Minimum N threshold.  SNPs between this and qual_threshold are reported as N')
    parser.add_argument('-y', '--mq_threshold', action='store', dest='mq_threshold', default=MQ_THRESHOLD, required=False, help='Optional: At least one position per group must have this minimum MQ threshold to be called.')
    parser.add_argument('-t', '--reference_type', action='store', dest='reference_type', required=False, help='Reference type group/directory with dependencies')
    parser.add_argument('-n', '--no_filters', action='store_true', dest='no_filters', default=False, help='Optional: turn off filters.  Default matches vsnp3_step2.py, which applies them.')
    parser.add_argument('--filter_density', action='store_true', dest='filter_density', help='Optional: Remove SNPs when density threshold is exceeded within specified window size')
    parser.add_argument('--density_threshold', action='store', dest='density_threshold', type=int, default=3, help='Optional: Minimum number of SNPs required to trigger density filtering (default: 3)')
    parser.add_argument('--density_window', action='store', dest='density_window', type=int, default=20, help='Optional: Window size in base pairs for density filtering (default: 20)')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling')
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.basename(__file__), __version__))
    args = parser.parse_args()

    if args.reference_type:
        ro = Ref_Options(args.reference_type)
        if ro.metadata and not args.metadata:
            args.metadata = ro.metadata
        if ro.defining_snps and not args.defining_snps:
            args.defining_snps = ro.defining_snps
        if ro.gbk and not args.gbk:
            args.gbk = ro.gbk

    group = Group(pickle_file=args.pickle_file, metadata=args.metadata, gbk_list=args.gbk, defining_snps=args.defining_snps, abs_pos=args.abs_pos, group=args.group, no_filters=args.no_filters, qual_threshold=int(args.qual_threshold), n_threshold=int(args.n_threshold), mq_threshold=int(args.mq_threshold), filter_density=args.filter_density, density_threshold=args.density_threshold, density_window=args.density_window, debug=args.debug)

# Created 2021 by Tod Stuber