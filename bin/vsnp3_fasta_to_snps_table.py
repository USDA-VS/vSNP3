#!/usr/bin/env python

from vsnp3_version import __version__

import os
import subprocess
import shutil
import sys
import re
import glob
import time
import random
import locale
from datetime import datetime, date, timedelta
from io import StringIO
import pandas as pd
import numpy as np
import multiprocessing
import threading
from concurrent import futures
import argparse
import textwrap
from Bio import SeqIO
from Bio import Phylo
from collections import defaultdict

import vsnp3_run

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")

# ---------------------------------------------------------------------------
# Bootstrap support values
#
# By default this script runs RAxML once and keeps the best tree only
# (RAxML_bestTree.raxml), so the trees it writes carry NO bootstrap support
# values. That is upstream vSNP3 behaviour, and it is the default here too.
#
# To get support values, set VSNP3_BOOTSTRAP to the number of replicates in
# the environment before running Step 2 — there is no command-line flag:
#
#     VSNP3_BOOTSTRAP=100 vsnp3_step2.py <usual arguments>
#
# or, for a whole shell session:
#
#     export VSNP3_BOOTSTRAP=100
#     vsnp3_step2.py <usual arguments>
#
# With N > 0, RAxML runs its rapid-bootstrap analysis (-f a -x 7777 -N N) and
# the tree kept is RAxML_bipartitions.raxml — the best-scoring ML tree of that
# analysis, with bootstrap proportions as internal node labels. Unset, 0, a
# negative or a non-numeric value all mean stock best-tree-only behaviour.
#
# The variable is read by bootstrap_replicates() below, from Tree.__init__ (the
# KL_BOOTSTRAP block) and from the step 2 summary, which states the tree mode
# and the value set. It is an environment variable rather than a flag so that
# vsnp3_group_on_defining_snps.py needs no changes — the value crosses the
# import chain and RAxML worker pools for free. Note it applies to EVERY tree
# the run builds (each group, plus the all-VCF tree with -a), and the
# bootstrap analysis is the expensive part: expect a much longer tree stage
# on large runs.
#
# That tree is not the default tree with labels added: -f a searches by a
# different procedure, so the topology, and with it the table row order, can
# differ from the best-tree-only run (on a 40-sample test set, 4 of 38 splits).
# The 3.37 validation covered the default only; leave VSNP3_BOOTSTRAP unset for
# validated runs until this mode has been validated in its own right.
# ---------------------------------------------------------------------------

# RAxML build names, most specific first.  Order is preserved from the original
# probe so the same binary is chosen on any given machine.
_RAXML_CANDIDATES = ('raxml', 'raxmlHPC-PTHREADS-AVX2', 'raxmlHPC-PTHREADS-AVX',
                     'raxmlHPC-PTHREADS', 'raxmlHPC-SSE3', 'raxmlHPC')
_RAXML = None


def _resolve_raxml():
    '''
    Name of the RAxML binary to use, resolved once per process.

    This replaces a six-deep try/except chain that discovered the binary by
    *executing* each candidate with no arguments and catching OSError.  Two
    problems with that: it launched up to seven processes on every Tree()
    construction -- once per group, so ~1,400 launches on a 197-group run -- and
    the successful candidate printed its usage text to stderr, which was not
    redirected, so RAxML usage appeared in the middle of step 2's output.

    shutil.which() answers the same question from PATH without running anything.
    '''
    global _RAXML
    if _RAXML is None:
        _RAXML = next((name for name in _RAXML_CANDIDATES if shutil.which(name)),
                      'raxmlHPC')
    return _RAXML


def _env_int(name, default):
    try:
        value = int(str(os.environ.get(name) or '').strip() or default)
    except ValueError:
        value = default
    return value


def bootstrap_replicates():
    '''
    (replicates, raw value) for VSNP3_BOOTSTRAP; see "Bootstrap support values"
    at the top of this file.  Unset, 0, negative or not a whole number all mean
    0 replicates, the best-tree-only default.  The raw value is returned too so
    the step 2 summary can show what was set, not only what it was read as.
    '''
    raw = os.environ.get('VSNP3_BOOTSTRAP')
    try:
        replicates = int(str(raw or '0').strip() or '0')
    except ValueError:
        replicates = 0
    return max(0, replicates), raw


def raxml_threads(binary=None):
    '''
    Threads one RAxML run is given (-T).  4 unless VSNP3_RAXML_THREADS says
    otherwise, and 1 for a build that is not PTHREADS, which ignores -T.
    '''
    binary = binary or _resolve_raxml()
    name = os.path.basename(binary)
    if 'PTHREADS' not in name and name != 'raxml':
        return 1
    return max(1, _env_int('VSNP3_RAXML_THREADS', 4))


_RAXML_SLOTS = None
_RAXML_SLOTS_LOCK = threading.Lock()
_RAXML_VERSION = None


def raxml_slots():
    '''
    Semaphore bounding how many RAxML runs are in flight at once.

    Step 2 builds every group's tree from a pool of cpu_count/1.2 threads,
    and each run asked for 4 threads, so a 16-core machine ran 13 RAxML
    processes x 4 = 52 threads at once.  RAxML's worker threads busy-wait,
    so oversubscription is not merely wasteful: spinning threads take cores
    from the ones with work.  The thread count per run is unchanged -- the
    same binary with the same -T gives the same tree -- but at most
    cores / threads runs are started at once (VSNP3_MAX_CPUS caps the cores
    counted).  The pool itself is not shrunk: a group whose RAxML has
    finished carries on with its tables and HTML tree while others wait.
    '''
    global _RAXML_SLOTS
    with _RAXML_SLOTS_LOCK:
        if _RAXML_SLOTS is None:
            cores = multiprocessing.cpu_count()
            cap = _env_int('VSNP3_MAX_CPUS', 0)
            if cap > 0:
                cores = min(cores, cap)
            _RAXML_SLOTS = threading.BoundedSemaphore(max(1, cores // raxml_threads()))
        return _RAXML_SLOTS


def _raxml_version(raxml):
    '''The banner line, read once per process rather than once per group.'''
    global _RAXML_VERSION
    if _RAXML_VERSION is None:
        try:
            banner = subprocess.check_output([raxml, '-v'], text=True).split('\n')
            _RAXML_VERSION = next((line for line in banner if 'RAxML' in line),
                                  banner[0] if banner else 'RAxML version unknown')
        except (subprocess.CalledProcessError, OSError, IndexError) as e:
            _RAXML_VERSION = f'RAxML version unknown ({e})'
    return _RAXML_VERSION


class Tree:
    '''
    '''

    def __init__(self, fasta_alignments=None, write_path=None, tree_name=None, debug=False,):
        raxml = _resolve_raxml()
        self.raxml = raxml
        # Concurrency is governed by _raxml_slots() (VSNP3_MAX_CPUS), not here.
        self.debug = debug
        
        self.beginTime = datetime.now()
        self.startTime = datetime.now()
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')

        if not write_path:
            write_path = os.getcwd()
        if not tree_name:
            tree_name = "raxml"
        self.write_path = write_path
        self.tree_name = tree_name

        self.raxml_version = _raxml_version(raxml)

        # Checked, and stderr kept.  This was os.system(... > /dev/null 2>&1) with the
        # return value discarded, so the only sign of failure was a FileNotFoundError
        # on the rename below -- which then fell through and left self.newick pointing
        # at a path that does not exist, to fail confusingly much later inside a
        # process pool.
        newick = os.path.join(write_path, '{}_{}.tre'.format(tree_name, self.st))
        # --- BEGIN KL_BOOTSTRAP ---
        # Bootstrap replicates, from the environment (0/unset = best tree only, the
        # upstream behaviour). vsnp_gui's Step 2 form sets this.
        _kl_boot, _ = bootstrap_replicates()
        # -T is still passed to a non-PTHREADS build, as it always was; that
        # build ignores it.  The slot count accounts for that (raxml_threads).
        threads = str(max(1, _env_int('VSNP3_RAXML_THREADS', 4)))
        if _kl_boot > 0:
            command = [raxml, '-f', 'a', '-s', fasta_alignments, '-n', 'raxml', '-m', 'GTRCATI',
                       '-o', 'root', '-w', write_path, '-p', '456123', '-x', '7777', '-N', str(_kl_boot), '-T', threads]
            best_tree = os.path.join(write_path, 'RAxML_bipartitions.raxml')
        else:
            command = [raxml, '-s', fasta_alignments, '-n', 'raxml', '-m', 'GTRCATI',
                       '-o', 'root', '-w', write_path, '-p', '456123', '-T', threads]
            best_tree = os.path.join(write_path, 'RAxML_bestTree.raxml')
        with raxml_slots():
            vsnp3_run.run(command)
        # --- END KL_BOOTSTRAP ---
        vsnp3_run.require_output(best_tree, what='RAxML best tree')
        os.rename(best_tree, newick)

        # RAxML_info is the only account RAxML gives of what it did; it used to be
        # deleted along with everything else, including on the failure path.
        info = os.path.join(write_path, 'RAxML_info.raxml')
        if os.path.exists(info):
            if self.debug:
                os.rename(info, os.path.join(write_path, f'RAxML_info_{tree_name}.txt'))
            else:
                os.remove(info)
        for each in sorted(glob.glob(os.path.join(write_path, 'RAxML*'))):
            os.remove(each)
        for reduced in sorted(glob.glob(os.path.join(write_path, '*.reduced'))):
            os.remove(reduced)
        self.newick = newick
    
    def checksum_match_to_text(self, tree):
        # read entire tree into variable as string obj
        with open(tree, 'rt') as open_tree:
            entire_file = open_tree.read()
            print (entire_file)
        with open("idtable.txt" , 'rt') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split("\t")
                org = str(line[1])
                org = org.rstrip()
                check = str(line[0])
                check = check.rstrip()
                print (org)
                print (check)
                entire_file = re.sub(check, "'" + org + "'", entire_file)
            f.close()
        outfile = "NAMES-UPDATED-" + os.path.basename(tree)
        write_out = open(outfile , 'wt')
        write_out.write(entire_file)
        write_out.close()
        return outfile

    def df_to_fasta(self, df, alignment_file):
        test_duplicates=[] # if duplicate name in alignment fasta raxml with error and exit
        with open(alignment_file, 'w') as write_out:
            for index, row in df.iterrows():
                test_duplicates.append(row.name)
                if test_duplicates.count(row.name) < 2:
                    print('>{}'.format(row.name), file=write_out)
                    for pos in row:
                        print(pos, end='', file=write_out)
                    print("", file=write_out)

    def get_parsimonious_pos(self, in_df):
        try:
            ref_series = in_df.loc['reference_seq']
            in_df = in_df.drop(['reference_seq']) #in all_vcf reference_seq needs to be removed
        except KeyError:
            print('Check that there is a "reference_seq" named')
            sys.exit(1)
        # print('in_df size: {}'.format(in_df.shape))
        parsimony = in_df.loc[:, (in_df != in_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = in_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        return out_df


class Tables:

    def __init__(self, fasta_alignments=None, df_alignments=None, tree=None, gbk=None, mq=None, write_path=None, groupings_dict=None, show_groups=False, table_name=None, debug=False, sample_coverage_dict=None):
        self.fasta_alignments = fasta_alignments
        self.df_alignments = df_alignments
        self.tree = tree
        self.gbk = gbk
        self.mq = mq
        self.show_groups = show_groups
        self.debug = debug
        self.st = datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H-%M-%S')
        self.groupings_dict = groupings_dict
        self.sample_coverage_dict = sample_coverage_dict or {}
                
        if self.groupings_dict:
            group_vcfs_dict = defaultdict(list) #invert the key, values
            for group, dataframes in self.groupings_dict.items():
                for vcf in dataframes.keys():
                    group_vcfs_dict[vcf].append(group)
            self.group_vcfs_dict = dict(sorted(group_vcfs_dict.items())) #sorts on key(vcf name)

        if not write_path:
            write_path = os.getcwd()
        if not table_name:
            table_name = "generic"
        self.write_path = write_path
        self.table_name = table_name

    @staticmethod
    def leaf_order(tree_path, table_name=''):
        '''
        'root' followed by the tree's other leaves in file order, read with
        Bio.Phylo.  The table rows follow this order.
        '''
        try:
            newick = Phylo.read(tree_path, 'newick')
        except Exception as e:                                  # noqa: BLE001
            raise RuntimeError(f'[{table_name}] could not read the RAxML tree '
                               f'{os.path.basename(str(tree_path))}: {type(e).__name__}: {e}') from e
        leaves = [clade.name for clade in newick.get_terminals() if clade.name]
        return ['root'] + [name for name in leaves if name != 'root']

    def build_tables(self, ):
        if self.df_alignments is None:
            reformated = os.path.join(self.write_path, 'reformated.fasta')
            with open(reformated, 'w') as reformat:
                sequence = SeqIO.parse(self.fasta_alignments, "fasta")
                for each in sequence:
                    print('>{}\n{}'.format(each.description, each.seq), file=reformat)
            fasta_df = pd.read_csv(reformated, header=None, sep='^')
            os.remove(reformated)
            seq = fasta_df.iloc[1::2].reset_index(drop=True)
            header = fasta_df.iloc[0::2].reset_index(drop=True).replace(to_replace=r'>', value='', regex=True)
            seq = seq.rename({0:"seq"}, axis='columns')
            header = header.rename({0:"header"}, axis='columns')
            fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
            fasta_df = fasta_df.rename({0: 'header', 1: 'seq'}, axis='columns')
            seq = fasta_df['seq'].apply(lambda x: pd.Series(list(x)))
            fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
            fasta_df = fasta_df.set_index(0)
        else:
            fasta_df = self.df_alignments
        
        # Clean Unicode characters in the FASTA DataFrame index to match the main script processing
        import unicodedata
        
        def clean_unicode(text):
            # Normalize Unicode and convert to ASCII
            normalized = unicodedata.normalize('NFKD', str(text))
            ascii_text = normalized.encode('ascii', 'ignore').decode('ascii')
            return ascii_text
        
        # Apply the same Unicode cleaning to FASTA sample names that was applied in main script
        fasta_df.index = fasta_df.index.map(clean_unicode)
        
        # Leaf order from the tree, read with Bio.Phylo.  This was a chain of
        # regex substitutions over the Newick text; the one meant to strip
        # branch lengths, r'[0-9].*\.[0-9].*\n', matched anywhere in a line, so
        # a sample name containing <digit>.<digit> was deleted (22-12345.1_zc.vcf
        # vanished) or cut and glued to the next leaf (NZ_CP0123.1.vcf and
        # B_zc.vcf became NZ_CPB_zc.vcf).  The filter below then kept only names
        # present in the alignment, so those samples left every table and the
        # HTML tree with no message.  get_terminals() walks the tree in file
        # order, so the rows come out in the order the regex gave for names it
        # did not damage; underscores and quotes are handled by the parser.
        sample_order = self.leaf_order(self.tree, self.table_name)
        
        # Clean Unicode characters in sample_order to match the processing done in main script
        cleaned_sample_order = []
        for sample in sample_order:
            if sample == 'root':
                cleaned_sample_order.append(sample)
            else:
                cleaned_sample_order.append(clean_unicode(sample))
        
        sample_order = cleaned_sample_order
        
        # DEBUG: Print available samples and requested samples
        if self.debug:
            print("[{}] Available samples in fasta_df:".format(self.table_name))
            print(sorted(fasta_df.index.tolist()))
            print("[{}] Requested samples in sample_order:".format(self.table_name))
            print(sorted(sample_order))
        
        # Filter sample_order to only include samples that exist in fasta_df
        available_samples = set(fasta_df.index)
        filtered_sample_order = [sample for sample in sample_order if sample in available_samples]
        
        # A sample in the tree that is absent from the alignment is an error, not
        # something to guess around.  This used to substitute possible_matches[0],
        # chosen by case-insensitive substring match -- so '19-01' matched '19-011'
        # and the wrong isolate's SNP row was written into the table.  The near-miss
        # diagnostic is genuinely useful for tracking down name mismatches, so it is
        # kept; only the substitution is gone.
        missing_samples = sorted(set(sample_order) - available_samples)
        if missing_samples:
            detail = []
            for missing_sample in missing_samples:
                near = sorted(s for s in available_samples
                              if missing_sample.lower() in s.lower()
                              or s.lower() in missing_sample.lower())
                detail.append(f'  {missing_sample!r}'
                              + (f' -- similar names present: {near}' if near else
                                 ' -- no similar name present'))
            raise RuntimeError(
                f'[{self.table_name}] {len(missing_samples)} sample(s) in the tree are '
                f'not in the alignment:\n' + '\n'.join(detail)
                + '\nThis is usually a name mismatch between the metadata file and the '
                  'VCF file names, or a Unicode difference. Refusing to substitute a '
                  'similarly named sample, because that would put one isolate\'s SNP '
                  'calls on another isolate\'s row.')
        # The other direction.  A sample in the alignment but absent from the
        # tree would otherwise be silently left out of every table, because
        # the filter above keeps only tree names.  RAxML writes every taxon it
        # was given, so this can only be a parsing fault, and it must be loud.
        not_in_tree = sorted(available_samples - set(sample_order))
        if not_in_tree:
            raise RuntimeError(
                f'[{self.table_name}] {len(not_in_tree)} sample(s) in the alignment are '
                f'not in the tree {os.path.basename(self.tree)}: {not_in_tree[:10]}'
                + (' ...' if len(not_in_tree) > 10 else '')
                + '. No table was written for this group rather than one with rows missing.')

        # Remove duplicates while preserving order
        seen = set()
        final_sample_order = []
        for sample in filtered_sample_order:
            if sample not in seen:
                final_sample_order.append(sample)
                seen.add(sample)
        
        if self.debug:
            print("[{}] Final sample order: {}".format(self.table_name, final_sample_order))
        
        # Use the filtered sample order
        tree_order = fasta_df.loc[final_sample_order]
        tree_order2 = fasta_df.loc[final_sample_order]
        
        # Continue with the rest of the existing method...
        self.write_out_table(tree_order, 'sorted')

        # The three per-column statistics the cascade sorts use, from one pass
        # over the array.  They were three nested Python loops over every cell,
        # each re-indexing column.iloc[0] per element: 3.5 s at 100 samples x
        # 2,000 positions and 28 s at 400 x 4,000, per group.
        snp_per_column, snp_from_top, longest_run_start = self.snp_column_stats(tree_order)

        ## Sort bias to total number of SNPs per column
        # count number of SNPs in each column
        row1 = pd.Series(snp_per_column, tree_order.columns, name="snp_per_column")

        # get the snp count per column: rows below the root before the first SNP
        row2 = pd.Series(snp_from_top, tree_order.columns, name="snp_from_top")
        tree_order = pd.concat([tree_order, pd.DataFrame([row1])])
        tree_order = pd.concat([tree_order, pd.DataFrame([row2])])
        tree_order = tree_order.T
        tree_order = tree_order.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order = tree_order.T

        # remove snp_per_column and snp_from_top rows
        cascade_order_df = tree_order[:-2]
        self.write_out_table(cascade_order_df, 'cascade1')
        del row1
        del row2
        del snp_from_top
        del cascade_order_df

        ## Start 2nd cascading table: sort bias to longest continues vertical SNP count per column
        row1 = pd.Series(snp_per_column, tree_order2.columns, name="snp_per_column")
        # where the longest unbroken vertical run of SNPs starts, per column
        row2 = pd.Series(longest_run_start, tree_order2.columns, name="snp_from_top")
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row1])])
        tree_order2 = pd.concat([tree_order2, pd.DataFrame([row2])])
        tree_order2 = tree_order2.T
        tree_order2 = tree_order2.sort_values(['snp_from_top', 'snp_per_column'], ascending=[True, False])
        tree_order2 = tree_order2.T

        # remove snp_per_column and snp_from_top rows
        cascade_order_df = tree_order2[:-2]
        self.write_out_table(cascade_order_df, 'cascade2')


    @staticmethod
    def snp_column_stats(frame):
        '''
        Per-column statistics of an alignment frame whose first row is the root.

        A cell is a SNP when it differs from the root call and is not '-'.
        Returns three lists, one value per column:
          snp_per_column    SNP cells in the column, root row included
          snp_from_top      rows below the root before the first SNP cell
                            (every row when there is none)
          longest_run_start 0-based row below the root where the longest
                            unbroken run of SNP cells begins; the first such
                            run on a tie; 0 when there is none
        Comparisons are made on the cell objects, as the loops made them.
        '''
        values = frame.to_numpy(dtype=object)
        n_rows, n_cols = values.shape
        if n_rows == 0 or n_cols == 0:
            return [0] * n_cols, [max(n_rows - 1, 0)] * n_cols, [0] * n_cols
        root = values[0]
        diff = (values != root[None, :]) & (values != '-')
        snp_per_column = diff.sum(axis=0)

        below = diff[1:]
        n_below = below.shape[0]
        if n_below:
            snp_from_top = np.where(below.any(axis=0), below.argmax(axis=0), n_below)
        else:
            snp_from_top = np.zeros(n_cols, dtype=np.int64)

        longest = np.zeros(n_cols, dtype=np.int64)
        if n_below:
            padded = np.zeros((n_below + 2, n_cols), dtype=bool)
            padded[1:-1] = below
            # run r starts at row s (below coordinates) and ends before row e
            s_cols, s_rows = np.nonzero((~padded[:-1] & padded[1:]).T)
            e_cols, e_rows = np.nonzero((padded[:-1] & ~padded[1:]).T)
            if len(s_rows):
                lengths = e_rows - s_rows
                colmax = np.zeros(n_cols, dtype=np.int64)
                np.maximum.at(colmax, s_cols, lengths)
                is_max = lengths == colmax[s_cols]
                cols_m, first = np.unique(s_cols[is_max], return_index=True)
                longest[cols_m] = s_rows[is_max][first]
        return snp_per_column.tolist(), snp_from_top.tolist(), longest.tolist()

    ###Break and write out table
    # In Tables class, modify write_out_table method:
    def write_out_table(self, df, table_type=None):
        if hasattr(self.mq, 'abs_pos'):
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.mq)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T

        if hasattr(self.gbk, 'abs_pos') and not self.gbk.empty:
            df_temp = df.T.reset_index()
            df_temp = df_temp.merge(self.gbk)
            df_temp = df_temp.set_index('abs_pos')
            df = df_temp.T
        else:
            df = pd.concat([df, pd.Series(name='no annotations').to_frame().T])

        # Add coverage data (only once)
        if self.sample_coverage_dict:
            coverage_data = {
                sample: self.sample_coverage_dict.get(sample, 0.0)
                for sample in df.index 
                if sample not in ['root', 'no annotations', 'Average Coverage']
            }
            if coverage_data:  # Only add if we have coverage data
                df['Average Coverage'] = pd.Series(coverage_data)

        if self.show_groups:
            # Add grouping information
            joined_data = {key: '; '.join(map(str, value)) for key, value in self.group_vcfs_dict.items()}
            new_series = pd.Series(joined_data)
            if 'Grouping' not in df.columns:
                df.insert(0, 'Grouping', new_series)

        # Continue with existing chunking code...
        max_size=10000
        count=0
        chunk_start=0
        chunck_end=0
        column_count = df.shape[1]
        if column_count > max_size:
            while column_count > max_size:
                count += 1
                chunck_end += max_size
                df_split = df.iloc[:, chunk_start:chunck_end].copy()
                # Test df_split, not df: 'Grouping' was inserted into df above, so
                # this guard was always False and only the first chunk -- the one
                # whose slice happens to include column 0 -- carried the column,
                # while excel_formatter styled every chunk as though it did.
                if self.show_groups and 'Grouping' not in df_split.columns:
                    df_split.insert(0, 'Grouping', new_series)
                self.excel_formatter(self.json_round_trip(df_split),
                                os.path.join(self.write_path, '{}_{}table{}-{}.xlsx'.format(self.table_name, table_type + '_', count, self.st)))
                chunk_start += max_size
                column_count -= max_size
            count += 1
            df_split = df.iloc[:, chunk_start:].copy()
            if self.show_groups and 'Grouping' not in df_split.columns:
                df_split.insert(0, 'Grouping', new_series)
            self.excel_formatter(self.json_round_trip(df_split),
                            os.path.join(self.write_path, '{}_{}table{}-{}.xlsx'.format(self.table_name, table_type + '_', count, self.st)))
            self.multiple_excel_files = True
            if table_type == 'cascade1':
                # No HTML tree is built from a split table.  Said here so the
                # step 2 summary can list the group, instead of the tree being
                # silently absent for exactly the largest groups.
                self.html_tree_skipped = (f'table split into {count} workbooks '
                                          f'({df.shape[1]:,} columns, more than {max_size:,} '
                                          f'positions); no HTML tree is built from a split table')
        else:
            if self.show_groups and 'Grouping' not in df.columns:
                df.insert(0, 'Grouping', new_series)
            table_df = self.json_round_trip(df)
            write_to = os.path.join(self.write_path, '{}_{}_table-{}.xlsx'.format(self.table_name, table_type, self.st))
            self.excel_formatter(table_df, write_to)
            self.multiple_excel_files = False
            if table_type == 'cascade1':
                self.table_to_tree_path = write_to
                # The frame html_tree would get by reading the workbook back,
                # so step 2 can hand it over instead of reading the workbook.
                self.table_to_tree_df = self.as_read_back(table_df)

    # The JSON round trip every table went through on its way to Excel.  It is
    # kept because it is not a no-op: floats are rounded to 10 significant
    # digits, an all-numeric index becomes integers, and both reach the
    # workbook.  It just no longer goes through a file in the output directory.
    @staticmethod
    def json_round_trip(df):
        return pd.read_json(StringIO(df.to_json(orient='split')), orient='split')

    # What pandas.read_excel reads as missing by default.  Mirrored in
    # as_read_back so a frame handed to html_tree matches one read from disk.
    EXCEL_NA_VALUES = frozenset([
        '-1.#IND', '1.#QNAN', '1.#IND', '-1.#QNAN', '#N/A N/A', '#N/A', 'N/A',
        'n/a', 'NA', '<NA>', '#NA', 'NULL', 'null', 'NaN', '-NaN', 'nan', '-nan',
        'None', ''])

    def as_read_back(self, table_df):
        '''
        table_df as vsnp3_html_tree.read_table returns the workbook that
        excel_formatter writes from it: the index as a first column named
        'Unnamed: 0', the default NA strings read as missing everywhere except
        in that first column, and the Average Coverage column as
        excel_formatter rewrites it (rounded to an integer, blank on the root,
        MQ and 'no annotations' rows).

        The first column holds sample names, and 'NA', 'None', 'null' or 'nan'
        is a name there, not a missing value.  Read as missing, such a sample
        was labelled "nan" in the HTML tree and two of them shared one clade
        key, which showed as support above 100%.
        '''
        out = table_df.reset_index()
        out.columns = ['Unnamed: 0'] + list(table_df.columns)
        values = out.to_numpy(dtype=object)
        na = pd.Series(values.ravel()).isin(list(self.EXCEL_NA_VALUES)).to_numpy().reshape(values.shape)
        na[:, 0] = False
        if na.any():
            for col in np.flatnonzero(na.any(axis=0)):
                column = out.iloc[:, col].astype(object)
                column[na[:, col]] = np.nan
                out.iloc[:, col] = column
        if 'Average Coverage' in out.columns:
            coverage = pd.to_numeric(out['Average Coverage'], errors='coerce').astype(float)
            coverage[coverage.notna()] = np.rint(coverage[coverage.notna()])
            blank = out['Unnamed: 0'].isin(['root', 'MQ', 'no annotations'])
            coverage[blank] = np.nan
            out['Average Coverage'] = coverage
        return out

    @staticmethod
    def _write_cell(ws, row, col, value, cell_format=None):
        '''
        One cell as pandas' xlsxwriter engine wrote it: missing values are
        left blank, integers and floats are numbers, booleans are booleans,
        anything else is its str().  Strings go through ws.write() only when
        xlsxwriter would treat them specially (formula or URL forms).

        cell_format is for the index column, which pandas wrote with its
        header style: a missing index value is then a formatted blank cell.
        '''
        if value is None or (isinstance(value, float) and value != value):
            if cell_format is not None:
                ws.write_blank(row, col, None, cell_format)
            return
        if isinstance(value, str):
            if value == '':
                if cell_format is not None:
                    ws.write_blank(row, col, None, cell_format)
                return
            if value.startswith('=') or value.startswith('{='):
                # Never a formula: a sample or metadata name that starts with
                # '=' went through ws.write(), which made it a formula cell
                # that Excel then "repaired".
                ws.write_string(row, col, value, cell_format)
            elif ':' in value:
                ws.write(row, col, value, cell_format)
            else:
                ws.write_string(row, col, value, cell_format)
        elif isinstance(value, (bool, np.bool_)):
            ws.write_boolean(row, col, bool(value), cell_format)
        elif isinstance(value, (int, np.integer)):
            ws.write_number(row, col, int(value), cell_format)
        elif isinstance(value, (float, np.floating)):
            value = float(value)
            if value != value:
                return
            if value == float('inf'):
                ws.write_string(row, col, 'inf', cell_format)
            elif value == float('-inf'):
                ws.write_string(row, col, '-inf', cell_format)
            else:
                ws.write_number(row, col, value, cell_format)
        elif isinstance(value, (datetime, date, timedelta)):
            ws.write(row, col, value, cell_format)
        else:
            try:
                if pd.isna(value):
                    if cell_format is not None:
                        ws.write_blank(row, col, None, cell_format)
                    return
            except (TypeError, ValueError):
                pass
            ws.write_string(row, col, str(value), cell_format)

    def excel_formatter(self, df_json, write_to, group=None):
        """Format Excel output with proper styling and conditional formatting."""
        import xlsxwriter
        st = self.st
        if isinstance(df_json, pd.DataFrame):
            table_df = df_json
        else:
            table_df = pd.read_json(df_json, orient='split')
        wb = xlsxwriter.Workbook(write_to)
        wb.use_zip64()
        ws = wb.add_worksheet('Sheet1')

        # Define formats
        formatA = wb.add_format({'bg_color': '#58FA82'})
        formatG = wb.add_format({'bg_color': '#F7FE2E'})
        formatC = wb.add_format({'bg_color': '#0000FF'})
        formatT = wb.add_format({'bg_color': '#FF0000'})
        formatnormal = wb.add_format({'bg_color': '#FDFEFE'})
        formatlowqual = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formathighqual = wb.add_format({'font_color': '#000000', 'bg_color': '#FDFEFE'})
        formatambigous = wb.add_format({'font_color': '#C70039', 'bg_color': '#E2CFDD'})
        formatN = wb.add_format({'bg_color': '#E2CFDD'})
        format_rotation = wb.add_format({})
        format_rotation.set_rotation(90)
        formatannotation = wb.add_format({'font_color': '#0A028C', 'rotation': '-90', 'align': 'top'})
        # Column A (sample names, root, MQ and the annotation label) as
        # to_excel wrote the index: pandas' header style, bold with a thin
        # border, centred, top-aligned.  A cell format takes precedence over a
        # row format, so the annotation label keeps this rather than taking the
        # annotation row's rotated blue.
        formatindex = wb.add_format({'bold': True, 'border': 1, 'align': 'center', 'valign': 'top'})

        # Format for integer display in coverage column
        formatcoverage = wb.add_format({
            'num_format': '0',  # Integer format
            'align': 'right'
        })

        rows, cols = table_df.shape
        end_col = cols

        # Set starting column based on whether groups are shown
        if self.show_groups:
            start_col = 2
        else:
            start_col = 1

        # Basic column formatting
        ws.set_column(0, 0, 30)  # Sample name column
        ws.set_column(1, cols, 2.1)  # Position columns

        # Get Coverage column position if present
        coverage_col = None
        if 'Average Coverage' in table_df.columns:
            coverage_col = table_df.columns.get_loc('Average Coverage') + 1
            ws.set_column(coverage_col, coverage_col, 4.2)  # Double width for coverage column

        # Freeze panes
        ws.freeze_panes(2, start_col)

        # Last column for conditional formatting: stop before coverage column, or include all SNP columns
        if coverage_col is not None:
            format_end_col = coverage_col - 1
        else:
            format_end_col = end_col

        # Apply quality-based formatting (excluding Average Coverage column)
        ws.conditional_format(rows - 2, start_col, rows - 1, format_end_col,
                            {'type': 'cell', 'criteria': '<', 'value': 55, 'format': formatlowqual})

        # Apply reference-based formatting (excluding Average Coverage column)
        if self.show_groups:
            ws.conditional_format(2, start_col, rows - 2, format_end_col,
                                {'type': 'cell', 'criteria': '==', 'value': 'C$2', 'format': formatnormal})
        else:
            ws.conditional_format(2, start_col, rows - 2, format_end_col,
                                {'type': 'cell', 'criteria': '==', 'value': 'B$2', 'format': formatnormal})

        # Apply nucleotide-based formatting
        nucleotides = [
            ('A', formatA), ('G', formatG), ('C', formatC), ('T', formatT),
            ('S', formatambigous), ('Y', formatambigous), ('R', formatambigous),
            ('W', formatambigous), ('K', formatambigous), ('M', formatambigous),
            ('N', formatN), ('-', formatN)
        ]

        for nuc, format_obj in nucleotides:
            ws.conditional_format(2, start_col, rows - 2, format_end_col,
                                {'type': 'text',
                                'criteria': 'containing',
                                'value': nuc,
                                'format': format_obj})

        # Write column headers with rotation, the index column and the cells.
        # pandas' to_excel did the index and cells, then the headers were
        # rewritten rotated; it spent most of its time building a cell object
        # per value.  Same cells, written directly: 7.7 s -> 3.4 s per 400 x
        # 4,000 table.
        for col_num, col_name in enumerate(table_df.columns):
            ws.write(0, col_num + 1, col_name, format_rotation)
        write_cell = self._write_cell
        body = table_df.to_numpy(dtype=object)
        for row_num, (index_value, row_values) in enumerate(zip(table_df.index, body), 1):
            write_cell(ws, row_num, 0, index_value, formatindex)
            for col_num, value in enumerate(row_values.tolist(), 1):
                write_cell(ws, row_num, col_num, value)

        # Handle Average Coverage column
        coverage_col_name = 'Average Coverage'
        if coverage_col_name in table_df.columns:
            coverage_col = table_df.columns.get_loc(coverage_col_name) + 1
            coverage_data = table_df[coverage_col_name]

            # First, write blank cells for all special rows and the annotation row
            for row_idx in range(1, rows + 2):  # Include both regular rows and annotation row
                try:
                    if row_idx == rows + 1:  # This is the annotation row
                        ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))
                    else:
                        row_name = table_df.index[row_idx-1]
                        if row_name in ['root', 'MQ', 'no annotations']:
                            ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))
                except IndexError:
                    ws.write_blank(row_idx, coverage_col, None, wb.add_format({'border': 0}))

            # Then write coverage values for regular rows
            for row_idx in range(1, rows):
                row_name = table_df.index[row_idx-1]
                if row_name not in ['root', 'MQ', 'no annotations']:
                    try:
                        value = coverage_data.iloc[row_idx-1]
                        if pd.notna(value):
                            # Convert to integer and write
                            int_value = int(round(float(value)))
                            ws.write(row_idx, coverage_col, int_value, formatcoverage)
                        else:
                            ws.write(row_idx, coverage_col, 'N/A', formatcoverage)
                    except (IndexError, ValueError, TypeError):
                        ws.write(row_idx, coverage_col, 'N/A', formatcoverage)

        # Set annotation row
        ws.set_row(rows, 400, formatannotation)
        
        wb.close()

class Hash_Names:

    def __init__(self, fasta_alignments=None, debug=False,):
        self.fasta_alignments = fasta_alignments
        self.debug = debug

    def hash_fasta(self,):
        unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
        hashed_fasta = self.fasta_alignments.replace('.fasta', '_hashed.fasta')
        self.hashed_fasta = hashed_fasta
        checksum_dict = {}
        record_iterator = SeqIO.parse(self.fasta_alignments, "fasta")
        outfasta = open(hashed_fasta , 'at')
        for fasta_file in record_iterator:
            if fasta_file.description == "root":
                print('>{}\n{}'.format(fasta_file.description, fasta_file.seq), file=outfasta)
                checksum_dict.update({fasta_file.description:fasta_file.description})
            else:
                unique_number = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(16))
                if fasta_file.description in checksum_dict.values():
                    dup_header = ''.join(random.choice('0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ') for i in range(2))
                    checksum_dict.update({unique_number:fasta_file.description + "-DUPLICATE_HEADER_NAME-" + dup_header})
                else:
                    checksum_dict.update({unique_number:fasta_file.description})
                print ('>{}\n{}'.format(unique_number, fasta_file.seq), file=outfasta)
        outfasta.close()

        idtable = open("idtable.txt" , 'wt')
        for key, value in checksum_dict.items():
            print('{}\t{}'.format(key, value), file=idtable)
        idtable.close()
        self.idtable = "idtable.txt"
        return hashed_fasta

    def hash_tree_to_original_name(self, tree):
        # read entire tree into variable as string obj
        with open(tree, 'rt') as open_tree:
            entire_file = open_tree.read()
        with open(self.idtable , 'rt') as f:
            for line in f:
                line = line.strip('\n')
                line = line.split("\t")
                org = str(line[1])
                org = org.rstrip()
                check = str(line[0])
                check = check.rstrip()
                entire_file = re.sub(check, "'" + org + "'", entire_file)
            f.close()
        original_tree_names = tree.replace('.tre', '_original_names.tre')
        write_out = open(original_tree_names , 'wt')
        write_out.write(entire_file)
        write_out.close()
        return original_tree_names


class Parsimonious:
    def __init__(self, fasta_alignments=None, debug=False,):

        with open('reformated.fasta', 'w') as reformat:
            sequence = SeqIO.parse(fasta_alignments, "fasta")
            for each in sequence:
                print('>{}\n{}'.format(each.description, each.seq), file=reformat)

        fasta_df = pd.read_csv('reformated.fasta', header=None, sep='^')
        seq = fasta_df.iloc[1::2].reset_index(drop=True)
        header = fasta_df.iloc[0::2].reset_index(drop=True).replace(to_replace=r'>', value='', regex=True)
        seq = seq.rename({0:"seq"}, axis='columns')
        header = header.rename({0:"header"}, axis='columns')
        fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
        fasta_df = fasta_df.rename({0: 'header', 1: 'seq'}, axis='columns')
        seq = fasta_df['seq'].apply(lambda x: pd.Series(list(x)))
        fasta_df = pd.concat([header, seq], axis='columns', ignore_index=True)
        fasta_df = fasta_df.set_index(0)
        try:
            ref_series = fasta_df.loc['root']
            fasta_df = fasta_df.drop(['root']) #in all_vcf reference_seq needs to be removed
        except KeyError:
            print('Check that there is a "root" named')
            sys.exit(1)
        parsimony = fasta_df.loc[:, (fasta_df != fasta_df.iloc[0]).any()]
        parsimony_positions = list(parsimony)
        parse_df = fasta_df[parsimony_positions]
        ref_df = ref_series.to_frame()
        ref_df = ref_df.T
        out_df = pd.concat([parse_df, ref_df], join='inner')
        parsimonious_fasta = fasta_alignments.replace('.fasta', '_parsimonious.fasta')
        with open(parsimonious_fasta, 'w') as write_out:
            for index, row in out_df.iterrows():
                print('>{}'.format(row.name), file=write_out)
                for pos in row:
                    print(pos, end='', file=write_out)
                print("", file=write_out)
        self.parsimonious_fasta = parsimonious_fasta

        if not debug:
            os.remove('reformated.fasta')


if __name__ == '__main__':
    # Set multiprocessing start method for Python 3.12 compatibility
    multiprocessing.set_start_method('spawn', force=True)

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
        ---------------------------------------------------------
        Usage:
        vsnp3_fasta_to_snps_table.py -f *fasta -pn
        ---------------------------------------------------------

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True, help='Provide an alignment file in FASTA format')
    parser.add_argument('-p', '--parsimonious', action='store_true', dest='parsimonious', help='Only keep parsimonious SNPs from FASTA alignment file.  This is different than the uninformative SNPs removed via vSNP pipeline.  This is to be used when just working with an aligned FASTA file.')
    parser.add_argument('--show_groups', action='store_true', dest='show_groups', help='Show group names in SNP table')
    parser.add_argument('-n', '--hash_names', action='store_true', dest='hash_names', help='Hash FASTA names to rid of any RAxML illegal characters')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', help='Optional: Keep debugging files and run without pooling')
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.abspath(__file__), __version__))
    args = parser.parse_args()
    print('\n{} SET ARGUMENTS:'.format(os.path.basename(__file__)))
    print(args)

    initial_fasta = args.fasta

    if args.parsimonious:
        parsimonious = Parsimonious(args.fasta)
        args.fasta = parsimonious.parsimonious_fasta

    if args.hash_names:
        hash_names = Hash_Names(args.fasta)
        args.fasta = hash_names.hash_fasta()

    tree = Tree(fasta_alignments=args.fasta, debug=args.debug)

    if args.hash_names:
        tree.newick = hash_names.hash_tree_to_original_name(tree.newick)
    
    if args.parsimonious:
        fasta_alignments=parsimonious.parsimonious_fasta
    else:
        fasta_alignments=initial_fasta
    tables = Tables(fasta_alignments=fasta_alignments, tree=tree.newick, debug=False)
    tables.build_tables()