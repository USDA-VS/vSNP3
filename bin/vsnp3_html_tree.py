#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from collections import defaultdict
import os
import re
import math
import sys

from vsnp3_version import __version__

# Bootstrap resampling seed.  Fixed so support values are reproducible: RAxML is
# already seeded with -p 456123, so this was the last source of run-to-run
# variation in a published tree.  Overridable with --bootstrap_seed.
BOOTSTRAP_SEED = 456123

# SNP table columns are abs_pos strings, '<chrom>:<1-based position>'.
POSITION_COLUMN = re.compile(r'^.+:\d+$')

# Columns that sit alongside the positions and are not sequence data.
NON_POSITION_COLUMNS = ('Average Coverage', 'Grouping')


class NoPositionColumnsError(RuntimeError):
    '''The SNP table yielded no position columns, so no tree can be built.'''


def select_position_columns(df):
    '''
    First column (sample names) plus every position column.

    This used to keep only columns starting with 'NC_'.  Any reference whose
    contigs are not RefSeq NC_ accessions therefore yielded zero position
    columns and a star tree reporting "Total SNPs in tree: 0" as a normal
    result.  Nine of the shipped references are affected, among them
    Mycobacterium_orygis (CP063804), mtbc0_v1.1 (MTBC0), LT708304,
    para-CP033688 and four Brucella genomes on NZ_ or KN accessions.

    Raises rather than returning an empty frame: a tree with no SNPs in it is
    the plausible-looking wrong answer, which is worse than no tree.
    '''
    df_clean = df.copy()
    for column in NON_POSITION_COLUMNS:
        if column in df_clean.columns:
            df_clean = df_clean.drop(column, axis=1)

    name_column = df_clean.columns[0]
    positions = [col for col in df_clean.columns[1:]
                 if isinstance(col, str) and POSITION_COLUMN.match(col)]
    if not positions:
        offered = [str(c) for c in df_clean.columns[1:6]]
        raise NoPositionColumnsError(
            'no SNP position columns found in the table. Expected column names of '
            "the form '<chrom>:<position>', for example 'NC_002945.4:1295549'; "
            f'saw {offered}. No tree was written, because one built from zero '
            'positions would report every sample as identical.')
    return df_clean[[name_column] + positions]


class Node:
    def __init__(self, name=None):
        self.name = name
        self.children = []
        self.snps = []  # List of SNP descriptions
        self.x = 0
        self.y = 0
        self.parent = None
        self.excel_index = None
        self.leaf_indexes = []  # To track leaf node positions
        self.unrooted_x = 0
        self.unrooted_y = 0
        self.bootstrap = None  # Store bootstrap support value


# Rows of a SNP table that are not samples (besides the root in row 0): the
# per-position mapping quality, and the annotation row, which is named
# 'annotation' when a GenBank file was used and 'no annotations' when not.
NON_SAMPLE_ROWS = ("MQ", "annotation", "no annotations")


class SnpMatrix:
    '''
    The SNP content of one alignment table, as a bit matrix, with the
    agglomerative tree builder that both layouts and the bootstrap use.

    The tree is built exactly as before: repeatedly merge the pair of nodes
    that share the most SNP descriptors ('<pos>:<root>-><base>'), the first
    such pair in list order winning a tie, until nothing is shared.  What
    changed is the cost.  The original kept one Python set per node and
    intersected every pair on every round, O(S^3) set operations per tree,
    and built the tree 102 times per group (two layouts and 100 bootstrap
    replicates): 3.7 s per tree at 200 samples, so a 200-sample group
    needed most of the 15-minute timeout and a 400-sample group never
    finished.

    Here every descriptor is a column of a boolean matrix, the pairwise
    shared counts come from one matrix product, and a merge only has to
    update the counts of the new node and the rows whose best partner it
    consumed.  Descriptors with the same sample pattern are collapsed into
    one weighted column, so a clonal group of S samples has about 2S
    columns however many SNPs it carries.  A bootstrap replicate is the
    same matrix with per-column weights counting how many sampled
    descriptors fall in each pattern, so it needs no resampled DataFrame.

    Every quantity the original tree depends on is reproduced: the pair
    order, the strictly-greater tie-break, the sorted descriptor lists on
    each node (sorted as strings, duplicates removed on merge, kept on a
    leaf that is never merged), excel_index and leaf_indexes.
    '''

    def __init__(self, df_clean):
        sequences = df_clean.iloc[:, 1:]
        self.positions = sequences.columns.tolist()
        root_sequence = sequences.iloc[0].tolist()
        # Rows that are not samples: the MQ row, the annotation row (named
        # 'annotation' with a GenBank file, 'no annotations' without one) and,
        # after .iloc[1:], the root.  Without a gbk the all-NaN 'no annotations'
        # row used to pass this filter and become a leaf that differed from the
        # root at every position.
        samples = df_clean[~df_clean.iloc[:, 0].isin(NON_SAMPLE_ROWS)].iloc[1:]
        self.sample_names = samples.iloc[:, 0].tolist()
        self.sample_index = samples.index.tolist()
        self.values = samples.iloc[:, 1:].to_numpy(dtype=object)
        root = np.array(root_sequence, dtype=object)
        n_samples, n_positions = self.values.shape

        # A descriptor is present where the call differs from the root and is
        # not a no-call.  Object comparison, so every value compares exactly
        # as it did in get_snp_profile().
        if n_samples and n_positions:
            diff = (self.values != root[None, :]) & (self.values != '-')
        else:
            diff = np.zeros((n_samples, n_positions), dtype=bool)
        rows, cols = np.nonzero(diff)
        cell = self.values[rows, cols]
        codes, uniques = pd.factorize(cell, use_na_sentinel=False)
        n_uniques = max(len(uniques), 1)
        key = cols.astype(np.int64) * n_uniques + codes
        ukey, inverse = np.unique(key, return_inverse=True)
        ucol = ukey // n_uniques
        ucode = ukey % n_uniques
        strings = np.array(
            ["{pos}:{root_base}->{base}".format(pos=self.positions[c], root_base=root[c],
                                                base=uniques[k])
             for c, k in zip(ucol.tolist(), ucode.tolist())], dtype=object)
        # Sorted, as sorted() sorts str, so a rank order is a string order.
        if len(strings):
            self.desc, first, inv2 = np.unique(strings, return_index=True, return_inverse=True)
        else:
            self.desc, first, inv2 = strings, np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
        inv2 = np.asarray(inv2).reshape(-1)
        rank = inv2[inverse] if len(inverse) else np.zeros(0, dtype=np.int64)
        n_desc = len(self.desc)
        self.desc_column = ucol[first] if len(first) else np.zeros(0, dtype=np.int64)

        presence = np.zeros((n_samples, n_desc), dtype=bool)
        presence[rows, rank] = True
        # Per-sample descriptor lists with multiplicity, in string order: a
        # position sampled twice by the bootstrap is listed twice, as
        # get_snp_profile() listed it.
        splits = np.searchsorted(rows, np.arange(1, n_samples))
        self.leaf_ranks = [np.sort(r) for r in np.split(rank, splits)]

        # Collapse descriptors with the same sample pattern into one column.
        if n_desc:
            packed = np.packbits(presence, axis=0).T
            upat, pat = np.unique(packed, axis=0, return_inverse=True)
            self.pat_of_desc = np.asarray(pat).reshape(-1)
            self.patterns = np.unpackbits(upat, axis=1, count=n_samples).astype(bool).T
        else:
            self.pat_of_desc = np.zeros(0, dtype=np.int64)
            self.patterns = np.zeros((n_samples, 0), dtype=bool)
        self.n_patterns = self.patterns.shape[1]
        self.weight = np.bincount(self.pat_of_desc, minlength=self.n_patterns)

    def bootstrap_weight(self, sampled_positions):
        '''Pattern weights for a replicate: descriptors whose position was drawn.'''
        sampled = set(np.asarray(sampled_positions).tolist())
        column_drawn = np.fromiter((p in sampled for p in self.positions),
                                   dtype=bool, count=len(self.positions))
        drawn = column_drawn[self.desc_column]
        return np.bincount(self.pat_of_desc[drawn], minlength=self.n_patterns)

    def leaf_snps(self, sample_row):
        return self.desc[self.leaf_ranks[sample_row]].tolist()

    def _snps_of(self, pattern_mask):
        return self.desc[np.flatnonzero(pattern_mask[self.pat_of_desc])].tolist()

    def _shared_counts(self, membership, weight, full=False):
        '''membership: nodes x patterns bool; returns nodes x nodes shared counts.'''
        n = membership.shape[0]
        out = np.zeros((n, n), dtype=np.int32)
        if membership.shape[1] == 0:
            return out
        if full:
            # Every sample row: the same transposed matrix serves all 100
            # bootstrap replicates, so it is built once.
            right = getattr(self, '_right_full', None)
            if right is None:
                right = self._right_full = membership.T.astype(np.float32)
        else:
            right = membership.T.astype(np.float32)
        wf = weight.astype(np.float32)
        step = 2048
        for start in range(0, n, step):
            block = membership[start:start + step].astype(np.float32) * wf
            out[start:start + step] = np.rint(block @ right).astype(np.int32)
        return out

    @staticmethod
    def _initial_best(cnt):
        '''
        For each row, the largest count among the rows after it and the first
        row reaching it (-1 for the last row, which has none).  Blockwise, so
        the masked copy never holds the whole matrix.
        '''
        n = cnt.shape[0]
        rowmax = np.full(n, -1, dtype=np.int64)
        rowarg = np.full(n, -1, dtype=np.int64)
        cols = np.arange(n)
        step = 1024
        for start in range(0, n, step):
            rows = np.arange(start, min(start + step, n))
            block = np.where(cols[None, :] > rows[:, None], cnt[start:start + step], -1)
            rowmax[rows] = block.max(axis=1)
            rowarg[rows] = block.argmax(axis=1)
        rowarg[rowmax < 0] = -1
        return rowmax, rowarg

    def agglomerate(self, leaf_nodes, root, sample_rows, weight, materialize=True):
        '''
        Merge leaf_nodes (one per entry of sample_rows, in that order) under
        root.  weight gives the number of descriptors each pattern counts for.
        With materialize=False the nodes carry no descriptor strings and no
        leaf_indexes; the bootstrap only reads the topology.
        '''
        n = len(leaf_nodes)
        if n == 0:
            return root
        if n == 1:
            root.children.append(leaf_nodes[0])
            leaf_nodes[0].parent = root
            return root

        full = (n == self.patterns.shape[0]
                and bool(np.array_equal(sample_rows, np.arange(n))))
        vec = self.patterns if full else self.patterns[sample_rows]
        cnt = self._shared_counts(vec, weight, full=full)   # slot x slot shared counts
        vec = vec.copy()                                     # slot -> pattern membership
        wf = weight.astype(np.float32)
        nodes = list(leaf_nodes)                             # slot -> Node
        order = np.arange(n)                                 # alive slots, list order
        # Best partner of each row among the rows after it in list order:
        # the value and the first slot reaching it, as the nested loop found them.
        rowmax, rowarg = self._initial_best(cnt)
        position = np.empty(n, dtype=np.int64)
        attached_to_root = False

        while len(order) > 1:
            k = int(np.argmax(rowmax[order]))
            i = int(order[k])
            max_shared = int(rowmax[i])
            if max_shared <= 0:
                for slot in order.tolist():
                    root.children.append(nodes[slot])
                    nodes[slot].parent = root
                attached_to_root = True
                break
            j = int(rowarg[i])

            node1, node2 = nodes[i], nodes[j]
            shared = vec[i] & vec[j]
            private = vec[i] & ~shared
            internal = Node()
            if materialize:
                internal.leaf_indexes = sorted(node1.leaf_indexes + node2.leaf_indexes)
                internal.snps = self._snps_of(shared)
                node1.snps = self._snps_of(private)
                node2.snps = self._snps_of(vec[j] & ~shared)
            internal.children = [node1, node2]
            node1.parent = internal
            node2.parent = internal

            # The new node takes slot i and goes to the end of the list.
            others = order[(order != i) & (order != j)]
            z = i
            order = np.append(others, z)

            # Counts against the new node, |k & i & j| for every other node k.
            # Either summed over the patterns the new node carries, or taken
            # from the row node i already has and reduced by what i keeps for
            # itself, |k & i| - |k & (i - j)|.  The shared set is the whole
            # ancestral path of the pair and grows with depth; the private
            # part is one edge.  Whichever is smaller is used; both are exact.
            if len(others):
                idx_shared = np.flatnonzero(shared)
                idx_private = np.flatnonzero(private)
                if len(idx_private) <= len(idx_shared):
                    c = cnt[i, others].astype(np.int64)
                    if len(idx_private):
                        c = c - np.rint(vec[np.ix_(others, idx_private)].astype(np.float32)
                                        @ wf[idx_private]).astype(np.int64)
                elif len(idx_shared):
                    c = np.rint(vec[np.ix_(others, idx_shared)].astype(np.float32)
                                @ wf[idx_shared]).astype(np.int64)
                else:
                    c = np.zeros(len(others), dtype=np.int64)
                cnt[z, others] = c
                cnt[others, z] = c
                # The new node is last, so it becomes a row's best partner only
                # when strictly better than what the row already had.
                better = c > rowmax[others]
                rowmax[others[better]] = c[better]
                rowarg[others[better]] = z
                # Rows whose best partner was consumed are rescanned over their
                # successors, first-reaching-the-maximum winning as before.
                need = others[(rowarg[others] == i) | (rowarg[others] == j)]
                if len(need):
                    position[order] = np.arange(len(order))
                    sub = cnt[np.ix_(need, order)]
                    after = np.arange(len(order))[None, :] > position[need][:, None]
                    sub = np.where(after, sub, -1)
                    a = sub.argmax(axis=1)
                    rowmax[need] = sub[np.arange(len(need)), a]
                    rowarg[need] = order[a]
            vec[i] = shared
            nodes[i] = internal
            rowmax[z] = -1
            rowarg[z] = -1

        if not attached_to_root and len(order) == 1:
            root.children.append(nodes[int(order[0])])
            nodes[int(order[0])].parent = root
        return root

    def build_rectangular(self, weight=None, materialize=True):
        '''One leaf per sample row, in table order.'''
        leaf_nodes = []
        for row, (idx, name) in enumerate(zip(self.sample_index, self.sample_names)):
            node = Node(name)
            if materialize:
                node.snps = self.leaf_snps(row)
            node.excel_index = idx
            node.leaf_indexes = [idx]
            leaf_nodes.append(node)
        root = Node("Reference")
        return self.agglomerate(leaf_nodes, root, np.arange(len(leaf_nodes)),
                                self.weight if weight is None else weight, materialize)

    def build_unrooted(self):
        '''One leaf per distinct sequence; identical samples share a leaf.'''
        sequence_groups = defaultdict(list)
        for row, (idx, name) in enumerate(zip(self.sample_index, self.sample_names)):
            sequence_groups[tuple(self.values[row])].append((idx, name))
        leaf_nodes, sample_rows = [], []
        processed_seqs = set()
        for row, (idx, name) in enumerate(zip(self.sample_index, self.sample_names)):
            seq_key = tuple(self.values[row])
            if seq_key in processed_seqs:
                continue
            node = Node(name)
            node.snps = self.leaf_snps(row)
            node.excel_index = idx
            node.leaf_indexes = [idx]
            identical_group = sequence_groups[seq_key]
            if len(identical_group) > 1:
                node.all_names = sorted(name for _, name in identical_group)
            leaf_nodes.append(node)
            sample_rows.append(row)
            processed_seqs.add(seq_key)
        root = Node("Reference")
        return self.agglomerate(leaf_nodes, root, np.asarray(sample_rows, dtype=np.int64),
                                self.weight, True)


def parse_args():
    parser = argparse.ArgumentParser(description='Create interactive phylogenetic tree from SNP alignment')
    parser.add_argument('-f', '--file', required=True, help='Input Excel file containing SNP alignment')
    parser.add_argument('-o', '--output', default='phylo_tree.html', 
                       help='Output HTML file name (default: phylo_tree.html)')
    parser.add_argument('-v', '--version', action='version', version='{}: version {}'.format(os.path.basename(__file__), __version__))
    parser.add_argument('--bootstrap_seed', type=int, default=BOOTSTRAP_SEED,
                        help='Seed for bootstrap resampling.  Fixed by default so support values are reproducible.')
    parser.add_argument('-b', '--bootstrap', type=int, default=100,
                       help='Number of bootstrap replicates (default: 100)')
    return parser.parse_args()

def get_snp_profile(seq, root_seq, positions):
    """Get SNP profile for a sequence compared to root"""
    snps = []
    for i, (base, root_base) in enumerate(zip(seq, root_seq)):
        if base != root_base and base != '-':
            snps.append("{pos}:{root_base}->{base}".format(pos=positions[i], root_base=root_base, base=base))
    return sorted(snps)

def find_identical_samples(sequences, positions, sequences_df):
    """Group identical sequences before tree building"""
    identical_groups = {}
    processed = set()
    
    # Compare each sequence to others
    for idx1, seq1 in sequences.iterrows():
        if idx1 in processed:
            continue
            
        group = [idx1]
        sequence1 = tuple(seq1)
        
        # Find all identical sequences
        for idx2, seq2 in sequences.iterrows():
            if idx2 != idx1 and idx2 not in processed:
                if tuple(seq2) == sequence1:
                    group.append(idx2)
        
        if len(group) > 1:
            group_key = frozenset(group)
            identical_groups[group_key] = [sequences_df.iloc[i].iloc[0] for i in group]  # Store names
            processed.update(group)
        else:
            processed.add(idx1)
    
    return identical_groups

def assign_unrooted_coordinates_improved(node, angle=0, angle_span=2*math.pi, parent=None, level=0, processed_names=None):
    """
    Assign coordinates for unrooted layout with improved reference branch positioning
    and better handling of grouped samples
    """
    if processed_names is None:
        processed_names = set()
    
    # For root node, initialize position at center
    if parent is None:
        node.unrooted_x = 0
        node.unrooted_y = 0
        
        # Special handling for root node's children
        if node.children:
            total_leaves = sum(count_leaves(child) for child in node.children)
            current_angle = -math.pi/6  # Start from slightly right of top
            
            # Process children with adjusted angles
            for child in node.children:
                child_leaves = count_leaves(child)
                
                # Adjust angle span based on whether it's the reference node
                if child is node.children[0]:  # Reference node
                    child_span = math.pi/6  # Small fixed span for reference
                    fixed_angle = -math.pi/2  # Point straight up
                    branch_length = 60  # Longer fixed length for reference branch
                    child.unrooted_x = node.unrooted_x + branch_length * math.cos(fixed_angle)
                    child.unrooted_y = node.unrooted_y + branch_length * math.sin(fixed_angle)
                else:
                    # For other children, distribute them in the remaining space
                    # Updated to handle the case when total_leaves equals 1
                    denominator = max(total_leaves - 1, 1)  # Avoid division by zero
                    child_span = (2*math.pi - math.pi/3) * (child_leaves / denominator)
                    child_angle = current_angle + child_span/2
                    branch_length = len(child.snps) * 15 if child.snps else 30
                    child.unrooted_x = node.unrooted_x + branch_length * math.cos(child_angle)
                    child.unrooted_y = node.unrooted_y + branch_length * math.sin(child_angle)
                
                # Process subtrees
                if child.children:
                    next_angle = fixed_angle if child is node.children[0] else child_angle
                    next_span = child_span
                    assign_unrooted_coordinates_improved(
                        child,
                        angle=next_angle,
                        angle_span=next_span,
                        parent=node,
                        level=level+1,
                        processed_names=processed_names
                    )
                
                if child is not node.children[0]:
                    current_angle += child_span
    else:
        # Calculate branch length based on SNPs
        branch_length = len(node.snps) * 15 if node.snps else 30
        node.unrooted_x = parent.unrooted_x + branch_length * math.cos(angle)
        node.unrooted_y = parent.unrooted_y + branch_length * math.sin(angle)

        if node.children:
            # Identify and group identical terminal nodes
            terminal_groups = {}
            non_terminal_children = []
            
            for child in node.children:
                if not child.children:  # Terminal node
                    if hasattr(child, 'all_names'):
                        group_key = tuple(sorted(child.all_names))
                        if group_key not in terminal_groups:
                            terminal_groups[group_key] = []
                        terminal_groups[group_key].append(child)
                    else:
                        group_key = (child.name,)
                        if group_key not in terminal_groups:
                            terminal_groups[group_key] = []
                        terminal_groups[group_key].append(child)
                else:
                    non_terminal_children.append(child)
            
            # Process groups and non-terminal children
            nodes_to_process = non_terminal_children + [group[0] for group in terminal_groups.values()]
            
            # Calculate total weight
            total_weight = sum(count_leaves(child) for child in nodes_to_process)
            
            if total_weight > 0:
                # Distribute angles
                current_angle = angle - angle_span/2
                for child in nodes_to_process:
                    weight = count_leaves(child)
                    child_span = (weight / total_weight) * angle_span
                    child_angle = current_angle + child_span/2
                    
                    # Position the child
                    if hasattr(child, 'all_names'):
                        # For grouped nodes, position them closer together
                        group_key = tuple(sorted(child.all_names))
                        group = terminal_groups[group_key]
                        base_x = node.unrooted_x + branch_length * math.cos(child_angle)
                        base_y = node.unrooted_y + branch_length * math.sin(child_angle)
                        
                        # Position all nodes in group at same base position
                        for i, grouped_node in enumerate(group):
                            grouped_node.unrooted_x = base_x
                            grouped_node.unrooted_y = base_y
                    else:
                        # Process non-grouped node normally
                        assign_unrooted_coordinates_improved(
                            child,
                            angle=child_angle,
                            angle_span=child_span,
                            parent=node,
                            level=level+1,
                            processed_names=processed_names
                        )
                    
                    current_angle += child_span

def calculate_bootstrap_support(df, n_replicates=100, bootstrap_seed=BOOTSTRAP_SEED, matrix=None):
    """
    Bootstrap support per clade, keyed as leaf_bits() keys them.

    Positions are drawn exactly as before (same generator, same call), but a
    replicate is built from the SnpMatrix of the full table with the drawn
    positions as weights, rather than from a resampled copy of the DataFrame.
    The set of descriptors a sample carries in the replicate is the same
    either way, so each replicate tree is the same tree.
    """
    sequences = df.iloc[:, 1:]
    positions = sequences.columns.tolist()
    if matrix is None:
        matrix = SnpMatrix(select_position_columns(df))
    bits_of = leaf_bits(matrix.sample_names)

    # Store all branch configurations
    branch_configs = defaultdict(int)
    rng = np.random.default_rng(bootstrap_seed)

    for _ in range(n_replicates):
        # Randomly sample positions with replacement
        bootstrap_positions = rng.choice(positions, size=len(positions), replace=True)
        bootstrap_root = matrix.build_rectangular(
            weight=matrix.bootstrap_weight(bootstrap_positions), materialize=False)
        get_branch_configs(bootstrap_root, branch_configs, bits_of)

    # Convert counts to percentages
    for config in branch_configs:
        branch_configs[config] = (branch_configs[config] / n_replicates) * 100

    return branch_configs

def leaf_bits(names):
    '''
    One bit per leaf name.  A clade is keyed by the integer whose set bits
    are its leaves.  This used to be '|'.join(sorted(leaf names)), rebuilt
    from a full leaf walk at every internal node of every tree -- quadratic
    in the sample count for a chain-like tree, and 102 trees per group.
    The key still identifies exactly the set of leaf names, so the support
    values are unchanged.
    '''
    return {name: 1 << i for i, name in enumerate(sorted(set(names), key=str))}


def _clade_bits(node, configs, bits_of):
    if not node.children:
        return bits_of[node.name]
    bits = 0
    for child in node.children:
        bits |= _clade_bits(child, configs, bits_of)
    if bits & (bits - 1) and configs is not None:      # more than one leaf
        configs[bits] += 1
    return bits


def get_branch_configs(node, configs, bits_of=None):
    """Count each clade (set of leaf names) of the tree in configs"""
    if bits_of is None:
        bits_of = leaf_bits(leaf.name for leaf in get_leaves(node))
    _clade_bits(node, configs, bits_of)

def get_leaves(node):
    """Get all leaf nodes in the subtree"""
    if not node.children:
        return [node]
    leaves = []
    for child in node.children:
        leaves.extend(get_leaves(child))
    return leaves

def add_bootstrap_to_node(node, bootstrap_values, bits_of=None):
    """Add bootstrap support values to tree nodes"""
    if bits_of is None:
        bits_of = leaf_bits(leaf.name for leaf in get_leaves(node))

    def visit(current):
        if not current.children:
            return bits_of[current.name]
        bits = 0
        for child in current.children:
            bits |= visit(child)
        if bits & (bits - 1):
            current.bootstrap = bootstrap_values.get(bits, 0)
        return bits
    visit(node)

def assign_coordinates(node, x=0, y=0, dy=50, y_positions=None):
    """Assign x,y coordinates to tree nodes"""
    if y_positions is None:
        y_positions = []
    
    node.x = x
    
    if not node.children:
        # Leaf node - position based on Excel order
        if y_positions:
            node.y = y_positions[-1] + dy
        else:
            node.y = y
        y_positions.append(node.y)
    else:
        # Sort children by maximum leaf index to reverse the order
        node.children.sort(key=lambda n: max(n.leaf_indexes) if n.leaf_indexes else 0, reverse=True)
        
        child_ys = []
        for child in node.children:
            # Branch length proportional to number of SNPs
            dx = len(child.snps) * 50  # Scale factor
            assign_coordinates(child, x + dx, y, dy, y_positions)
            child_ys.append(child.y)
        node.y = sum(child_ys) / len(child_ys) if child_ys else y
    
    return node

def get_snp_path(node):
    """Get the complete SNP path from root to node"""
    path = []
    current = node
    while current.parent is not None:
        path.extend(current.snps)
        current = current.parent
    return sorted(path)

def get_full_snp_profile(node):
    """Get the complete SNP profile including all ancestor SNPs"""
    profile = set()
    current = node
    while current:
        profile.update(current.snps)
        current = current.parent
    return sorted(list(profile))

def get_terminal_snps(node):
    """Get the complete set of SNPs from root to this node"""
    snps = set()
    current = node
    while current:
        snps.update(current.snps)
        current = current.parent
    return frozenset(snps)

def count_leaves(node):
    """
    Count the leaf nodes in the subtree.

    Iterative.  A chain-like group is as deep as it has samples, and the
    recursive form re-entered the interpreter through sum() over a generator
    at every level, so sys.setrecursionlimit() did not protect it: a chain of
    3,000 samples still raised RecursionError at the C-level limit.
    """
    total = 0
    stack = [node]
    while stack:
        current = stack.pop()
        if not current.children:
            total += 1
        else:
            stack.extend(current.children)
    return total

def create_hover_text(snps, bootstrap=None):
    """Create hover text for branches including bootstrap values"""
    text = "SNPs ({snp_count}):<br>".format(snp_count=len(snps)) + "<br>".join(snps)
    if bootstrap is not None:
        text = "Bootstrap: {bootstrap:.1f}%<br><br>".format(bootstrap=bootstrap) + text
    return text

def create_hover_points(snps, bootstrap=None, chunk_size=30):
    """Split SNPs into chunks for multiple hover points"""
    chunks = []
    for i in range(0, len(snps), chunk_size):
        chunk = snps[i:i + chunk_size]
        hover_text = "SNPs {start}-{end} of {total}:<br>".format(
            start=i+1,
            end=min(i+chunk_size, len(snps)),
            total=len(snps)
        ) + "<br>".join(chunk)
        if bootstrap is not None:
            hover_text = "Bootstrap: {bootstrap:.1f}%<br><br>".format(bootstrap=bootstrap) + hover_text
        chunks.append(hover_text)
    return chunks

def add_bootstrap_annotation(fig, x, y, bootstrap, offset=10):
    """Add bootstrap value annotation with CSS class for toggling"""
    if bootstrap is not None:
        fig.add_annotation(
            x=x,
            y=y + offset,
            text="{bootstrap:.0f}".format(bootstrap=bootstrap),
            showarrow=False,
            font=dict(size=10),
            yanchor='bottom',
            # Add name for identification in JavaScript
            name='bootstrap_annotation'
        )

def get_sample_profile(node):
    """Get the complete SNP profile for a sample, including ancestral SNPs"""
    profile = []
    current = node
    while current:
        profile.extend(current.snps)
        current = current.parent
    return frozenset(profile)

def group_identical_samples(nodes):
    """
    Group samples with identical SNP profiles into polytomies.
    Returns a dictionary mapping SNP profiles to lists of nodes.
    """
    groups = defaultdict(list)
    for node in nodes:
        if not node.children:  # Only group leaf nodes
            profile = get_sample_profile(node)
            groups[profile].append(node)
    return groups

def group_identical_terminals(node):
    """
    Group terminal nodes with identical SNP profiles into polytomies.
    Returns a dictionary mapping SNP profiles to lists of nodes.
    """
    groups = defaultdict(list)
    terminal_nodes = collect_terminal_nodes(node)
    for terminal in terminal_nodes:
        profile = get_terminal_profile(terminal)
        groups[profile].append(terminal)
    return groups

def get_terminal_profile(node):
    """Get the complete SNP profile for a terminal node, including ancestral SNPs"""
    profile = []
    current = node
    while current:
        profile.extend(current.snps)
        current = current.parent
    return frozenset(profile)

def collect_terminal_nodes(node):
    """Collect all terminal nodes from a subtree"""
    if not node.children:
        return [node]
    terminals = []
    for child in node.children:
        terminals.extend(collect_terminal_nodes(child))
    return terminals

def get_complete_profile(node):
    """Get the complete SNP profile as a frozenset for comparison"""
    profile = []
    current = node
    while current:
        profile.extend(current.snps)
        current = current.parent
    return frozenset(sorted(profile))

def get_complete_snp_profile(node):
    """Get complete SNP profile including all ancestral SNPs"""
    profile = []
    current = node
    while current and current.parent:  # Stop before root
        profile.extend(current.snps)
        current = current.parent
    return frozenset(sorted(profile))

def process_identical_samples(root):
    """Pre-process the tree to identify and mark identical samples"""
    # Collect all leaf nodes
    leaves = get_leaves(root)
    
    # Group leaves by their complete SNP profile
    groups = defaultdict(list)
    for leaf in leaves:
        profile = get_complete_snp_profile(leaf)
        groups[profile].append(leaf)
    
    # Store groups for later use
    root.identical_groups = groups
    return groups

def assign_unrooted_coordinates(node, start_angle=0, angle_range=2*math.pi, level=0, parent_x=0, parent_y=0):
    """Assign coordinates for unrooted layout"""
    snp_scale = 15
    total_leaves = count_leaves(node)
    
    if level == 0:
        # Root node at center
        node.unrooted_x = 0
        node.unrooted_y = 0
    else:
        # Calculate position
        angle = start_angle + (angle_range / 2)
        branch_length = len(node.snps) * snp_scale if node.snps else 30
        
        # Set position using polar coordinates
        node.unrooted_x = parent_x + branch_length * math.cos(angle)
        node.unrooted_y = parent_y + branch_length * math.sin(angle)
    
    current_angle = start_angle
    
    if node.children:
        # Group identical terminal nodes
        terminal_groups = {}
        for child in node.children:
            if not child.children:  # Only for terminal nodes
                snp_set = get_terminal_snps(child)
                if snp_set not in terminal_groups:
                    terminal_groups[snp_set] = []
                terminal_groups[snp_set].append(child)
        
        # Process children
        processed_terminals = set()
        
        for child in node.children:
            child_leaves = count_leaves(child)
            # Avoid division by zero
            if total_leaves > 0:
                child_angle_range = (angle_range * child_leaves) / total_leaves
            else:
                child_angle_range = angle_range
            
            if child.children:
                # Internal node - process normally
                next_angle = assign_unrooted_coordinates(
                    child,
                    current_angle,
                    child_angle_range,
                    level + 1,
                    node.unrooted_x,
                    node.unrooted_y
                )
                current_angle = next_angle
            else:
                # Terminal node - check if part of a group
                snp_set = get_terminal_snps(child)
                if snp_set not in processed_terminals:
                    processed_terminals.add(snp_set)
                    group = terminal_groups[snp_set]
                    
                    # Assign same coordinates to all nodes in group
                    next_angle = assign_unrooted_coordinates(
                        child,
                        current_angle,
                        child_angle_range,
                        level + 1,
                        node.unrooted_x,
                        node.unrooted_y
                    )
                    
                    # Copy coordinates to other nodes in group
                    for other in group[1:]:
                        other.unrooted_x = child.unrooted_x
                        other.unrooted_y = child.unrooted_y
                    
                    current_angle = next_angle
    
    return current_angle

def build_tree(df):
    """Alias for backward compatibility, redirects to build_tree_rectangular"""
    return build_tree_rectangular(df)

def build_tree_rectangular(df):
    """Build phylogenetic tree optimized for rectangular layout"""
    return SnpMatrix(select_position_columns(df)).build_rectangular()

def build_tree_unrooted(df):
    """Build phylogenetic tree optimized for unrooted layout"""
    return SnpMatrix(select_position_columns(df)).build_unrooted()

def create_tree_visualization(df, output_file, input_filename, sample_coverage_dict=None, bootstrap_seed=BOOTSTRAP_SEED):
    """Create interactive tree visualization using plotly"""
    df_for_tree = df.copy()
    coverage_dict = {}

    if sample_coverage_dict:
        coverage_dict = sample_coverage_dict
    elif 'Average Coverage' in df_for_tree.columns:
        try:
            mask = ~df_for_tree.iloc[:, 0].isin(("Reference",) + NON_SAMPLE_ROWS)
            sample_names = df_for_tree.loc[mask, df_for_tree.columns[0]]
            coverage_values = df_for_tree.loc[mask, 'Average Coverage']
            coverage_dict = dict(zip(sample_names, coverage_values))
            df_for_tree = df_for_tree.drop('Average Coverage', axis=1)
        except Exception as e:
            print("Warning: Could not extract coverage data: {error}".format(error=e))
            coverage_dict = {}
    
    # One matrix serves both layouts and every bootstrap replicate.
    matrix = SnpMatrix(select_position_columns(df_for_tree))
    # Every walk of the tree here is recursive, and a chain-like group is as
    # deep as it has samples, so the default limit of 1,000 frames failed
    # ("maximum recursion depth exceeded") on such groups past a few hundred
    # samples.  Python 3.11+ keeps pure-Python frames off the C stack.
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 20 * len(matrix.sample_names) + 1000))
    rect_root = matrix.build_rectangular()
    unrooted_root = matrix.build_unrooted()
    bits_of = leaf_bits(matrix.sample_names)

    sequences = df_for_tree.iloc[:, 1:]
    positions = sequences.columns.tolist()

    bootstrap_values = defaultdict(int)
    n_replicates = 100
    # Seeded: an unseeded np.random.choice made the published support values differ
    # on every run over identical data.  Two runs of the same 6-sample dataset gave
    # the same clade 100% and then 98%.
    rng = np.random.default_rng(bootstrap_seed)

    for _ in range(n_replicates):
        # The same draw as before.  The replicate tree is built from the
        # drawn positions as pattern weights rather than from a resampled
        # DataFrame: a sample carries the same descriptor set either way.
        bootstrap_positions = rng.choice(positions, size=len(positions), replace=True)
        bootstrap_root = matrix.build_rectangular(
            weight=matrix.bootstrap_weight(bootstrap_positions), materialize=False)
        get_branch_configs(bootstrap_root, bootstrap_values, bits_of)

    for config in bootstrap_values:
        bootstrap_values[config] = (bootstrap_values[config] / n_replicates) * 100

    add_bootstrap_to_node(rect_root, bootstrap_values, bits_of)
    add_bootstrap_to_node(unrooted_root, bootstrap_values, bits_of)
    
    fig_rect = go.Figure()
    fig_unrooted = go.Figure()
    
    assign_coordinates(rect_root)
    assign_unrooted_coordinates_improved(unrooted_root)
    
    add_rectangular_traces(fig_rect, rect_root, sample_coverage_dict=coverage_dict)
    add_unrooted_traces_improved(fig_unrooted, unrooted_root, sample_coverage_dict=coverage_dict)
    
    leaf_nodes = get_leaves(unrooted_root)
    x_coords = [node.unrooted_x for node in leaf_nodes]
    y_coords = [node.unrooted_y for node in leaf_nodes]
    
    max_range = 1.2 * max(
        max(abs(max(x_coords) if x_coords else 0), abs(min(x_coords) if x_coords else 0)),
        max(abs(max(y_coords) if y_coords else 0), abs(min(y_coords) if y_coords else 0))
    )
    
    # Calculate total SNPs in the tree
    def count_total_snps(node):
        total = len(node.snps)
        for child in node.children:
            total += count_total_snps(child)
        return total
    
    total_snps = count_total_snps(rect_root)
    
    # Add SNP scale bar for rectangular tree
    snp_scale = 50  # Must match the scale factor used in assign_coordinates
    scale_length = 5  # Number of SNPs represented by the scale bar
    
    # Calculate the minimum y position of all leaves and place scale below
    leaf_nodes_rect = get_leaves(rect_root)
    min_y = min(node.y for node in leaf_nodes_rect) if leaf_nodes_rect else 0
    scale_y = min_y - 100  # Place scale below the tree
    
    # Position scale in bottom left corner
    scale_x = min(node.x for node in leaf_nodes_rect) if leaf_nodes_rect else 0
    
    # Add scale bar with refined styling
    fig_rect.add_shape(
        type="line",
        x0=scale_x,
        y0=scale_y,
        x1=scale_x + scale_length * snp_scale,
        y1=scale_y,
        line=dict(color="black", width=1.5)
    )
    
    # Add refined scale labels and the total SNPs annotation.  Appended in one
    # assignment: each add_annotation() re-validates every annotation already
    # on the figure, 5 s for these six at 1,000 samples.
    scale_annotations = _FigureBatch()
    scale_annotations.add_annotation(
        x=scale_x - 5,
        y=scale_y,
        text="0",
        showarrow=False,
        font=dict(size=9),
        xanchor='right',
        yanchor='middle'
    )
    scale_annotations.add_annotation(
        x=scale_x + scale_length * snp_scale + 5,
        y=scale_y,
        text="{scale_length} SNPs".format(scale_length=scale_length),
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    scale_annotations.add_annotation(
        x=scale_x,
        y=scale_y - 25,
        text="Total SNPs in tree: {total_snps}".format(total_snps=total_snps),
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    scale_annotations.flush(fig_rect)
    
    # Add SNP scale bar for unrooted tree
    unrooted_snp_scale = 15  # Must match scale in assign_unrooted_coordinates
    unrooted_scale_x = -max_range + 30
    unrooted_scale_y = -max_range + 30
    
    # Add scale bar with refined styling
    fig_unrooted.add_shape(
        type="line",
        x0=unrooted_scale_x,
        y0=unrooted_scale_y,
        x1=unrooted_scale_x + scale_length * unrooted_snp_scale,
        y1=unrooted_scale_y,
        line=dict(color="black", width=1.5)
    )
    
    # Add refined scale labels and the total SNPs annotation for unrooted tree
    scale_annotations = _FigureBatch()
    scale_annotations.add_annotation(
        x=unrooted_scale_x - 5,
        y=unrooted_scale_y,
        text="0",
        showarrow=False,
        font=dict(size=9),
        xanchor='right',
        yanchor='middle'
    )
    scale_annotations.add_annotation(
        x=unrooted_scale_x + scale_length * unrooted_snp_scale + 5,
        y=unrooted_scale_y,
        text="{scale_length} SNPs".format(scale_length=scale_length),
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    scale_annotations.add_annotation(
        x=unrooted_scale_x,
        y=unrooted_scale_y - 25,
        text="Total SNPs in tree: {total_snps}".format(total_snps=total_snps),
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    scale_annotations.flush(fig_unrooted)
    
    # Get base filename for title
    base_filename = os.path.splitext(os.path.basename(input_filename))[0]
    
    # Configure layouts
    layout_common = dict(
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=100, r=100, t=100, b=100),
        hovermode='closest'
    )
    
    title_text = ("{base_filename}<br>"
                 "Hover over branches/dots to see SNPs and bootstrap values<br>"
                 "Hover over sample names to see coverage information").format(base_filename=base_filename)
    
    fig_rect.update_layout(
        title=dict(text=title_text, x=0.5),
        width=1200,
        height=800,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        **layout_common
    )
    
    fig_unrooted.update_layout(
        title=dict(text=title_text + " - Unrooted Layout", x=0.5),
        width=1200,
        height=1200,
        xaxis=dict(
            range=[-max_range, max_range],
            showgrid=False,
            zeroline=False,
            showticklabels=False,
            scaleanchor="y",
            scaleratio=1
        ),
        yaxis=dict(
            range=[-max_range, max_range],
            showgrid=False,
            zeroline=False,
            showticklabels=False
        ),
        **layout_common
    )
    
    # Create HTML with toggle buttons
    html_content = """
    <html>
    <head>
        <title>Phylogenetic Tree Visualization</title>
        <style>
            .button {{
                padding: 10px 20px;
                margin: 10px;
                cursor: pointer;
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 4px;
            }}
            .button:hover {{
                background-color: #45a049;
            }}
            .button.inactive {{
                background-color: #f44336;
            }}
            .button.inactive:hover {{
                background-color: #d32f2f;
            }}
            .control-panel {{
                text-align: center;
                margin: 20px 0;
                padding: 15px;
                background-color: #f5f5f5;
                border-radius: 8px;
                border: 1px solid #ddd;
            }}
            .control-group {{
                display: inline-block;
                margin: 0 20px;
            }}
            .control-label {{
                font-weight: bold;
                margin-right: 10px;
                display: inline-block;
            }}
        </style>
    </head>
    <body>
        <div class="control-panel">
            <div class="control-group">
                <span class="control-label">Tree Layout:</span>
                <button class="button" onclick="toggleView()" id="toggleButton">Unrooted</button>
            </div>
            <div class="control-group">
                <span class="control-label">Bootstrap Values:</span>
                <button class="button" onclick="toggleBootstrap()" id="bootstrapButton">Hide Bootstrap</button>
            </div>
        </div>
        <div id="rectangular-tree">{rectangular_tree_html}</div>
        <div id="unrooted-tree" style="display: none;">{unrooted_tree_html}</div>

        <script>
            let bootstrapVisible = true;

            function toggleView() {{
                var rect = document.getElementById('rectangular-tree');
                var unrooted = document.getElementById('unrooted-tree');
                var button = document.getElementById('toggleButton');

                if (rect.style.display === 'none') {{
                    rect.style.display = 'block';
                    unrooted.style.display = 'none';
                    button.textContent = 'Unrooted';
                }} else {{
                    rect.style.display = 'none';
                    unrooted.style.display = 'block';
                    button.textContent = 'Rectangular';
                }}
            }}

            function toggleBootstrap() {{
                console.log('Bootstrap toggle clicked, current state:', bootstrapVisible);
                var button = document.getElementById('bootstrapButton');
                bootstrapVisible = !bootstrapVisible;

                if (bootstrapVisible) {{
                    button.textContent = 'Hide Bootstrap';
                    button.classList.remove('inactive');
                }} else {{
                    button.textContent = 'Show Bootstrap';
                    button.classList.add('inactive');
                }}

                // Toggle bootstrap annotations in both tree views
                toggleBootstrapInPlot('rectangular-tree');
                toggleBootstrapInPlot('unrooted-tree');
            }}

            function findPlotlyDiv(containerId) {{
                // Find the actual Plotly div within the container
                var container = document.getElementById(containerId);
                if (!container) return null;

                // Look for div with class 'plotly-graph-div' or similar
                var plotlyDiv = container.querySelector('div[class*="plotly"]') ||
                               container.querySelector('.js-plotly-plot') ||
                               container.querySelector('div').children[0];

                return plotlyDiv;
            }}

            function toggleBootstrapInPlot(containerId) {{
                console.log('Attempting to toggle bootstrap in:', containerId);

                var plotDiv = findPlotlyDiv(containerId);
                if (!plotDiv) {{
                    console.log('Plot div not found for:', containerId);
                    return;
                }}

                console.log('Found plot div:', plotDiv);
                console.log('Plot data available:', !!plotDiv.data);
                console.log('Plot layout available:', !!plotDiv.layout);

                if (plotDiv.layout && plotDiv.layout.annotations) {{
                    console.log('Found', plotDiv.layout.annotations.length, 'annotations');

                    var updatedAnnotations = plotDiv.layout.annotations.map(function(annotation, index) {{
                        // Check if this is a bootstrap annotation (contains only digits)
                        var isBootstrap = annotation.text && /^\\d+$/.test(annotation.text.toString());

                        if (isBootstrap) {{
                            console.log('Found bootstrap annotation:', annotation.text, 'at index:', index);
                            return Object.assign({{}}, annotation, {{
                                visible: bootstrapVisible
                            }});
                        }}
                        return annotation;
                    }});

                    console.log('Updating annotations, bootstrap visible:', bootstrapVisible);

                    Plotly.relayout(plotDiv, {{'annotations': updatedAnnotations}})
                        .then(function() {{
                            console.log('Successfully updated annotations');
                        }})
                        .catch(function(error) {{
                            console.error('Error updating annotations:', error);
                        }});
                }} else {{
                    console.log('No annotations found in layout');
                }}
            }}

            // Initialize after page loads
            window.addEventListener('load', function() {{
                console.log('Page loaded, initializing...');

                // Wait for plots to fully render
                setTimeout(function() {{
                    console.log('Checking for plots...');
                    var rectPlot = findPlotlyDiv('rectangular-tree');
                    var unrootedPlot = findPlotlyDiv('unrooted-tree');

                    console.log('Rectangular plot found:', !!rectPlot);
                    console.log('Unrooted plot found:', !!unrootedPlot);

                    if (rectPlot && rectPlot.layout) {{
                        console.log('Rectangular plot annotations:', rectPlot.layout.annotations ? rectPlot.layout.annotations.length : 0);
                    }}
                    if (unrootedPlot && unrootedPlot.layout) {{
                        console.log('Unrooted plot annotations:', unrootedPlot.layout.annotations ? unrootedPlot.layout.annotations.length : 0);
                    }}
                }}, 2000);
            }});
        </script>
    </body>
    </html>
    """.format(
        # plotly.js is embedded once, with the first figure in the document;
        # the second figure reuses it.  Each file used to carry the bundle
        # twice (about 9 MB per group) and also loaded "plotly-latest" from the
        # CDN, a URL frozen at plotly.js 1.58 and deprecated, which needed
        # network access and was then overwritten by the embedded bundle.
        rectangular_tree_html=fig_rect.to_html(full_html=False, include_plotlyjs=True),
        unrooted_tree_html=fig_unrooted.to_html(full_html=False, include_plotlyjs=False)
    )
    
    # Save as interactive HTML
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

class _FigureBatch:
    '''
    Collects traces and annotations and adds them to the figure in one call
    each.  fig.add_annotation() re-validates and deep-copies every
    annotation already on the figure, so adding S labels one at a time is
    O(S^2): 13 s of a 100-sample tree, and sixteen times that at 400.  The
    figure ends up with the same traces and annotations in the same order.
    '''

    def __init__(self):
        self.traces = []
        self.annotations = []

    def add_trace(self, trace):
        self.traces.append(trace)

    def add_annotation(self, **kwargs):
        self.annotations.append(kwargs)

    def flush(self, fig):
        if self.traces:
            fig.add_traces(self.traces)
        if self.annotations:
            fig.layout.annotations = tuple(fig.layout.annotations) + tuple(self.annotations)


def add_rectangular_traces(fig, node, parent=None, sample_coverage_dict=None):
    """Add traces for rectangular layout with bootstrap values and coverage information"""
    if not isinstance(fig, _FigureBatch):
        batch = _FigureBatch()
        add_rectangular_traces(batch, node, parent, sample_coverage_dict)
        batch.flush(fig)
        return
    if parent:
        if node.snps:
            sorted_snps = sorted(node.snps)
            x_start, x_end = parent.x, node.x
            y_pos = node.y
            
            if len(sorted_snps) > 40:
                hover_chunks = create_hover_points(sorted_snps, node.bootstrap)
                
                fig.add_trace(go.Scatter(
                    x=[x_start, x_end],
                    y=[y_pos, y_pos],
                    mode='lines',
                    line=dict(color='rgb(40,40,40)', width=2),
                    hoverinfo='skip',
                    showlegend=False
                ))
                
                if node.bootstrap is not None:
                    x_mid = (x_start + x_end) / 2
                    add_bootstrap_annotation(fig, x_mid, y_pos, node.bootstrap)
                
                for i, hover_text in enumerate(hover_chunks):
                    x_pos = x_start + (x_end - x_start) * ((i + 1) / (len(hover_chunks) + 1))
                    fig.add_trace(go.Scatter(
                        x=[x_pos],
                        y=[y_pos],
                        mode='markers',
                        marker=dict(size=10, color='rgb(40,40,40)', opacity=0.5),
                        hoverinfo='text',
                        hovertext=hover_text,
                        hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace", align='left'),
                        showlegend=False
                    ))
            else:
                hover_text = create_hover_text(sorted_snps, node.bootstrap)
                fig.add_trace(go.Scatter(
                    x=[x_start, x_end],
                    y=[y_pos, y_pos],
                    mode='lines',
                    line=dict(color='rgb(40,40,40)', width=2),
                    hoverinfo='text',
                    hovertext=hover_text,
                    hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace", align='left'),
                    showlegend=False
                ))
                
                if node.bootstrap is not None:
                    x_mid = (x_start + x_end) / 2
                    add_bootstrap_annotation(fig, x_mid, y_pos, node.bootstrap)
        else:
            fig.add_trace(go.Scatter(
                x=[parent.x, node.x],
                y=[node.y, node.y],
                mode='lines',
                line=dict(color='rgb(40,40,40)', width=2),
                hoverinfo='skip',
                showlegend=False
            ))
        
        if parent.y != node.y:
            fig.add_trace(go.Scatter(
                x=[parent.x, parent.x],
                y=[parent.y, node.y],
                mode='lines',
                line=dict(color='rgb(40,40,40)', width=2),
                hoverinfo='skip',
                showlegend=False
            ))
    
    if node.name:
        # Add coverage information to hover text if available
        hover_text = node.name
        if sample_coverage_dict and node.name in sample_coverage_dict:
            hover_text = "{name}<br>Coverage: {coverage}X".format(
                name=node.name,
                coverage=sample_coverage_dict[node.name]
            )
        
        fig.add_annotation(
            x=node.x,
            y=node.y,
            text=node.name,
            showarrow=False,
            xanchor='left',
            xshift=10,
            font=dict(size=12),
            hovertext=hover_text,
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace"),
        )
    
    for child in node.children:
        add_rectangular_traces(fig, child, node, sample_coverage_dict)

def add_unrooted_traces_improved(fig, node, parent=None, label_offset_base=20, processed_positions=None, sample_coverage_dict=None):
    """Add traces for improved unrooted layout with coverage information"""
    if not isinstance(fig, _FigureBatch):
        batch = _FigureBatch()
        add_unrooted_traces_improved(batch, node, parent, label_offset_base,
                                     processed_positions, sample_coverage_dict)
        batch.flush(fig)
        return
    if processed_positions is None:
        processed_positions = set()

    if parent:
        if node.snps:
            sorted_snps = sorted(node.snps)
            x_start, x_end = parent.unrooted_x, node.unrooted_x
            y_start, y_end = parent.unrooted_y, node.unrooted_y
            
            if len(sorted_snps) > 40:
                hover_chunks = create_hover_points(sorted_snps, node.bootstrap)
                
                fig.add_trace(go.Scatter(
                    x=[x_start, x_end],
                    y=[y_start, y_end],
                    mode='lines',
                    line=dict(color='rgb(40,40,40)', width=2),
                    hoverinfo='skip',
                    showlegend=False
                ))
                
                for i, hover_text in enumerate(hover_chunks):
                    t = (i + 1) / (len(hover_chunks) + 1)
                    x_pos = x_start + (x_end - x_start) * t
                    y_pos = y_start + (y_end - y_start) * t
                    
                    fig.add_trace(go.Scatter(
                        x=[x_pos],
                        y=[y_pos],
                        mode='markers',
                        marker=dict(size=8, color='rgb(40,40,40)', opacity=0.5),
                        hoverinfo='text',
                        hovertext=hover_text,
                        hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace"),
                        showlegend=False
                    ))
            else:
                hover_text = create_hover_text(sorted_snps, node.bootstrap)
                fig.add_trace(go.Scatter(
                    x=[x_start, x_end],
                    y=[y_start, y_end],
                    mode='lines',
                    line=dict(color='rgb(40,40,40)', width=2),
                    hoverinfo='text',
                    hovertext=hover_text,
                    hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace"),
                    showlegend=False
                ))
        else:
            fig.add_trace(go.Scatter(
                x=[parent.unrooted_x, node.unrooted_x],
                y=[parent.unrooted_y, node.unrooted_y],
                mode='lines',
                line=dict(color='rgb(40,40,40)', width=2),
                hoverinfo='skip',
                showlegend=False
            ))

        if node.bootstrap is not None:
            x_mid = (parent.unrooted_x + node.unrooted_x) / 2
            y_mid = (parent.unrooted_y + node.unrooted_y) / 2
            dx = node.unrooted_x - parent.unrooted_x
            dy = node.unrooted_y - parent.unrooted_y
            length = math.sqrt(dx*dx + dy*dy)
            if length > 0:
                offset = 10
                x_offset = -dy/length * offset
                y_offset = dx/length * offset
                fig.add_annotation(
                    x=x_mid + x_offset,
                    y=y_mid + y_offset,
                    text="{bootstrap:.0f}".format(bootstrap=node.bootstrap),
                    showarrow=False,
                    font=dict(size=10),
                    bgcolor='rgba(255,255,255,0.8)',
                    # Add name for identification in JavaScript
                    name='bootstrap_annotation'
                )

    if node.name:
        pos_key = (round(node.unrooted_x, 2), round(node.unrooted_y, 2))
        if pos_key not in processed_positions:
            processed_positions.add(pos_key)
            
            if hasattr(node, 'all_names') and len(node.all_names) > 1:
                if parent:
                    dx = node.unrooted_x - parent.unrooted_x
                    dy = node.unrooted_y - parent.unrooted_y
                    angle = math.atan2(dy, dx)
                    
                    vertical_height = 15 * (len(node.all_names) - 1)
                    
                    fig.add_trace(go.Scatter(
                        x=[node.unrooted_x, node.unrooted_x],
                        y=[node.unrooted_y - vertical_height/2, node.unrooted_y + vertical_height/2],
                        mode='lines',
                        line=dict(color='rgb(40,40,40)', width=2),
                        hoverinfo='skip',
                        showlegend=False
                    ))
                    
                    if -math.pi/2 <= angle <= math.pi/2:
                        xanchor = 'left'
                        xshift = 5
                    else:
                        xanchor = 'right'
                        xshift = -5
                    
                    for i, name in enumerate(sorted(node.all_names)):
                        y_offset = node.unrooted_y + (i - (len(node.all_names)-1)/2) * 15
                        hover_text = name
                        if sample_coverage_dict and name in sample_coverage_dict:
                            hover_text = "{name}<br>Coverage: {coverage}X".format(
                                name=name,
                                coverage=sample_coverage_dict[name]
                            )
                        
                        fig.add_annotation(
                            x=node.unrooted_x,
                            y=y_offset,
                            text=name,
                            showarrow=False,
                            font=dict(size=11),
                            xanchor=xanchor,
                            yanchor='middle',
                            xshift=xshift,
                            bgcolor='rgba(255,255,255,0.8)',
                            borderpad=2,
                            hovertext=hover_text,
                            hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace"),
                        )
            else:
                if parent:
                    dx = node.unrooted_x - parent.unrooted_x
                    dy = node.unrooted_y - parent.unrooted_y
                    angle = math.atan2(dy, dx)
                    base_offset = label_offset_base
                    
                    label_x = node.unrooted_x + math.cos(angle) * base_offset
                    label_y = node.unrooted_y + math.sin(angle) * base_offset
                    
                    if -math.pi/2 <= angle <= math.pi/2:
                        xanchor = 'left'
                        xshift = 5
                    else:
                        xanchor = 'right'
                        xshift = -5
                else:
                    label_x = node.unrooted_x
                    label_y = node.unrooted_y
                    xanchor = 'center'
                    xshift = 0
                
                hover_text = node.name
                if sample_coverage_dict and node.name in sample_coverage_dict:
                    hover_text = "{name}<br>Coverage: {coverage}X".format(
                        name=node.name,
                        coverage=sample_coverage_dict[node.name]
                    )
                
                fig.add_annotation(
                    x=label_x,
                    y=label_y,
                    text=node.name,
                    showarrow=False,
                    font=dict(size=11),
                    xanchor=xanchor,
                    yanchor='middle',
                    xshift=xshift,
                    bgcolor='rgba(255,255,255,0.8)',
                    borderpad=2,
                    hovertext=hover_text,
                    hoverlabel=dict(bgcolor="white", font_size=12, font_family="monospace"),
                )
    
    for child in node.children:
        add_unrooted_traces_improved(fig, child, node, label_offset_base, processed_positions, sample_coverage_dict)

def read_table(input_file):
    '''
    The SNP table workbook as html_tree reads it: pd.read_excel, except that
    the first column, the sample names, is read without the default NA
    strings.  'NA', 'None', 'null', 'nan' and 'N/A' are names there; read as
    missing, such a sample was labelled "nan" in the tree and two of them
    shared one clade key, which showed as support above 100%.  The other
    columns are read exactly as before.
    '''
    df = pd.read_excel(input_file, engine='openpyxl')
    names = pd.read_excel(input_file, engine='openpyxl', keep_default_na=False)
    df[df.columns[0]] = names.iloc[:, 0].to_numpy(dtype=object)
    return df


def html_tree(file_path=None, sample_coverage_dict=None, bootstrap_seed=BOOTSTRAP_SEED, df=None,
              raise_errors=False):
    """
    Main function to create phylogenetic tree visualization with coverage information

    df: the table as read_table(file_path) would return it.  Step 2 passes
    the frame it has just written (Tables.table_to_tree_df) so the workbook
    is not read back through openpyxl, which took longer than everything
    else in this function once the tree builder was fixed: 5.6 s for a
    400-sample table.  file_path still names the output.

    raise_errors: let an exception propagate instead of printing it and
    returning 1.  Step 2 sets it so the failure, and its message, reach the
    process that records it.
    """
    if file_path:
        input_file = file_path
        n_bootstrap = 100
    else:
        args = parse_args()
        input_file = args.file
        n_bootstrap = args.bootstrap
        bootstrap_seed = args.bootstrap_seed

    try:
        if df is None:
            df = read_table(input_file)
        input_name = input_file.rsplit('.', 1)[0]
        output_file = "{input_name}_position_tree.html".format(input_name=input_name)
        create_tree_visualization(df, output_file, input_file, sample_coverage_dict,
                                  bootstrap_seed=bootstrap_seed)
        return 0

    except Exception as e:
        if raise_errors:
            raise
        print("Error: {error}".format(error=e))
        import traceback
        traceback.print_exc()  # Added traceback for better error reporting
        return 1
    
if __name__ == "__main__":
    exit(html_tree())