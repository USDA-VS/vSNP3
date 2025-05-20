#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from collections import defaultdict
import os
import math

__version__ = "3.30"

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

def parse_args():
    parser = argparse.ArgumentParser(description='Create interactive phylogenetic tree from SNP alignment')
    parser.add_argument('-f', '--file', required=True, help='Input Excel file containing SNP alignment')
    parser.add_argument('-o', '--output', default='phylo_tree.html', 
                       help='Output HTML file name (default: phylo_tree.html)')
    parser.add_argument('-b', '--bootstrap', type=int, default=100,
                       help='Number of bootstrap replicates (default: 100)')
    return parser.parse_args()

def get_snp_profile(seq, root_seq, positions):
    """Get SNP profile for a sequence compared to root"""
    snps = []
    for i, (base, root_base) in enumerate(zip(seq, root_seq)):
        if base != root_base and base != '-':
            snps.append(f"{positions[i]}:{root_base}->{base}")
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

def calculate_bootstrap_support(df, n_replicates=100):
    """Calculate bootstrap support values for tree branches"""
    sequences = df.iloc[:, 1:]
    positions = sequences.columns.tolist()
    
    # Store all branch configurations
    branch_configs = defaultdict(int)
    
    for _ in range(n_replicates):
        # Randomly sample positions with replacement
        bootstrap_positions = np.random.choice(positions, size=len(positions), replace=True)
        bootstrap_df = pd.concat([df.iloc[:, 0], df[bootstrap_positions]], axis=1)
        
        # Build bootstrap tree
        bootstrap_root = build_tree(bootstrap_df)
        
        # Get branch configurations from bootstrap tree
        get_branch_configs(bootstrap_root, branch_configs)
    
    # Convert counts to percentages
    for config in branch_configs:
        branch_configs[config] = (branch_configs[config] / n_replicates) * 100
        
    return branch_configs

def get_branch_configs(node, configs):
    """Extract branch configurations from tree"""
    if node.children:
        # Sort leaf names to ensure consistent configuration strings
        leaf_names = sorted([leaf.name for leaf in get_leaves(node)])
        if len(leaf_names) > 1:  # Only count non-trivial branches
            config = '|'.join(leaf_names)
            configs[config] += 1
        
        for child in node.children:
            get_branch_configs(child, configs)

def get_leaves(node):
    """Get all leaf nodes in the subtree"""
    if not node.children:
        return [node]
    leaves = []
    for child in node.children:
        leaves.extend(get_leaves(child))
    return leaves

def add_bootstrap_to_node(node, bootstrap_values):
    """Add bootstrap support values to tree nodes"""
    if node.children:
        leaf_names = sorted([leaf.name for leaf in get_leaves(node)])
        if len(leaf_names) > 1:
            config = '|'.join(leaf_names)
            node.bootstrap = bootstrap_values.get(config, 0)
        for child in node.children:
            add_bootstrap_to_node(child, bootstrap_values)

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
    """Count the number of leaf nodes in the subtree"""
    if not node.children:
        return 1
    return sum(count_leaves(child) for child in node.children)

def create_hover_text(snps, bootstrap=None):
    """Create hover text for branches including bootstrap values"""
    text = f"SNPs ({len(snps)}):<br>" + "<br>".join(snps)
    if bootstrap is not None:
        text = f"Bootstrap: {bootstrap:.1f}%<br><br>" + text
    return text

def create_hover_points(snps, bootstrap=None, chunk_size=30):
    """Split SNPs into chunks for multiple hover points"""
    chunks = []
    for i in range(0, len(snps), chunk_size):
        chunk = snps[i:i + chunk_size]
        hover_text = f"SNPs {i+1}-{min(i+chunk_size, len(snps))} of {len(snps)}:<br>" + "<br>".join(chunk)
        if bootstrap is not None:
            hover_text = f"Bootstrap: {bootstrap:.1f}%<br><br>" + hover_text
        chunks.append(hover_text)
    return chunks

def add_bootstrap_annotation(fig, x, y, bootstrap, offset=10):
    """Add bootstrap value annotation"""
    if bootstrap is not None:
        fig.add_annotation(
            x=x,
            y=y + offset,
            text=f"{bootstrap:.0f}",
            showarrow=False,
            font=dict(size=10),
            yanchor='bottom'
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
    # Drop non-SNP columns
    df_clean = df.copy()
    
    # Remove 'Average Coverage' column if it exists
    if 'Average Coverage' in df_clean.columns:
        df_clean = df_clean.drop('Average Coverage', axis=1)
    
    # Keep first column (names) and all NC columns
    cols_to_keep = [df_clean.columns[0]] + [col for col in df_clean.columns if col.startswith('NC_')]
    df_clean = df_clean[cols_to_keep]
    
    # Get sequences
    sequences = df_clean.iloc[:, 1:]
    positions = sequences.columns.tolist()
    root_sequence = sequences.iloc[0].tolist()
    
    # Get samples (excluding root, MQ, and annotation)
    samples = df_clean[~df_clean.iloc[:, 0].isin(["MQ", "annotation"])].iloc[1:]
    
    # Create leaf nodes for each sample
    leaf_nodes = []
    for idx, row in samples.iterrows():
        node = Node(row.iloc[0])
        node.snps = get_snp_profile(row.iloc[1:], root_sequence, positions)
        node.excel_index = idx
        node.leaf_indexes = [idx]
        leaf_nodes.append(node)
    
    # Create root node
    root = Node("Reference")
    
    def get_shared_snps(nodes):
        if not nodes:
            return set()
        snp_sets = [set(node.snps) for node in nodes]
        return set.intersection(*snp_sets)
    
    # Build tree bottom-up based on shared SNPs
    nodes = leaf_nodes.copy()
    while len(nodes) > 1:
        max_shared = 0
        best_pair = None
        
        for i, node1 in enumerate(nodes):
            for j, node2 in enumerate(nodes[i+1:], i+1):
                shared = get_shared_snps([node1, node2])
                if len(shared) > max_shared:
                    max_shared = len(shared)
                    best_pair = (node1, node2)
        
        if not best_pair or max_shared == 0:
            for node in nodes:
                root.children.append(node)
                node.parent = root
            break
        
        node1, node2 = best_pair
        internal = Node()
        internal.snps = sorted(get_shared_snps([node1, node2]))
        internal.leaf_indexes = sorted(node1.leaf_indexes + node2.leaf_indexes)
        
        shared_snps = set(internal.snps)
        node1.snps = sorted(set(node1.snps) - shared_snps)
        node2.snps = sorted(set(node2.snps) - shared_snps)
        
        internal.children = [node1, node2]
        node1.parent = internal
        node2.parent = internal
        
        nodes.remove(node1)
        nodes.remove(node2)
        nodes.append(internal)
    
    if len(nodes) == 1:
        root.children.append(nodes[0])
        nodes[0].parent = root
    
    return root

def build_tree_unrooted(df):
    """Build phylogenetic tree optimized for unrooted layout"""
    # Drop non-SNP columns
    df_clean = df.copy()
    
    # Remove 'Average Coverage' column if it exists
    if 'Average Coverage' in df_clean.columns:
        df_clean = df_clean.drop('Average Coverage', axis=1)
    
    # Keep first column (names) and all NC columns
    cols_to_keep = [df_clean.columns[0]] + [col for col in df_clean.columns if col.startswith('NC_')]
    df_clean = df_clean[cols_to_keep]
    
    # Get sequences
    sequences = df_clean.iloc[:, 1:]
    positions = sequences.columns.tolist()
    root_sequence = sequences.iloc[0].tolist()
    
    # Get samples (excluding root, MQ, and annotation)
    samples = df_clean[~df_clean.iloc[:, 0].isin(["MQ", "annotation"])].iloc[1:]
    
    # First, collect identical sequence information
    sequence_groups = defaultdict(list)
    for idx, row in samples.iterrows():
        seq_key = tuple(row.iloc[1:])
        sequence_groups[seq_key].append((idx, row.iloc[0]))
    
    # Create leaf nodes
    leaf_nodes = []
    processed_seqs = set()
    
    for idx, row in samples.iterrows():
        seq_key = tuple(row.iloc[1:])
        if seq_key not in processed_seqs:
            node = Node(row.iloc[0])
            node.snps = get_snp_profile(row.iloc[1:], root_sequence, positions)
            node.excel_index = idx
            node.leaf_indexes = [idx]
            
            # Handle identical sequences
            identical_group = sequence_groups[seq_key]
            if len(identical_group) > 1:
                node.all_names = sorted(name for _, name in identical_group)
            
            leaf_nodes.append(node)
            processed_seqs.add(seq_key)
    
    # Create root node
    root = Node("Reference")
    
    def get_shared_snps(nodes):
        if not nodes:
            return set()
        snp_sets = [set(node.snps) for node in nodes]
        return set.intersection(*snp_sets)
    
    # Build tree bottom-up based on shared SNPs
    nodes = leaf_nodes.copy()
    while len(nodes) > 1:
        max_shared = 0
        best_pair = None
        
        for i, node1 in enumerate(nodes):
            for j, node2 in enumerate(nodes[i+1:], i+1):
                shared = get_shared_snps([node1, node2])
                if len(shared) > max_shared:
                    max_shared = len(shared)
                    best_pair = (node1, node2)
        
        if not best_pair or max_shared == 0:
            for node in nodes:
                root.children.append(node)
                node.parent = root
            break
        
        node1, node2 = best_pair
        internal = Node()
        internal.snps = sorted(get_shared_snps([node1, node2]))
        internal.leaf_indexes = sorted(node1.leaf_indexes + node2.leaf_indexes)
        
        shared_snps = set(internal.snps)
        node1.snps = sorted(set(node1.snps) - shared_snps)
        node2.snps = sorted(set(node2.snps) - shared_snps)
        
        internal.children = [node1, node2]
        node1.parent = internal
        node2.parent = internal
        
        nodes.remove(node1)
        nodes.remove(node2)
        nodes.append(internal)
    
    if len(nodes) == 1:
        root.children.append(nodes[0])
        nodes[0].parent = root
    
    return root

def create_tree_visualization(df, output_file, input_filename, sample_coverage_dict=None):
    """Create interactive tree visualization using plotly"""
    df_for_tree = df.copy()
    coverage_dict = {}

    if sample_coverage_dict:
        coverage_dict = sample_coverage_dict
    elif 'Average Coverage' in df_for_tree.columns:
        try:
            mask = ~df_for_tree.iloc[:, 0].isin(["Reference", "MQ", "annotation"])
            sample_names = df_for_tree.loc[mask, df_for_tree.columns[0]]
            coverage_values = df_for_tree.loc[mask, 'Average Coverage']
            coverage_dict = dict(zip(sample_names, coverage_values))
            df_for_tree = df_for_tree.drop('Average Coverage', axis=1)
        except Exception as e:
            print(f"Warning: Could not extract coverage data: {e}")
            coverage_dict = {}
    
    rect_root = build_tree_rectangular(df_for_tree)
    unrooted_root = build_tree_unrooted(df_for_tree)
    
    sequences = df_for_tree.iloc[:, 1:]
    positions = sequences.columns.tolist()
    
    def build_bootstrap_tree(bootstrap_df):
        """Helper function to build a tree for bootstrap calculations"""
        return build_tree_rectangular(bootstrap_df)
    
    bootstrap_values = defaultdict(int)
    n_replicates = 100
    
    for _ in range(n_replicates):
        bootstrap_positions = np.random.choice(positions, size=len(positions), replace=True)
        bootstrap_df = pd.concat([df_for_tree.iloc[:, 0], df_for_tree[bootstrap_positions]], axis=1)
        bootstrap_root = build_bootstrap_tree(bootstrap_df)
        get_branch_configs(bootstrap_root, bootstrap_values)
    
    for config in bootstrap_values:
        bootstrap_values[config] = (bootstrap_values[config] / n_replicates) * 100
    
    add_bootstrap_to_node(rect_root, bootstrap_values)
    add_bootstrap_to_node(unrooted_root, bootstrap_values)
    
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
    
    # Add refined scale labels
    fig_rect.add_annotation(
        x=scale_x - 5,
        y=scale_y,
        text="0",
        showarrow=False,
        font=dict(size=9),
        xanchor='right',
        yanchor='middle'
    )
    fig_rect.add_annotation(
        x=scale_x + scale_length * snp_scale + 5,
        y=scale_y,
        text=f"{scale_length} SNPs",
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    
    # Add total SNPs annotation
    fig_rect.add_annotation(
        x=scale_x,
        y=scale_y - 25,
        text=f"Total SNPs in tree: {total_snps}",
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    
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
    
    # Add refined scale labels
    fig_unrooted.add_annotation(
        x=unrooted_scale_x - 5,
        y=unrooted_scale_y,
        text="0",
        showarrow=False,
        font=dict(size=9),
        xanchor='right',
        yanchor='middle'
    )
    fig_unrooted.add_annotation(
        x=unrooted_scale_x + scale_length * unrooted_snp_scale + 5,
        y=unrooted_scale_y,
        text=f"{scale_length} SNPs",
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    
    # Add total SNPs annotation for unrooted tree
    fig_unrooted.add_annotation(
        x=unrooted_scale_x,
        y=unrooted_scale_y - 25,
        text=f"Total SNPs in tree: {total_snps}",
        showarrow=False,
        font=dict(size=9),
        xanchor='left',
        yanchor='middle'
    )
    
    # Get base filename for title
    base_filename = os.path.splitext(os.path.basename(input_filename))[0]
    
    # Configure layouts
    layout_common = dict(
        showlegend=False,
        plot_bgcolor='white',
        margin=dict(l=100, r=100, t=100, b=100),
        hovermode='closest'
    )
    
    title_text = (f"{base_filename}<br>"
                 "Hover over branches/dots to see SNPs and bootstrap values<br>"
                 "Hover over sample names to see coverage information")
    
    fig_rect.update_layout(
        title=dict(text=title_text, x=0.5),
        width=1200,
        height=800,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        **layout_common
    )
    
    # Create HTML with toggle button
    html_content = f"""
    <html>
    <head>
        <title>Phylogenetic Tree Visualization</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
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
        </style>
    </head>
    <body>
        <div style="text-align: center;">
            <button class="button" onclick="toggleView()" id="toggleButton">Unrooted</button>
        </div>
        <div id="rectangular-tree">{fig_rect.to_html(full_html=False)}</div>
        <div id="unrooted-tree" style="display: none;">{fig_unrooted.to_html(full_html=False)}</div>
        
        <script>
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
        </script>
    </body>
    </html>
    """
    
    # Save as interactive HTML
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

def add_rectangular_traces(fig, node, parent=None, sample_coverage_dict=None):
    """Add traces for rectangular layout with bootstrap values and coverage information"""
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
            hover_text = f"{node.name}<br>Coverage: {sample_coverage_dict[node.name]}X"
        
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
                    text=f"{node.bootstrap:.0f}",
                    showarrow=False,
                    font=dict(size=10),
                    bgcolor='rgba(255,255,255,0.8)'
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
                            hover_text = f"{name}<br>Coverage: {sample_coverage_dict[name]}X"
                        
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
                    hover_text = f"{node.name}<br>Coverage: {sample_coverage_dict[node.name]}X"
                
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

def html_tree(file_path=None, sample_coverage_dict=None):
    """Main function to create phylogenetic tree visualization with coverage information"""
    if file_path:
        input_file = file_path
        n_bootstrap = 100
    else:
        args = parse_args()
        input_file = args.file
        n_bootstrap = args.bootstrap
    
    try:
        df = pd.read_excel(input_file, engine='openpyxl')  # Updated to specify engine
        input_name = input_file.rsplit('.', 1)[0]
        output_file = f"{input_name}_position_tree.html"
        create_tree_visualization(df, output_file, input_file, sample_coverage_dict)
        return 0
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()  # Added traceback for better error reporting
        return 1
    
if __name__ == "__main__":
    exit(html_tree())
    
    fig_unrooted.update_layout(
        title=dict(text=f"{title_text} - Unrooted Layout", x=0.5),
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