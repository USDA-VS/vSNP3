#!/usr/bin/env python3
"""
Benchmark the "FASTAs, RAxML trees, tables and HTML trees" phase of vSNP3 step 2
on a synthetic group, sub-step by sub-step, against a chosen bin/ directory.

  bench_tree_phase.py gen  --out DIR --samples S --positions P [--seed N]
  bench_tree_phase.py run  --bin BIN --data DIR --out DIR [--raxml] [--html]
                           [--html_from_df] [--profile] [--show_groups] [--dp]

`gen` writes the same inputs raxml_table_build() sees: the group alignment frame
(a random clonal tree with a few per cent of no-calls, N and IUPAC codes), the
FASTA, the annotation and MQ frames, a coverage dict and a groupings dict.
`run` imports Tree/Tables/html_tree from BIN and times each sub-step; it also
checks that the frame Tables hands to html_tree equals what html_tree would read
from the workbook.  RAxML is run once per data dir and the .tre cached, so table
and HTML timings for two BIN versions (say, a pristine copy and a working copy)
use the same tree.  Compare two output directories with the regression tools:
convert each with snapshot.py and run compare.py, or diff the HTML after
normalizing plotly's uuid and the run timestamp.

Typical use, from the repository root:

  python3 tests/regression/bench_tree_phase.py gen --out /tmp/g400 --samples 400 --positions 4000
  python3 tests/regression/bench_tree_phase.py run --bin bin --data /tmp/g400 --out /tmp/g400_out \
      --raxml --html --html_from_df
"""
import argparse
import cProfile
import importlib
import json
import os
import pickle
import pstats
import shutil
import sys
import time
from collections import defaultdict

import numpy as np
import pandas as pd

CHROM = 'NC_002945.4'
GENOME = 4_349_904


def gen(args):
    rng = np.random.default_rng(args.seed)
    S, P = args.samples, args.positions
    os.makedirs(args.out, exist_ok=True)
    names = [f'S{i:04d}' for i in range(S)]

    # Random bifurcating tree over the samples: each internal node is an edge
    # that carries mutations shared by every leaf below it.
    nodes = [[n] for n in names]            # leaf sets
    edges = [list(x) for x in nodes]        # terminal edges first
    while len(nodes) > 1:
        i, j = sorted(rng.choice(len(nodes), size=2, replace=False))
        merged = nodes[i] + nodes[j]
        edges.append(merged)
        nodes = [n for k, n in enumerate(nodes) if k not in (i, j)] + [merged]

    positions = np.sort(rng.choice(np.arange(1, GENOME), size=P, replace=False))
    abs_pos = [f'{CHROM}:{p}' for p in positions]
    bases = np.array(list('ACGT'))
    root = rng.choice(bases, size=P)

    # Every position mutates on exactly one edge (clonal, no homoplasy).  Weight
    # internal edges by subtree size so deep edges carry more SNPs, like real data.
    weights = np.array([len(e) for e in edges], dtype=float)
    weights[:S] = 1.5                       # terminal edges: a few private SNPs each
    weights /= weights.sum()
    edge_of_pos = rng.choice(len(edges), size=P, p=weights)
    alt = np.array([rng.choice([b for b in 'ACGT' if b != r]) for r in root])

    idx = {n: k for k, n in enumerate(names)}
    mat = np.tile(root, (S, 1)).astype('<U1')
    for p in range(P):
        for leaf in edges[edge_of_pos[p]]:
            mat[idx[leaf], p] = alt[p]

    # Noise: no coverage, N, mixed calls.
    noise = rng.random((S, P))
    mat[noise < 0.010] = '-'
    mat[(noise >= 0.010) & (noise < 0.013)] = 'N'
    iupac = np.array(list('RYSWKM'))
    mixed = (noise >= 0.013) & (noise < 0.016)
    mat[mixed] = rng.choice(iupac, size=int(mixed.sum()))

    group_df = pd.DataFrame(mat, index=names, columns=abs_pos)
    group_df.loc['root'] = root            # dict_to_dataframe puts root last
    group_df = group_df.astype(object)
    group_df.columns.name = 'abs_pos'      # as dict_to_dataframe leaves it
    group_df.to_pickle(os.path.join(args.out, 'group_df.pkl'))

    with open(os.path.join(args.out, 'bench-group.fasta'), 'w') as fh:
        for name, row in group_df.iterrows():
            print(f'>{name}', file=fh)
            print(''.join(row.values), file=fh)

    # Annotation and MQ frames cover more positions than the group (run-wide).
    extra = [f'{CHROM}:{p}' for p in rng.choice(np.arange(1, GENOME), size=P // 2, replace=False)]
    all_pos = sorted(set(abs_pos) | set(extra), key=lambda s: int(s.split(':')[1]))
    ann = pd.DataFrame({'abs_pos': all_pos,
                        'annotation': [f'A->T, gene{k % 500}:L{k}F, hypothetical protein, nonsynonymous, codon_pos_{k % 3 + 1}'
                                       if k % 7 else 'position not annotated'
                                       for k in range(len(all_pos))]})
    ann.to_pickle(os.path.join(args.out, 'annotation.pkl'))
    mq = pd.DataFrame({'abs_pos': all_pos, 'MQ': rng.uniform(40, 60, size=len(all_pos))})
    mq.to_pickle(os.path.join(args.out, 'mq.pkl'))
    cov = {n: int(rng.integers(20, 200)) for n in names}
    with open(os.path.join(args.out, 'coverage.json'), 'w') as fh:
        json.dump(cov, fh)
    groupings = {'bench-group': {n: pd.DataFrame() for n in names},
                 'bench-sub': {n: pd.DataFrame() for n in names[: S // 3]}}
    with open(os.path.join(args.out, 'groupings.pkl'), 'wb') as fh:
        pickle.dump(groupings, fh)
    print(f'generated {S} samples x {P} positions -> {args.out}')


class Timer:
    def __init__(self):
        self.t = defaultdict(float)
        self.n = defaultdict(int)

    def wrap(self, obj, attr, label):
        orig = getattr(obj, attr)
        timer = self

        depth = [0]

        def wrapped(*a, **k):
            depth[0] += 1
            t0 = time.perf_counter()
            try:
                return orig(*a, **k)
            finally:
                depth[0] -= 1
                if depth[0] == 0:          # outermost call only (recursive functions)
                    timer.t[label] += time.perf_counter() - t0
                    timer.n[label] += 1
        setattr(obj, attr, wrapped)
        return orig


def run(args):
    sys.path.insert(0, os.path.abspath(args.bin))
    for m in list(sys.modules):
        if m.startswith('vsnp3'):
            del sys.modules[m]
    fts = importlib.import_module('vsnp3_fasta_to_snps_table')
    vht = importlib.import_module('vsnp3_html_tree')
    assert fts.__file__.startswith(os.path.abspath(args.bin)), fts.__file__

    data = args.data
    out = os.path.abspath(args.out)
    if os.path.exists(out):
        shutil.rmtree(out)
    os.makedirs(out)
    group = 'bench-group'
    group_df = pd.read_pickle(os.path.join(data, 'group_df.pkl'))
    ann = pd.read_pickle(os.path.join(data, 'annotation.pkl'))
    mq = pd.read_pickle(os.path.join(data, 'mq.pkl'))
    cov = json.load(open(os.path.join(data, 'coverage.json')))
    groupings = pickle.load(open(os.path.join(data, 'groupings.pkl'), 'rb'))
    fasta = os.path.join(out, f'{group}-bench.fasta')
    shutil.copy(os.path.join(data, 'bench-group.fasta'), fasta)

    timer = Timer()
    times = {}

    # --- RAxML -------------------------------------------------------------
    cached_tree = os.path.join(data, 'cached.tre')
    if args.raxml or not os.path.exists(cached_tree):
        t0 = time.perf_counter()
        tree = fts.Tree(fasta_alignments=fasta, write_path=out, tree_name=group)
        times['raxml'] = time.perf_counter() - t0
        shutil.copy(tree.newick, cached_tree)
        newick = tree.newick
    else:
        newick = os.path.join(out, f'{group}_cached.tre')
        shutil.copy(cached_tree, newick)
        times['raxml'] = None

    # --- Tables ------------------------------------------------------------
    timer.wrap(fts.Tables, 'write_out_table', 'tables.write_out_table')
    timer.wrap(fts.Tables, 'excel_formatter', 'tables.excel_formatter')
    tables = fts.Tables(df_alignments=group_df, tree=newick, gbk=ann, mq=mq,
                        write_path=out, groupings_dict=groupings,
                        show_groups=args.show_groups, table_name=group, debug=False,
                        sample_coverage_dict=cov if args.dp else {})
    prof = cProfile.Profile() if args.profile else None
    t0 = time.perf_counter()
    if prof:
        prof.runcall(tables.build_tables)
    else:
        tables.build_tables()
    times['tables.build_tables'] = time.perf_counter() - t0
    if prof:
        prof.dump_stats(os.path.join(out, 'tables.prof'))
    times['tables.write_out_table'] = timer.t['tables.write_out_table']
    # Does the frame handed to html_tree equal what html_tree would read from
    # the workbook?  read_table where BIN has it (it keeps NA-like sample
    # names), plain read_excel for an older BIN.
    tree_df = getattr(tables, 'table_to_tree_df', None)
    if tree_df is not None:
        t0 = time.perf_counter()
        if hasattr(vht, 'read_table'):
            read_back = vht.read_table(tables.table_to_tree_path)
        else:
            read_back = pd.read_excel(tables.table_to_tree_path, engine='openpyxl')
        times['check.read_excel'] = time.perf_counter() - t0
        problems = []
        if list(read_back.columns) != list(tree_df.columns):
            problems.append(f'columns differ: {list(read_back.columns)[:3]} vs {list(tree_df.columns)[:3]}')
        if read_back.shape != tree_df.shape:
            problems.append(f'shape {read_back.shape} vs {tree_df.shape}')
        else:
            dt = [(c, str(read_back[c].dtype), str(tree_df[c].dtype)) for c in read_back.columns
                  if str(read_back[c].dtype) != str(tree_df[c].dtype)]
            if dt:
                problems.append(f'{len(dt)} dtype differences, e.g. {dt[:3]}')
            a = read_back.to_numpy(dtype=object); b = tree_df.to_numpy(dtype=object)
            neq = ~((a == b) | (pd.isna(a) & pd.isna(b)))
            # the MQ row is read back as int where integral; html_tree never reads it
            mq_rows = (read_back.iloc[:, 0] == 'MQ').to_numpy()
            neq[mq_rows, :] = False
            n = int(neq.sum())
            if n:
                r, c = [x[0] for x in neq.nonzero()]
                problems.append(f'{n} cells differ, e.g. [{r},{c}] {a[r, c]!r} vs {b[r, c]!r}')
            # types of the cells html_tree reads: first column and sample rows
            ta = {type(x).__name__ for x in a[:, 0]}; tb = {type(x).__name__ for x in b[:, 0]}
            if ta != tb:
                problems.append(f'first-column cell types {ta} vs {tb}')
        times['check.tree_df_matches_read_excel'] = 'yes' if not problems else 'NO: ' + '; '.join(problems)
    if args.html and tree_df is not None and args.html_from_df:
        html_df = tree_df
    else:
        html_df = None
    times['tables.excel_formatter'] = timer.t['tables.excel_formatter']
    times['tables.python_loops'] = times['tables.build_tables'] - timer.t['tables.write_out_table']

    # --- HTML tree ---------------------------------------------------------
    if args.html:
        timer.wrap(vht, 'build_tree_rectangular', 'html.build_tree_rectangular')
        timer.wrap(vht, 'build_tree_unrooted', 'html.build_tree_unrooted')
        timer.wrap(vht, 'add_rectangular_traces', 'html.add_rectangular_traces')
        timer.wrap(vht, 'add_unrooted_traces_improved', 'html.add_unrooted_traces')
        timer.wrap(vht, 'get_branch_configs', 'html.get_branch_configs')
        timer.wrap(vht, 'add_bootstrap_to_node', 'html.add_bootstrap_to_node')
        timer.wrap(vht, 'assign_coordinates', 'html.assign_coordinates')
        timer.wrap(vht, 'assign_unrooted_coordinates_improved', 'html.assign_unrooted_coordinates')
        timer.wrap(vht.go.Figure, 'to_html', 'html.to_html')
        timer.wrap(vht.pd, 'read_excel', 'html.read_excel')
        for extra in ('_bootstrap_support', 'bootstrap_support', 'calculate_bootstrap_support'):
            if hasattr(vht, extra):
                timer.wrap(vht, extra, 'html.' + extra)
        html_kwargs = {'df': html_df} if html_df is not None else {}
        prof = cProfile.Profile() if args.profile else None
        t0 = time.perf_counter()
        if prof:
            rc = prof.runcall(vht.html_tree, tables.table_to_tree_path, **html_kwargs)
        else:
            rc = vht.html_tree(tables.table_to_tree_path, **html_kwargs)
        times['html.total'] = time.perf_counter() - t0
        if prof:
            prof.dump_stats(os.path.join(out, 'html.prof'))
        assert rc == 0, 'html_tree returned nonzero'
        for k, v in timer.t.items():
            if k.startswith('html.'):
                times[k] = v
        times['html.build_tree_rectangular.calls'] = timer.n['html.build_tree_rectangular']

    with open(os.path.join(out, 'times.json'), 'w') as fh:
        json.dump(times, fh, indent=1)
    width = max(len(k) for k in times)
    print(f'--- {args.bin}  data={data}')
    for k, v in times.items():
        if isinstance(v, float):
            print(f'  {k:<{width}}  {v:8.2f} s')
        else:
            print(f'  {k:<{width}}  {v}')
    if args.profile:
        for name in ('tables.prof', 'html.prof'):
            p = os.path.join(out, name)
            if os.path.exists(p):
                print(f'\n=== {name} top 25 by cumulative')
                pstats.Stats(p).sort_stats('cumulative').print_stats(25)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    g = sub.add_parser('gen')
    g.add_argument('--out', required=True)
    g.add_argument('--samples', type=int, required=True)
    g.add_argument('--positions', type=int, required=True)
    g.add_argument('--seed', type=int, default=1)
    r = sub.add_parser('run')
    r.add_argument('--bin', required=True)
    r.add_argument('--data', required=True)
    r.add_argument('--out', required=True)
    r.add_argument('--raxml', action='store_true')
    r.add_argument('--html', action='store_true')
    r.add_argument('--profile', action='store_true')
    r.add_argument('--show_groups', action='store_true')
    r.add_argument('--dp', action='store_true')
    r.add_argument('--html_from_df', action='store_true', help='hand the in-memory frame to html_tree')
    args = ap.parse_args()
    {'gen': gen, 'run': run}[args.cmd](args)


if __name__ == '__main__':
    main()
