#!/usr/bin/env python3
"""
Compare two step2 snapshots at the level that matters scientifically.

A byte diff is the wrong gate.  Sorting the VCF glob made sample order
deterministic, which legitimately reorders FASTA records and cascade-table
columns without changing a single call.  So compare:

  - SNP tables  as {(sheet, row label, position) -> value}
  - FASTA       as {sample -> sequence}
  - newick      as the set of leaf names, plus the raw string
  - bootstrap   as the ordered list of support values

usage: compare.py <snapshot_a> <snapshot_b>
exit 0 when the data agrees, 1 when it does not.
"""
import os
import re
import sys


def table_cells(path):
    sheets, cur = {}, None
    for line in open(path).read().splitlines():
        if line.startswith('### sheet:'):
            cur = line.split()[2]
            sheets[cur] = []
        elif cur is not None:
            sheets[cur].append(line.split('\t'))
    out = {}
    for sheet, rows in sheets.items():
        if not rows:
            continue
        positions = rows[0][1:]
        for row in rows[1:]:
            for pos, val in zip(positions, row[1:]):
                out[(sheet, row[0], pos)] = val
    return out


def fasta_records(path):
    recs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
        elif name:
            recs[name] = recs.get(name, '') + line
    return recs


def newick_leaves(path):
    body = open(path).read()
    return set(re.findall(r'[(,]([^(),:]+):', body))


def main():
    a_root, b_root = sys.argv[1], sys.argv[2]
    names = sorted(set(os.listdir(a_root)) | set(os.listdir(b_root)))
    problems, ordering = [], []

    for name in names:
        a, b = os.path.join(a_root, name), os.path.join(b_root, name)
        if name == 'MANIFEST':
            continue
        if not os.path.exists(a) or not os.path.exists(b):
            problems.append(f'{name}: present in only one snapshot')
            continue

        if name.endswith('.tsv'):
            ca, cb = table_cells(a), table_cells(b)
            only = set(ca) ^ set(cb)
            diff = {k for k in set(ca) & set(cb) if ca[k] != cb[k]}
            if only:
                problems.append(f'{name}: {len(only)} cells present in only one side'
                                + ''.join(f'\n      {k}' for k in sorted(only)[:5]))
            if diff:
                problems.append(f'{name}: {len(diff)} cells changed'
                                + ''.join(f'\n      {k}: {ca[k]!r} -> {cb[k]!r}'
                                          for k in sorted(diff)[:8]))
            if not only and not diff and open(a).read() != open(b).read():
                ordering.append(f'{name}: same cells, different row/column order')

        elif name.endswith('.fasta'):
            fa, fb = fasta_records(a), fasta_records(b)
            if fa != fb:
                for k in sorted(set(fa) | set(fb)):
                    if fa.get(k) != fb.get(k):
                        problems.append(f'{name}: {k}: {fa.get(k)} -> {fb.get(k)}')
            elif list(fa) != list(fb):
                ordering.append(f'{name}: same records, different order')

        elif name.endswith('.tre'):
            la, lb = newick_leaves(a), newick_leaves(b)
            if la != lb:
                problems.append(f'{name}: leaf set differs: {la ^ lb}')
            elif open(a).read() != open(b).read():
                ordering.append(f'{name}: same leaves, different newick string')

        elif name.endswith('.bootstrap'):
            va, vb = open(a).read().split(), open(b).read().split()
            if va != vb:
                problems.append(f'{name}: bootstrap support differs: {va} -> {vb}')

        elif name.endswith('.ERROR'):
            problems.append(f'{name}: snapshot recorded an error')

    for line in ordering:
        print(f'  order  {line}')
    for line in problems:
        print(f'  DATA   {line}')
    print()
    if problems:
        print(f'VERDICT: {len(problems)} data difference(s)')
        return 1
    print('VERDICT: data identical'
          + (f' ({len(ordering)} presentation-order difference(s))' if ordering else ''))
    return 0


if __name__ == '__main__':
    sys.exit(main())
