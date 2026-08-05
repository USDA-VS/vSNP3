#!/usr/bin/env python3
"""
Normalize a vSNP3 step2 output tree into a diffable snapshot.

Output filenames and in-file content carry run timestamps, so a raw `diff -r`
between two runs is useless.  This writes one text file per artifact with the
timestamp stripped from the name, Excel converted to TSV, and newick/fasta
passed through verbatim.  Two runs that agree scientifically produce byte
identical snapshots.

usage: snapshot.py <step2_output_dir> <snapshot_dir> [--elide-environment]

--elide-environment additionally drops the host and dependency-version report
from the summary HTML.  Leave it off when comparing two runs that are supposed
to share a stack -- there, a version that moved is a real finding.  Turn it on
when comparing runs on deliberately different stacks, e.g. the numpy 1 and
numpy 2 legs of CI: those differ in the version block by construction, and
without this the genuine question -- did any *data* move -- is buried.
"""
import os
import re
import sys
import hashlib

import pandas as pd

TS = re.compile(r'[-_]\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}')
# plotly stamps a fresh uuid4 on every graph div, so it differs between runs
# even when the figure is identical.  Not a determinism signal.
UUID = re.compile(r'[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}')
BOOTSTRAP = re.compile(r'Bootstrap: ([0-9.]+)%')

# vsnp3_step2.py writes everything from "System Information:" to </body> as
# environment metadata: OS, CPU, interpreter, and a version line per dependency.
# The section is the last thing in the document, so it can be cut as one span
# rather than by matching each line -- which would rot the moment a line is
# added to the report.
ENVIRONMENT_BLOCK = re.compile(
    r'<h2>System Information:</h2>.*?(?=</body>)', re.DOTALL)


def norm_name(rel):
    return TS.sub('', rel)


def dump_excel(src, dst):
    sheets = pd.read_excel(src, sheet_name=None, engine='openpyxl', header=None)
    with open(dst, 'w') as out:
        for name in sorted(sheets):
            df = sheets[name]
            print(f'### sheet: {name}  shape={df.shape}', file=out)
            df.to_csv(out, sep='\t', index=False, header=False, na_rep='<NA>')


def dump_text(src, dst, elide_environment=False):
    with open(src, 'r', errors='replace') as f:
        body = f.read()
    if elide_environment:
        body = ENVIRONMENT_BLOCK.sub('<h2>System Information:</h2>\n<elided>\n', body)
    # the HTML summary and validation log embed the run timestamp inline
    body = TS.sub('', body)
    body = re.sub(r'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}', '', body)
    body = re.sub(r'runtime: [0-9:.]+', 'runtime: <elided>', body)
    body = re.sub(r'Time: [0-9:.]+', 'Time: <elided>', body)
    body = re.sub(r'Total run time: [0-9.]+ \w+', 'Total run time: <elided>', body)
    body = UUID.sub('<uuid>', body)
    with open(dst, 'w') as out:
        out.write(body)
    # Bootstrap support is the one number in the HTML tree that is supposed to
    # be reproducible, so pull it out where a diff can actually see it.
    values = BOOTSTRAP.findall(body)
    if values:
        with open(dst + '.bootstrap', 'w') as out:
            for v in values:
                print(v, file=out)


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = {a for a in sys.argv[1:] if a.startswith('--')}
    unknown = flags - {'--elide-environment'}
    if len(args) != 2 or unknown:
        sys.exit(f'usage: snapshot.py <step2_output_dir> <snapshot_dir> '
                 f'[--elide-environment]' +
                 (f'\nunknown option(s): {" ".join(sorted(unknown))}' if unknown else ''))
    src_root, dst_root = args
    elide_environment = '--elide-environment' in flags
    os.makedirs(dst_root, exist_ok=True)
    manifest = []

    for dirpath, _dirnames, filenames in os.walk(src_root):
        for fn in sorted(filenames):
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, src_root)
            ext = os.path.splitext(fn)[1].lower()
            if ext == '.zip':
                continue                      # archive of the inputs, not an output
            out_rel = norm_name(rel).replace(os.sep, '__')
            try:
                if ext == '.xlsx':
                    dump_excel(full, os.path.join(dst_root, out_rel + '.tsv'))
                    out_rel += '.tsv'
                elif ext in ('.tre', '.fasta', '.txt', '.html', '.newick'):
                    dump_text(full, os.path.join(dst_root, out_rel), elide_environment)
                else:
                    continue
            except Exception as e:                      # noqa: BLE001 - snapshot tool
                with open(os.path.join(dst_root, out_rel + '.ERROR'), 'w') as out:
                    out.write(f'{type(e).__name__}: {e}\n')
                out_rel += '.ERROR'
            digest = hashlib.md5(
                open(os.path.join(dst_root, out_rel), 'rb').read()).hexdigest()
            manifest.append(f'{digest}  {out_rel}')

    with open(os.path.join(dst_root, 'MANIFEST'), 'w') as out:
        for line in sorted(manifest):
            print(line, file=out)
    print(f'{len(manifest)} artifacts -> {dst_root}')


if __name__ == '__main__':
    main()
