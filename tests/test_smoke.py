#!/usr/bin/env python3
"""
Cheap checks that would have caught several defects at the moment they were typed.

    python3 tests/test_smoke.py

No external data and no bioinformatics tools required, so this can run on every
push.  It is deliberately not a unit test of the science -- see test_annotation.py
for that.
"""

import ast
import os
import re
import string
import subprocess
import sys
import py_compile
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, 'bin')

# Modules with no CLI: importable libraries and the version/constant carriers.
NO_CLI = {'vsnp3_run.py',
          'vsnp3_file_setup.py', 'vsnp3_input_validator.py', 'krona_lca_all.py',
          'defining_snps_print.py'}

# Third-party imports that may legitimately be absent in a minimal environment.
OPTIONAL_IMPORTS = {'seaborn', 'matplotlib', 'plotly', 'dask', 'distro',
                    'regex', 'Bio', 'openpyxl',
                    'xlsxwriter', 'humanize'}

PASS, FAIL = [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f'\n        {detail}' if detail and not ok else ''))
    return ok


def scripts():
    return sorted(f for f in os.listdir(BIN) if f.endswith('.py'))


def section_compile():
    """
    Compile every module and treat a SyntaxWarning as a failure.  Three invalid
    escape sequences were sitting in the tree; those become SyntaxError in a
    future Python, and there was nothing to notice them.
    """
    print('\n[1] every module compiles, and with no SyntaxWarning')
    import warnings
    for name in scripts():
        path = os.path.join(BIN, name)
        with open(path) as f:
            source = f.read()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always', SyntaxWarning)
            try:
                compile(source, path, 'exec')
            except SyntaxError as e:
                check(f'{name} compiles cleanly', False, str(e))
                continue
        problems = [f'line {w.lineno}: {w.message}' for w in caught
                    if issubclass(w.category, SyntaxWarning)]
        check(f'{name} compiles cleanly', not problems, '; '.join(problems))


def section_format_strings():
    """
    Every literal format string must actually parse.

    vsnp3_excel_merge_defining_snps.py shipped '{}_{_{}.xlsx' in both of its output
    filenames, which raises ValueError: unexpected '{' in field name.  The tool did
    all its work and then died at write time, so Merge_Defining_SNPs had never
    completed a single run.  This check costs under a second.
    """
    print('\n[2] literal format strings and f-strings parse')
    bad = []
    for name in scripts():
        path = os.path.join(BIN, name)
        tree = ast.parse(open(path).read(), path)
        for node in ast.walk(tree):
            # "literal".format(...)
            if (isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Attribute)
                    and node.func.attr == 'format'
                    and isinstance(node.func.value, ast.Constant)
                    and isinstance(node.func.value.value, str)):
                try:
                    list(string.Formatter().parse(node.func.value.value))
                except ValueError as e:
                    bad.append(f'{name}:{node.lineno}: {e}: {node.func.value.value!r}')
    check('no malformed format strings', not bad, '\n        '.join(bad))


def section_cli():
    print('\n[3] --help and --version work, and every version agrees')
    versions = {}
    for name in scripts():
        if name in NO_CLI:
            continue
        path = os.path.join(BIN, name)
        for flag in ('--help', '--version'):
            proc = subprocess.run([sys.executable, path, flag],
                                  capture_output=True, text=True, timeout=120)
            out = (proc.stdout + proc.stderr).strip()
            ok = proc.returncode == 0
            if not ok and any(m in out for m in ('ModuleNotFoundError', 'ImportError')):
                missing = out.rsplit("'", 2)[-2] if "'" in out else ''
                if missing in OPTIONAL_IMPORTS:
                    print(f'  SKIP  {name} {flag}   (optional dependency {missing} absent)')
                    continue
            check(f'{name} {flag}', ok, out[-300:])
            if flag == '--version' and ok:
                token = out.split()
                if token:
                    versions[name] = token[-1]
                else:
                    check(f'{name} --version prints something', False,
                          'exited 0 but produced no output; if this is a library '
                          'module, add it to NO_CLI')
    distinct = sorted(set(versions.values()))
    check(f'all {len(versions)} scripts report one version', len(distinct) <= 1,
          f'found {distinct}')


def section_packaging():
    """
    conda_build/build.sh installs with `mv bin/*.py $PREFIX/bin`, so anything
    matching that glob is shipped to every user.  Tests must not be.
    """
    print('\n[4] packaging')
    stray = [f for f in scripts() if f.startswith('test_')]
    check('no test files in bin/', not stray, f'would be installed: {stray}')

    version_module = os.path.join(BIN, 'vsnp3_version.py')
    check('bin/vsnp3_version.py exists', os.path.exists(version_module))

    # Every script must import the version rather than carrying a literal.
    literals = []
    for name in scripts():
        if name == 'vsnp3_version.py':
            continue
        body = open(os.path.join(BIN, name)).read()
        if '__version__ = "' in body or "__version__ = '" in body:
            literals.append(name)
    check('no duplicated version literals', not literals, f'{literals}')


def section_sibling_imports():
    """
    bin/ is not a package; scripts resolve each other because sys.path[0] is the
    directory holding the invoked script.  Confirm that still holds from elsewhere,
    which is what an installed $PREFIX/bin/vsnp3_step2.py relies on.
    """
    print('\n[5] sibling imports resolve when invoked by path from another cwd')
    with tempfile.TemporaryDirectory() as cwd:
        for entry in ('vsnp3_step1.py', 'vsnp3_step2.py'):
            proc = subprocess.run([sys.executable, os.path.join(BIN, entry), '--version'],
                                  capture_output=True, text=True, cwd=cwd, timeout=120)
            check(f'{entry} runs from an unrelated cwd', proc.returncode == 0,
                  (proc.stdout + proc.stderr)[-300:])


def section_report():
    """
    Render the step 1 report and convert it.

    The report used to be a LaTeX file nothing checked, which is how it came to
    be produced by a pdflatex call whose exit status was discarded.  Both halves
    are exercised here instead: the HTML always, the PDF whenever the renderer
    is installed.  jinja2 and weasyprint are optional at import time, so a
    minimal environment skips rather than fails.
    """
    print('\n[6] step 1 report')
    sys.path.insert(0, BIN)
    try:
        from vsnp3_file_setup import HtmlReport, redact_paths
    except ImportError as e:
        print(f'  SKIP  report render   ({e})')
        return

    # A report is sent to whoever needs the result, so it must not carry the
    # directory layout of the machine that produced it.  These are real log lines.
    call = ('SYSTEM CALL: bwa mem -M -R "@RG\\tID:B18-0389\\tSM:B18-0389" -t 8 '
            '/home/tstuber/analysis/script_test_files/step1/test/B18-0389/'
            'NC_006932-NC_006933.fasta '
            '/home/tstuber/analysis/script_test_files/step1/test/B18-0389/'
            'B18-0389_S3_L001_R1_cut.fastq.gz > B18-0389.sam -- 2026-08-04_15:10:14')
    scrubbed = redact_paths(call)
    check('absolute paths are replaced by the filename',
          '/home/tstuber' not in scrubbed
          and 'NC_006932-NC_006933.fasta' in scrubbed
          and 'B18-0389_S3_L001_R1_cut.fastq.gz' in scrubbed, scrubbed)
    check('the rest of the command line survives redaction',
          'bwa mem -M -R' in scrubbed and '@RG\\tID:B18-0389' in scrubbed
          and '> B18-0389.sam' in scrubbed and '2026-08-04_15:10:14' in scrubbed,
          scrubbed)
    check('/dev/null is left alone',
          redact_paths('samtools fastq 2> /dev/null') .endswith('2> /dev/null'))
    check('a URL is not mistaken for a path',
          redact_paths('htslib http://www.htslib.org/ ok')
          == 'htslib http://www.htslib.org/ ok',
          redact_paths('htslib http://www.htslib.org/ ok'))
    check('a relative path is left alone',
          redact_paths('spades -o spades_assembly/scaffolds.fasta')
          == 'spades -o spades_assembly/scaffolds.fasta')
    stats = {'sample': 'TEST-001', 'date': '2026-01-01',
             'Reference': 'NC_002945.4', 'Reference Length': '4,349,904',
             'Average Depth': '88.2', 'Genome with Coverage': '99.4%',
             'Quality SNPs': '412', 'FASTQ Usability': 'Acceptable',
             'FASTQ_R1': 'TEST-001_R1.fastq.gz', 'FASTQ_R2': 'TEST-001_R2.fastq.gz',
             'R1 Read Count': '1,000,000', 'R2 Read Count': '1,000,000',
             'Spoligotype SB Number': 'SB0140', 'Groups': 'Mbovis-01'}
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmp:
        try:
            os.chdir(tmp)
            report = HtmlReport('TEST-001')
            report.add_notice('a notice the reader needs to see')
            report.add_table('Sourmash reference similarity',
                             ('Similarity', 'ID'), [('99.8%', 'NC_002945.4')])
            report.program_versions = f'bwa 0.7.17\n{call}'
            try:
                html_path = report.write(stats)
            except ImportError as e:
                print(f'  SKIP  report render   ({e})')
                return
            html = open(html_path, encoding='utf-8').read()
            check('HTML report contains the sample and a tile',
                  'TEST-001' in html and 'Average depth' in html)
            check('HTML report contains the notice and the extra table',
                  'a notice the reader needs to see' in html
                  and 'NC_002945.4' in html)
            # The report is read next to its own PDF; a dark-scheme override made
            # one of the two near-black and the other white.
            check('HTML report does not follow a dark colour scheme',
                  not re.search(r'@media[^{]*prefers-color-scheme', html))
            check('HTML report writes no LaTeX beside itself',
                  not [f for f in os.listdir('.') if f.endswith(('.tex', '.aux'))],
                  f'found {sorted(os.listdir("."))}')
            check('rendered report discloses no absolute path',
                  '/home/tstuber' not in html and '/data/' not in html,
                  'an absolute path reached the rendered report')
            try:
                import weasyprint                              # noqa: F401
            except ImportError as e:
                print(f'  SKIP  report PDF   ({e})')
                return
            pdf_path = report.to_pdf()
            ok = check('PDF renders from the HTML report',
                       bool(pdf_path) and os.path.exists(pdf_path or ''))
            if ok:
                with open(pdf_path, 'rb') as f:
                    head = f.read(5)
                check('PDF is a real PDF and not an empty file',
                      head == b'%PDF-' and os.path.getsize(pdf_path) > 1000,
                      f'starts {head!r}, {os.path.getsize(pdf_path)} bytes')
        finally:
            os.chdir(cwd)


def section_grouping():
    """
    A group in which every position carries the same call in every sample.

    There is nothing to align, so no tree can be built.  make_groupings used to
    represent each sample's empty result as a bare DataFrame with no columns, and
    dict_to_fasta then raised KeyError: 'abs_pos' -- which reads as one malformed
    VCF, names whichever sample happened to be iterated first, and propagated out
    of Group.__init__ so every other group's tree was lost with it.
    """
    print('\n[7] step 2 grouping: a group with no informative positions')
    sys.path.insert(0, BIN)
    try:
        import pandas as pd
        from vsnp3_group_on_defining_snps import Group
    except ImportError as e:
        print(f'  SKIP  grouping   ({e})')
        return

    def frame(alts):
        positions = (100, 200, 300)
        return pd.DataFrame({
            'abs_pos': [f'NC_002945v4:{p}' for p in positions],
            'CHROM': ['NC_002945v4'] * 3, 'POS': list(positions),
            'REF': ['A', 'C', 'G'], 'ALT': alts,
            'QUAL': [900.0] * 3, 'MQ': [60.0] * 3, 'AC': [2] * 3,
        })

    group = Group.__new__(Group)          # the methods, without the pipeline
    group.debug = False
    group.n_threshold, group.qual_threshold = 50, 150
    group.ambigious_lookup = {'AG': 'R'}
    informative = {'sampleA': frame(['T', 'T', 'A']),
                   'sampleB': frame(['T', 'G', 'A']),
                   'sampleC': frame(['C', 'T', 'A'])}
    uninformative = {'04-1161-fixed': frame(['T', 'T', 'A']),
                     '08-1137_TAM_TX_Fed': frame(['T', 'T', 'A'])}
    group.dataframes_names_updated = {**informative, **uninformative}

    finished = {}
    for item in (('Mbovis-01', informative), ('Mbovis-02', uninformative)):
        parsimony_sample_dict, name = group.make_groupings(item)
        finished[name] = parsimony_sample_dict
    # The filter Group.__init__ applies to the same dict.
    finished = {i: j for i, j in finished.items() if j != {}}

    check('a group with no informative positions is dropped, not fatal',
          'Mbovis-02' not in finished, f'survived as {finished.get("Mbovis-02")}')
    check('the informative group is unaffected', 'Mbovis-01' in finished)

    with tempfile.TemporaryDirectory() as tmp:
        fasta = os.path.join(tmp, 'Mbovis-01.fasta')
        try:
            group.dict_to_fasta(finished['Mbovis-01'], fasta)
            records = open(fasta).read().split('>')[1:]
            seqs = dict(r.split('\n')[:2] for r in records)
            ok = True
        except Exception as e:                                   # noqa: BLE001
            seqs, ok = {}, False
            check('the informative group still writes an alignment', False, repr(e))
        if ok:
            check('the informative group still writes an alignment',
                  seqs == {'sampleA': 'TT', 'sampleB': 'TG',
                           'sampleC': 'CT', 'root': 'AC'}, str(seqs))
            check('every sequence is the same length, as RAxML requires',
                  len({len(s) for s in seqs.values()}) == 1, str(seqs))


def section_sample_names():
    """
    Sample-name resolution and duplicate handling.

    resolve_sample_name was lifted out of a nested try/except chain so that
    collisions can be settled before any sample is processed; the substitutions
    are cumulative and the first pattern's '.' is a regex any-character, both of
    which are easy to change by accident.  Checked against the original chain.
    """
    print('\n[8] step 2 sample names')
    sys.path.insert(0, BIN)
    try:
        import pandas as pd
        from vsnp3_group_on_defining_snps import Group
    except ImportError as e:
        print(f'  SKIP  sample names   ({e})')
        return

    meta = pd.DataFrame({'file_name': ['a', 'b.vcf', 'c', 'd', '99-0100', 'e_zc'],
                         'metadata':  ['A', 'B', 'C', 'D', 'ZZ', 'E']})
    cases = {
        'a': 'A',                          # direct hit
        'c_zc.vcf': 'C',                   # .vcf then _zc stripped, cumulatively
        'd_zc_Val_TS.vcf': 'D',            # _zc_* stripped
        '99-0100_zc_Val_TS.vcf': 'ZZ',
        'e_zc.vcf': 'E',
        'unknown_zc.vcf': 'unknown',       # no hit: fully stripped basename
    }
    ok = True
    for name, want in cases.items():
        got = Group.resolve_sample_name(name, meta, True)
        if got != want:
            ok = False
            check(f'resolve_sample_name({name!r})', False, f'want {want!r}, got {got!r}')
    if ok:
        check(f'resolve_sample_name maps {len(cases)} name shapes as before', True)
    check('no metadata leaves the basename untouched',
          Group.resolve_sample_name('x_zc.vcf', meta, False) == 'x_zc.vcf')

    # One file matching several worksheet rows keeps one name, so no count changes
    # and nothing else in a run would mention it.  The row count is what lets the
    # caller report it at all.
    ambiguous = pd.DataFrame({'file_name': ['dup', 'dup', 'solo'],
                              'metadata':  ['first', 'second', 'only']})
    name, key, hits = Group.resolve_sample_name_detail('dup.vcf', ambiguous, True)
    check('a file matching two worksheet rows takes the first and reports the count',
          (name, key, hits) == ('first', 'dup', 2), f'got {(name, key, hits)}')
    name, key, hits = Group.resolve_sample_name_detail('solo.vcf', ambiguous, True)
    check('an unambiguous match reports one row',
          (name, key, hits) == ('only', 'solo', 1), f'got {(name, key, hits)}')
    name, key, hits = Group.resolve_sample_name_detail('nope.vcf', ambiguous, True)
    check('an unmatched file reports no key and no rows',
          (name, key, hits) == ('nope', None, 0), f'got {(name, key, hits)}')


def main():
    print('=' * 78)
    print('vSNP3 smoke tests')
    print('=' * 78)
    section_compile()
    section_format_strings()
    section_cli()
    section_packaging()
    section_sibling_imports()
    section_report()
    section_grouping()
    section_sample_names()
    print('\n' + '=' * 78)
    print(f'PASS {len(PASS)}   FAIL {len(FAIL)}')
    if FAIL:
        for name in FAIL:
            print(f'  - {name}')
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
