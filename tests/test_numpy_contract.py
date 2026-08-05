#!/usr/bin/env python3
"""
The numeric contract: what must not change when numpy or pandas moves.

    python3 tests/test_numpy_contract.py

vSNP3 supports a range of numpy and pandas rather than a single version, so
something has to hold that range honest.  This file is that something.  It is
hermetic -- no reference set, no bioinformatics tools -- so it runs on every leg
of the CI matrix, and every check here is one that a version bump could break
silently.

Background, because the reasoning is not obvious from the assertions:

  * The bootstrap in vsnp3_html_tree.py uses np.random.default_rng.  numpy's
    Generator, unlike the legacy RandomState, carries NO stream-compatibility
    guarantee across versions.  If the stream ever shifts, every published
    bootstrap support value shifts with it and nothing downstream would say so.
    Section 1 freezes the stream.

  * vsnp3_assembly.py's N50 indexes a numpy array.  Getting that wrong is not
    hypothetical: `int()` on a 1-D array was a DeprecationWarning through numpy
    2.2 and became a TypeError in 2.4, and the assembly path is not covered by
    any end-to-end test.  Section 2 covers it.

  * Section 4 checks that the interpreter actually running is inside the window
    meta.yaml declares, so a CI leg cannot silently resolve outside it.
"""

import os
import sys
import hashlib
import re
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, 'bin')
sys.path.insert(0, BIN)

PASS, FAIL = [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f'\n        {detail}' if detail and not ok else ''))
    return ok


# ---------------------------------------------------------------------------
# [1] the seeded bootstrap stream
# ---------------------------------------------------------------------------

# sha256 over 100 replicate draws of 2000 string positions at BOOTSTRAP_SEED.
# Measured identical on numpy 1.26.0, 2.2.6, 2.4.3 and 2.4.6.  choice() with
# replace=True and no p= delegates to integers(), which is PCG64 plus Lemire
# rejection sampling -- integer arithmetic, so it is also platform independent.
#
# If this fires, do NOT re-baseline it.  It means a numpy upgrade changed the
# bootstrap support values in every tree the pipeline has ever published, and
# that is a decision to make deliberately, not a constant to edit.
BOOTSTRAP_FINGERPRINT = '48abac2f1cf3608c38ced5a8a3e63f3695a4f60ed8773a14479734e79f6bda30'


def section_rng():
    print('\n[1] the seeded bootstrap draws the same positions on every numpy')
    import numpy as np
    from vsnp3_html_tree import BOOTSTRAP_SEED

    check('BOOTSTRAP_SEED is still 456123', BOOTSTRAP_SEED == 456123,
          f'seed is now {BOOTSTRAP_SEED}; the fingerprint below was taken at 456123')

    # Mirrors calculate_bootstrap_support() in vsnp3_html_tree.py: one
    # default_rng seeded once, then n_replicates draws of len(positions).
    positions = ['NC_002945.4:%d' % p for p in range(1, 2001)]
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    blob = b''.join(rng.choice(positions, size=len(positions), replace=True).tobytes()
                    for _ in range(100))
    got = hashlib.sha256(blob).hexdigest()
    check(f'bootstrap stream unchanged under numpy {np.__version__}',
          got == BOOTSTRAP_FINGERPRINT,
          f'expected {BOOTSTRAP_FINGERPRINT}\n        got      {got}\n'
          f'        A numpy upgrade has changed every published support value. '
          f'Do not edit the constant -- investigate.')


# ---------------------------------------------------------------------------
# [2] N50 / L50
# ---------------------------------------------------------------------------

def _n50(lengths):
    """Run the real Assemble.stats() over a temp FASTA and return (n50, l50)."""
    import contextlib
    import io
    from vsnp3_assembly import Assemble
    with tempfile.TemporaryDirectory() as d:
        fasta = os.path.join(d, 'contigs.fasta')
        with open(fasta, 'w') as f:
            for i, n in enumerate(lengths):
                # SPAdes-style header so the coverage parse takes its normal path.
                f.write(f'>NODE_{i+1}_length_{n}_cov_10.0\n{"A" * n}\n')
        # Skip __init__: it wants FASTQs and a Setup working directory, and
        # stats() reaches for neither.
        obj = Assemble.__new__(Assemble)
        with contextlib.redirect_stdout(io.StringIO()):
            obj.stats(FASTA=fasta)
        return obj.n50, obj.l50


def section_n50():
    """
    N50 here is "the first contig, largest first, at which the cumulative length
    reaches half the assembly", and L50 is its 1-based rank.  That is the
    definition the code implements; these cases pin it down so a rewrite of the
    indexing cannot quietly change it.
    """
    print('\n[2] N50 indexes the numpy array correctly')

    # Distinct lengths, cumsum strictly increasing: one match.  csum is
    # [5000, 8000, 9000, 9500] against half=4750, so the first contig already
    # crosses it.
    n50, l50 = _n50([5000, 3000, 1000, 500])
    check('N50 over distinct contig lengths', (n50, l50) == (5000, 1),
          f'got n50={n50} l50={l50}, expected 5000 / 1')

    # Equal lengths: csum [1000, 2000, 3000, 4000] against half=2000, so the
    # crossing is the second contig.  An off-by-one in the index shows up here.
    n50, l50 = _n50([1000, 1000, 1000, 1000])
    check('N50 over equal contig lengths', (n50, l50) == (1000, 2),
          f'got n50={n50} l50={l50}, expected 1000 / 2')

    # Single contig: the whole assembly is its own N50.
    n50, l50 = _n50([4242])
    check('N50 over a single contig', (n50, l50) == (4242, 1),
          f'got n50={n50} l50={l50}, expected 4242 / 1')

    # The regression case.  A zero-length record makes the last two cumsum
    # entries equal, so np.where returns TWO indices and the matched value is
    # that duplicate.  int() on a size-2 array is a TypeError on every numpy
    # ever released -- this is what the previous `int(ind[0])` did, and it is
    # why the scalar has to be indexed out first.
    n50, l50 = _n50([1000, 0])
    check('N50 survives a zero-length contig (two-index np.where)',
          (n50, l50) == (1000, 1), f'got n50={n50} l50={l50}, expected 1000 / 1')


# ---------------------------------------------------------------------------
# [3] the zero-coverage array path
# ---------------------------------------------------------------------------

def section_zero_coverage():
    """
    vsnp3_zero_coverage.py replaced four whole-genome dicts with int32 arrays
    (1,482 MB -> 17 MB).  The arithmetic that made that possible is what is
    checked here: the int64 sum guard against int32 overflow, and flatnonzero
    returning 1-based positions after the +1.
    """
    print('\n[3] the coverage arrays sum and index as intended')
    import numpy as np

    # sum(dtype=np.int64) on an int32 array.  On LP64 Linux and macOS numpy
    # already promotes an int32 sum to the platform int, so the explicit dtype
    # is belt-and-braces there; it is load-bearing where the default int is
    # 32-bit.  Either way the total must be exact above 2**31, which a
    # high-coverage bacterial genome reaches.
    array = np.full(1000, 3_000_000, dtype=np.int32)
    guarded = int(array.sum(dtype=np.int64))
    check('int64 accumulator sums past 2**31 exactly', guarded == 3_000_000_000,
          f'got {guarded}')

    # flatnonzero over ==0, then +1 for 1-based coordinates.
    array = np.array([5, 0, 0, 7, 0], dtype=np.int32)
    zeros = [int(offset) + 1 for offset in np.flatnonzero(array == 0)]
    check('uncovered positions are 1-based', zeros == [2, 3, 5], f'got {zeros}')

    # A contig with no reads at all must still be counted as uncovered, which is
    # why the arrays are sized from the FASTA rather than from samtools output.
    empty = np.zeros(4, dtype=np.int32)
    check('a contig with no reads is fully uncovered',
          len(np.flatnonzero(empty == 0)) == 4)


# ---------------------------------------------------------------------------
# [4] the running versions are inside the declared window
# ---------------------------------------------------------------------------

def _declared(package):
    """Pull a run-requirement spec out of conda_build/meta.yaml."""
    meta = open(os.path.join(REPO, 'conda_build', 'meta.yaml')).read()
    run = meta.split('requirements:', 1)[-1].split('test:', 1)[0]
    for line in run.splitlines():
        line = line.strip()
        if line.startswith('- ') and re.match(rf'-\s+{re.escape(package)}([\s>=<]|$)', line):
            return line[2:].strip()
    return None


def _satisfies(version, spec):
    """
    Evaluate a conda spec like 'numpy >=1.26.4,<3' against a version string.
    Conda pads the shorter side with zeros, so 1.26.4 > 1.26 -- which is the
    whole reason '<=1.26' was a trap.  Tuple comparison reproduces that.
    """
    def parts(v):
        return tuple(int(x) if x.isdigit() else 0
                     for x in re.split(r'[._-]', v.split('+')[0])[:4])

    def pad(a, b):
        n = max(len(a), len(b))
        return a + (0,) * (n - len(a)), b + (0,) * (n - len(b))

    for clause in spec.split(','):
        m = re.match(r'\s*(>=|<=|==|!=|>|<)?\s*([0-9][0-9._a-z-]*)\s*$', clause)
        if not m:
            return None  # a form this checker does not model; do not guess
        op, want = m.group(1) or '==', m.group(2)
        got, exp = pad(parts(version), parts(want))
        if not {'>=': got >= exp, '<=': got <= exp, '==': got == exp,
                '!=': got != exp, '>': got > exp, '<': got < exp}[op]:
            return False
    return True


def section_declared_window():
    print('\n[4] the interpreter and libraries are inside the declared window')
    import numpy as np
    import pandas as pd

    for package, version in (('python', '.'.join(str(x) for x in sys.version_info[:3])),
                             ('numpy', np.__version__),
                             ('pandas', pd.__version__)):
        line = _declared(package)
        if not check(f'meta.yaml declares {package}', line is not None,
                     'no run requirement found'):
            continue
        spec = line[len(package):].strip()
        ok = _satisfies(version, spec)
        if ok is None:
            check(f'{package} spec is in a form this test understands', False, line)
            continue
        check(f'{package} {version} satisfies "{line}"', ok,
              f'this environment is outside the window meta.yaml declares.\n'
              f'        In CI that means a matrix leg resolved wrong. On a\n'
              f'        development machine it usually means the environment\n'
              f'        predates a change to the recipe -- rebuild it:\n'
              f'            conda env create -f conda_build/vsnp3.yaml --force')


def main():
    print('=' * 78)
    print('vSNP3 numeric contract')
    print('=' * 78)
    try:
        import numpy, pandas  # noqa: F401
    except ImportError as e:
        # numpy and pandas are hard requirements, not optional imports -- see
        # OPTIONAL_IMPORTS in test_smoke.py.  Absent means a broken environment.
        print(f'\n### numpy/pandas not importable: {e}')
        return 1
    section_rng()
    section_n50()
    section_zero_coverage()
    section_declared_window()
    print('\n' + '=' * 78)
    print(f'PASS {len(PASS)}   FAIL {len(FAIL)}')
    if FAIL:
        for name in FAIL:
            print(f'  - {name}')
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
