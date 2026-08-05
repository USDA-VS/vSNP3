#!/usr/bin/env python3

'''
The version and the SNP calling thresholds: the one place either is written.

Kept deliberately free of imports.  Every script in bin/ imports from this module,
so anything heavy added here becomes an import-time cost for the whole suite --
which is also why the version does not live in vsnp3_file_setup.py.

Version.  Before this module the literal was duplicated in 29 files and had
drifted: three said 3.36 while twenty-six said 3.35, so a run using the corrected
annotation code was stamped 3.35 in both the step1 run log and the step2 HTML
report.  Run `vsnp3_version.py --check` to confirm a whole install agrees.

Thresholds.  These were written out in vsnp3_zero_coverage.py,
vsnp3_group_reporter.py, vsnp3_group_on_defining_snps.py, vsnp3_step2.py and
vsnp3_vcf_merge_to_fasta.py, and had already drifted -- the SNP QUAL cutoff was 150
in three files and 300 in a fourth.  vsnp3_group_reporter.py also carried a comment
asserting the values "must be hardcoded", which is why step1's group call could not
honour step2's -w/-x/-y.

Numbers outside Python -- conda_build/meta.yaml, the READMEs, the internal/*.sh
wrappers -- cannot import this.  internal/bump_version.sh rewrites them together
with this file, and tests/test_release_consistency.py fails if any disagree.
'''

import os
import sys

__version__ = "3.36"


# --------------------------------------------------------------- SNP thresholds

# A position at or above this QUAL is called as the ALT it was reported as.
QUAL_THRESHOLD = 150

# Between N_THRESHOLD and QUAL_THRESHOLD the call is uncertain and becomes N.
# Below N_THRESHOLD there is not enough evidence to contradict the reference.
N_THRESHOLD = 50

# At least one position per group must reach this mapping quality.
MQ_THRESHOLD = 56

# freebayes allele count: 2 is homozygous ALT, 1 is a heterozygous / mixed site
# that gets an IUPAC ambiguity code.
AC_HOMOZYGOUS = 2
AC_HETEROZYGOUS = 1

# Deliberately stricter than QUAL_THRESHOLD, and deliberately NOT unified with it.
# vsnp3_vcf_merge_to_fasta.py emits a consensus FASTA: a wrong base there is baked
# permanently into a published sequence, with no MQ averaging and no cross-sample
# corroboration to catch it.  In the SNP table a position must clear thresholds in
# at least one sample and is visible in context, so a lower cutoff is defensible
# there and is not here.
CONSENSUS_QUAL_THRESHOLD = 300

# Reads averaging longer than this, with no R2, indicate long-read sequencing.
# Illumina 2x300 averages ~300 nt; ONT averages well over 1000, so the gap is wide.
# This was previously tested against the MAXIMUM read length, which meant a single
# long read in an Illumina library rerouted the whole sample to the nanopore caller.
NANOPORE_AVG_READ_LEN_CUTOFF = 701

# QUAL substituted for records whose FORMAT is GT only, when the caller passes
# --assume_gt_only_quality.  Such records (assembly- or consensus-derived VCFs)
# carry no quality data at all, so without a substitute every threshold rejects
# them.  Written into the VCF header as ##vsnp3_synthesized_quality so that a
# reader can tell a synthesized value from a called one.
SYNTHESIZED_QUAL = 999

# Added to every nanopore QUAL by vsnp3_alignment_vcf.py.  Recorded here because it
# makes nanopore QUAL values incomparable to the thresholds above, and to raw caller
# output, unless the reader knows the offset was applied.
NANOPORE_QUAL_OFFSET = 100


# ------------------------------------------------------------------ install check

def installed_versions(directory=None):
    '''
    {script name: version string or None} for every vsnp3_*.py beside this file.

    Reads the literal with a regex rather than importing, so a script whose
    third-party dependencies are missing still reports its version instead of
    raising -- the point is to check the install, and a half-installed environment
    is exactly when that matters.
    '''
    import re
    directory = directory or os.path.dirname(os.path.realpath(__file__))
    pattern = re.compile(r'^__version__\s*=\s*["\']([^"\']+)["\']', re.M)
    found = {}
    for name in sorted(os.listdir(directory)):
        if not name.startswith('vsnp3_') or not name.endswith('.py'):
            continue
        try:
            with open(os.path.join(directory, name), errors='replace') as f:
                body = f.read()
        except OSError:
            found[name] = None
            continue
        literal = pattern.search(body)
        if literal:
            found[name] = literal.group(1)
        elif 'from vsnp3_version import' in body or 'import vsnp3_version' in body:
            found[name] = __version__       # imports the single source of truth
        else:
            found[name] = None
    return found


def _check(directory=None):
    versions = installed_versions(directory)
    if not versions:
        print('No vsnp3_*.py scripts found beside this file.')
        return 1
    width = max(len(n) for n in versions)
    disagree, unknown = {}, []
    for name, version in versions.items():
        print(f'  {name:<{width}}  {version or "unknown"}')
        if version is None:
            unknown.append(name)
        elif version != __version__:
            disagree[name] = version
    print()
    where = directory or os.path.dirname(os.path.realpath(__file__))
    print(f'{len(versions)} scripts in {where}')
    if disagree:
        print(f'MISMATCH: expected {__version__}, but {len(disagree)} disagree:')
        for name, version in sorted(disagree.items()):
            print(f'  {name}: {version}')
        return 1
    if unknown:
        print(f'{len(unknown)} script(s) report no version: {", ".join(unknown)}')
        return 1
    print(f'All {len(versions)} scripts report {__version__}.')
    return 0


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Report the vSNP3 version, or check that a whole install agrees.')
    parser.add_argument('--check', action='store_true',
                        help='List every vsnp3_*.py beside this file with its version and '
                             'exit non-zero if any disagree.  Works against a repo checkout '
                             'or an installed $PREFIX/bin.')
    parser.add_argument('--directory', default=None,
                        help='Directory to check instead of the one holding this file.')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    if args.check:
        sys.exit(_check(args.directory))
    print(f'{os.path.basename(__file__)}: version {__version__}')
