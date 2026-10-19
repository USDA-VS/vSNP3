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

__version__ = "3.37"


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

# The mixed-sample signature, used by the step 1 report and the step 2 summary.
# A clonal isolate makes few confident mixed calls; a mixture of two strains makes
# one at every site where the strains differ.  A "confident mixed call" is a
# single-base AC=1 record with QUAL >= QUAL_THRESHOLD and MQ >= MQ_THRESHOLD.  A
# sample is flagged when it has at least MIXED_MIN_AC1 of them and they outnumber
# its SNPs that pass the step 2 thresholds by more than MIXED_AC1_RATIO to 1.
# The 3.37 validation: known mixtures had 20-87 confident mixed calls per passing
# SNP (most had none passing); clonal samples at most 0.41 (11 validation samples,
# 30 simulated isolates).  The limits were first 10 and 1.0.  On a B. suis run
# (Oct 2026) that flagged 13 samples, 7 of them at 1.0-2.8 calls per passing SNP,
# which says too little to act on; the other 6 were at 5.5-13.6.  Only the clearest
# are flagged now: a ratio of 5 is 12 times the clonal maximum and a quarter of the
# lowest known mixture, and a minimum of 20 is that mixture's count at one passing
# SNP.  The minimum also keeps a handful of stray calls from being flagged.
MIXED_MIN_AC1 = 20
MIXED_AC1_RATIO = 5.0


def mixed_signature(passing_snps, confident_ac1):
    '''True when a sample's calls look like a mixture of strains (see MIXED_MIN_AC1).'''
    return confident_ac1 >= MIXED_MIN_AC1 and confident_ac1 > MIXED_AC1_RATIO * passing_snps

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


# ------------------------------------------------------ long-read (nanopore) calls
#
# Nanopore reads are called by bcftools mpileup/call (vsnp3_longread.py), whose
# QUAL cannot exceed about 228 -- its genotype likelihoods are capped at PL 255 --
# so it never meant the same as freebayes' QUAL and the thresholds above could not
# read it.  Until 3.37 every nanopore QUAL had 100 added to it (NANOPORE_QUAL_OFFSET,
# removed).  Now QUAL is replaced by the evidence QUAL: the Phred-scaled posterior
# probability that the site is reference-only, from the read counts bcftools used,
# the per-read error rate measured from the sample's own alignment, and the variant
# prior below.  It is uncapped, so it is on the freebayes scale and the thresholds
# keep their meaning.  The parameters are written into every nanopore VCF header.

# Prior probability that a site is a variant: the 1e-3 that freebayes (--theta)
# and bcftools (-P) also use.
LONGREAD_VARIANT_PRIOR = 0.001

# The per-read error rate is measured from the BAM (samtools stats: mismatches and
# indels per aligned base) and clamped to this range.  At the floor a clean read
# is worth 27 (Phred), close to what freebayes gives a Q30 Illumina read; at the
# ceiling 4.8.  The default is used only when samtools stats cannot be read.
LONGREAD_ERROR_RATE_MIN = 0.002
LONGREAD_ERROR_RATE_MAX = 0.25
LONGREAD_ERROR_RATE_DEFAULT = 0.05

# A call whose ALT reads all lie on one strand, although both strands are
# sequenced at the site, is limited to the N band.  Applied only with at least
# LONGREAD_STRAND_MIN_ALT ALT reads and LONGREAD_STRAND_MIN_PER_STRAND reads on
# each strand, so a thin site is judged by its QUAL alone.
LONGREAD_STRAND_MIN_ALT = 4
LONGREAD_STRAND_MIN_PER_STRAND = 2

# A call in or beside a reference homopolymer run of this many bases or more is
# limited to the N band: homopolymers are where nanopore reads mis-measure and
# where their alignment is least reliable.  Measured on the shipped references,
# runs of 6 or more touch 0.18% of M. bovis AF2122 and 0.41% of B. suis 1; runs
# of 5 or more, 1.1% and 1.8%.
LONGREAD_HOMOPOLYMER_MIN = 6

# bcftools' own QUAL must exceed this for a record to be kept at all: the gate
# (`QUAL>20`) the nanopore path has always applied.
LONGREAD_CALLER_MIN_QUAL = 20


# ------------------------------------------------------ step 1 usability verdicts

# The FASTQ and Reference "usability" verdicts shown as tiles in the step 1
# HTML report.  A sample is Poor if any 'poor' limit is crossed, Questionable
# if any 'questionable' limit is crossed, otherwise Acceptable, and the report
# names the measurement that decided it.  They are not written to the stats
# workbook (see HTML_ONLY_KEYS in vsnp3_file_setup.py).
#
# Loosened in 3.37.  The previous limits called a sample Questionable below
# 70x depth or above 10% duplicate reads, and Poor below 99% of the reference
# covered (zero-coverage > 1%), which flagged routine, usable runs.
FASTQ_USABILITY = {
    'poor':         {'depth': 20, 'duplicate_ratio': 0.80, 'q20_percent': 50.0},
    'questionable': {'depth': 40, 'duplicate_ratio': 0.30, 'q30_percent': 70.0},
}
# genome_coverage is the fraction of reference positions with any coverage;
# zero_coverage_percent is the same measurement from the other side, kept so
# either limit can be tightened on its own.
REFERENCE_USABILITY = {
    'poor':         {'genome_coverage': 0.90, 'zero_coverage_percent': 10.0, 'good_snp_percent': 0.09},
    'questionable': {'genome_coverage': 0.97, 'zero_coverage_percent': 3.0, 'good_snp_percent': 0.03},
}


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
