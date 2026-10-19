#!/usr/bin/env python3

'''
Long-read (Oxford Nanopore) SNP calling for vSNP3 step 1.

vSNP3 calls Illumina reads with freebayes and reads every VCF, in step 1 and
step 2, through the same thresholds: a single-base AC=2 record is a SNP at
QUAL > QUAL_THRESHOLD (150) and MQ >= MQ_THRESHOLD (56); between N_THRESHOLD
(50) and 150 it is N; below 50 the sample is taken to match the reference.
freebayes' QUAL is the Phred-scaled posterior probability that the site is a
variant, and it grows with the evidence: about 30 per clean Illumina read, into
the thousands at normal depth.

Nanopore reads were called with bcftools mpileup/call, whose QUAL is the same
quantity but cannot exceed about 228, because the genotype likelihoods behind
it are capped at PL 255.  So a strong nanopore call sat at 225-228, a moderate
one anywhere below, and nothing could reach the 300 consensus threshold of
vsnp3_vcf_merge_to_fasta.py.  The remedy until now was to add 100 to every
nanopore QUAL, which moved the bands but did not make them mean the same thing.

This module makes a nanopore VCF that the thresholds read correctly:

  1. bcftools mpileup is run with its nanopore profile (-X ont: no BAQ, bases at
     Q5 or better, base quality capped at 30, no indel calling) and bcftools
     call -m with ploidy 2, so AC=1 and AC=2 mean what they mean for freebayes
     (mixed site, homozygous ALT).  bcftools decides which sites are variants;
     its own QUAL is kept in INFO/VSNP_RAWQ and its genotype, where the model
     below decides otherwise, in INFO/VSNP_CALLER_GT.

  2. QUAL is replaced by the evidence QUAL: the Phred-scaled posterior
     probability that the site is reference-only, computed from the
     high-quality read counts bcftools used (DP4), a per-read error rate
     measured from the sample's own alignment (samtools stats: mismatches and
     indels per aligned base), and the standard variant prior of 0.001.  It is
     not capped, so it is on the freebayes scale: a clean site scores about
     10*log10((1-e)/e) per supporting read -- 18 at e = 1.5%, 27 at e = 0.2% --
     less the prior, and thresholds keep the meaning they have for Illumina.
     A nanopore read carries less evidence than an Illumina read, so more
     nanopore reads are needed for the same confidence, which is as it should be.
     QUAL is never below N_THRESHOLD: bcftools reported the site as a variant
     (its own QUAL above LONGREAD_CALLER_MIN_QUAL), so however thin the
     evidence, the sample is uncertain there, and step 2 shows N rather than
     the reference base.  Without the floor a 10x nanopore sample showed the
     reference base at true SNPs supported by four or five reads.

     The genotype comes from the same model: the site is mixed (AC=1, GT 0/1)
     when the mixed hypothesis is more likely than homozygous ALT under the
     measured error rate, and homozygous (AC=2, GT 1/1) otherwise.  bcftools
     trusts each base's own quality, which is capped at 30, so four
     reference-looking reads among thirty made it call a clonal sample mixed
     at a true SNP (a systematic error, the usual nanopore kind); at a 2%
     error rate the homozygous explanation is far more likely.  A mixed call
     needs a minor allele of roughly a fifth of the reads at 30x; a 60/40
     mixture is still mixed at every site where the strains differ.

  3. Two nanopore error modes that read counts cannot see are placed in the N
     band (QUAL between N_THRESHOLD and QUAL_THRESHOLD - 1, so step 2 reads the
     sample as N there; FILTER and INFO/VSNP_FLAG say why): ALT support on one
     strand only when both strands are sequenced, and a site in or beside a
     reference homopolymer run of LONGREAD_HOMOPOLYMER_MIN bases or more, where
     alignment of the indel-prone reads is least reliable.  N is the safe
     outcome: the sample is uncertain there, and never a false SNP or a false
     reference call.

Every parameter is recorded in the VCF header (##vsnp3_platform,
##vsnp3_longread_caller, ##vsnp3_longread_model), so a VCF says how it was
made and step 2 can name the long-read samples in its summary.

    vsnp3_longread.py -r reference.fasta -b sample_nodup.bam
    vsnp3_longread.py -r reference.fasta -b sample_nodup.bam --vcf raw_bcftools.vcf -o out.vcf

The first form calls and models; the second applies the model to a VCF
bcftools call has already written.
'''

from vsnp3_version import __version__
from vsnp3_version import (QUAL_THRESHOLD, N_THRESHOLD, LONGREAD_VARIANT_PRIOR,
                             LONGREAD_ERROR_RATE_MIN, LONGREAD_ERROR_RATE_MAX,
                             LONGREAD_ERROR_RATE_DEFAULT, LONGREAD_STRAND_MIN_ALT,
                             LONGREAD_STRAND_MIN_PER_STRAND, LONGREAD_HOMOPOLYMER_MIN,
                             LONGREAD_CALLER_MIN_QUAL)

import os
import re
import sys
import math
import argparse
import textwrap
import subprocess

import vsnp3_run

# What `bcftools mpileup -X ont` expands to (bcftools 1.17 and later).  Used
# verbatim on a bcftools that predates the profiles, so the same calling runs
# everywhere.  No BAQ, bases at Q5 or better, base quality capped at 30, no
# indel calling.
ONT_PROFILE_FLAGS = ['-B', '-Q', '5', '--max-BQ', '30', '-I']
# bcftools' default per-file depth limit is 250, which silently subsamples a
# deeply sequenced site.  freebayes has no such limit.
MPILEUP_MAX_DEPTH = 10000

FLAG_STRAND = 'strand'
FLAG_HOMOPOLYMER = 'homopolymer'
FILTER_NAMES = {FLAG_STRAND: 'vsnp_strand', FLAG_HOMOPOLYMER: 'vsnp_homopolymer'}


# ------------------------------------------------------------------ the model

def per_read_weight(error_rate):
    '''Phred-scaled evidence one clean read adds for ALT over REF: 10*log10((1-e)/e).'''
    e = min(max(float(error_rate), 1e-6), 0.5)
    return 10.0 * math.log10((1.0 - e) / e)


def evidence_qual(ref_fwd, ref_rev, alt_fwd, alt_rev, error_rate, prior=LONGREAD_VARIANT_PRIOR):
    '''
    The evidence QUAL: -10*log10 P(the site is reference-only | the reads).

    r reads show the reference base and a the alternate base, each read wrong
    with probability e (a read showing the other base at a true homozygous site).
    Three hypotheses with the usual prior -- reference (1 - prior), mixed and
    homozygous ALT (prior/2 each):

        L(ref) = (1-e)^r * e^a      L(mixed) = 0.5^(r+a)      L(alt) = e^r * (1-e)^a

    For a clean site this is close to (a - r) * per_read_weight(e) less about 33
    for the prior, which is the shape of freebayes' QUAL on Illumina reads.  A
    site with no reads scores 0.  Computed in log space, so depth is unlimited.
    '''
    r = int(ref_fwd) + int(ref_rev)
    a = int(alt_fwd) + int(alt_rev)
    if r + a <= 0:
        return 0.0
    e = min(max(float(error_rate), 1e-6), 0.5)
    prior = min(max(float(prior), 1e-9), 0.5)
    log_e, log_1e, log_half = math.log10(e), math.log10(1.0 - e), math.log10(0.5)
    lp_ref = math.log10(1.0 - prior) + r * log_1e + a * log_e
    lp_het = math.log10(prior / 2.0) + (r + a) * log_half
    lp_alt = math.log10(prior / 2.0) + r * log_e + a * log_1e
    top = max(lp_ref, lp_het, lp_alt)
    log_total = top + math.log10(10.0 ** (lp_ref - top) + 10.0 ** (lp_het - top)
                                 + 10.0 ** (lp_alt - top))
    qual = -10.0 * (lp_ref - log_total)
    return round(max(0.0, qual), 1)


def mixed_by_evidence(ref_fwd, ref_rev, alt_fwd, alt_rev, error_rate):
    '''
    True when the mixed hypothesis (half the reads each) explains the counts
    better than homozygous ALT does, each read wrong with probability e.  The
    two share the variant prior, so this is a plain likelihood comparison.
    '''
    r = int(ref_fwd) + int(ref_rev)
    a = int(alt_fwd) + int(alt_rev)
    if r + a <= 0:
        return False
    e = min(max(float(error_rate), 1e-6), 0.5)
    ll_het = (r + a) * math.log10(0.5)
    ll_alt = r * math.log10(e) + a * math.log10(1.0 - e)
    return ll_het > ll_alt


def homopolymer_touching(seq, index0):
    '''
    Length of the longest run of one base that contains the 0-based position
    index0 or ends immediately beside it.  A call inside a run, or one base off
    its end (the base the run would absorb if a read mis-measured it), both count.
    '''
    n = len(seq)
    if n == 0 or index0 < 0 or index0 >= n:
        return 0
    longest = 0
    for i in (index0 - 1, index0, index0 + 1):
        if i < 0 or i >= n:
            continue
        base = seq[i]
        left = i
        while left > 0 and seq[left - 1] == base:
            left -= 1
        right = i
        while right + 1 < n and seq[right + 1] == base:
            right += 1
        longest = max(longest, right - left + 1)
    return longest


def one_strand_only(ref_fwd, ref_rev, alt_fwd, alt_rev,
                    min_alt=LONGREAD_STRAND_MIN_ALT,
                    min_per_strand=LONGREAD_STRAND_MIN_PER_STRAND):
    '''
    True when the ALT allele is seen on one strand only although both strands
    are sequenced at the site (at least min_per_strand reads on each) and there
    are enough ALT reads (min_alt) for that to be more than chance.
    '''
    ref_fwd, ref_rev, alt_fwd, alt_rev = (int(ref_fwd), int(ref_rev), int(alt_fwd), int(alt_rev))
    if alt_fwd + alt_rev < min_alt:
        return False
    if ref_fwd + alt_fwd < min_per_strand or ref_rev + alt_rev < min_per_strand:
        return False
    return alt_fwd == 0 or alt_rev == 0


# ------------------------------------------------------------ measurements

def bcftools_version():
    try:
        out = subprocess.run(['bcftools', '--version'], capture_output=True, text=True,
                             check=False).stdout
        first = out.strip().splitlines()[0] if out.strip() else ''
        return first.split()[1] if len(first.split()) > 1 else 'unknown'
    except (OSError, subprocess.SubprocessError):
        return 'unknown'


def mpileup_profile_args():
    '''
    (arguments, description) for the nanopore profile: `-X ont` when this
    bcftools has the profiles (1.17 and later), the same flags written out when
    it does not.  Either way the run log says which.
    '''
    try:
        out = subprocess.run(['bcftools', 'mpileup', '-X', 'list'], capture_output=True,
                             text=True, check=False)
        text = out.stdout + out.stderr
        if re.search(r'^\s*ont\s*$', text, re.M) or re.search(r'^ont\b', text, re.M):
            return ['-X', 'ont'], '-X ont (' + ' '.join(ONT_PROFILE_FLAGS) + ')'
    except (OSError, subprocess.SubprocessError):
        pass
    return list(ONT_PROFILE_FLAGS), ' '.join(ONT_PROFILE_FLAGS) + ' (this bcftools has no -X profiles)'


def measure_error_rate(bam):
    '''
    (error rate, description) for the alignment: samtools stats' "error rate",
    mismatches and indels per aligned base from the NM tags, clamped to
    [LONGREAD_ERROR_RATE_MIN, LONGREAD_ERROR_RATE_MAX].  Counting indels makes
    the rate larger than the substitution rate alone, so the evidence QUAL is
    conservative.  When samtools stats cannot be read the default is used, and
    the description says so.
    '''
    try:
        proc = vsnp3_run.run(['samtools', 'stats', bam])
        text = vsnp3_run._text(proc.stdout)
    except vsnp3_run.ExternalToolError as e:
        return LONGREAD_ERROR_RATE_DEFAULT, (f'default {LONGREAD_ERROR_RATE_DEFAULT} '
                                             f'(samtools stats failed: {str(e).splitlines()[0]})')
    m = re.search(r'^SN\terror rate:\t([0-9.eE+-]+)', text, re.M)
    if not m:
        return LONGREAD_ERROR_RATE_DEFAULT, (f'default {LONGREAD_ERROR_RATE_DEFAULT} '
                                             f'(no error rate line in samtools stats)')
    measured = float(m.group(1))
    rate = min(max(measured, LONGREAD_ERROR_RATE_MIN), LONGREAD_ERROR_RATE_MAX)
    detail = f'{measured:.4f} measured by samtools stats (mismatches and indels per aligned base)'
    if rate != measured:
        detail += f', clamped to {rate}'
    return rate, detail


# ------------------------------------------------------------------ calling

def call_variants(reference, bam, raw_vcf, threads=8, log=None):
    '''
    bcftools mpileup (nanopore profile) | bcftools call -mv, ploidy 2, to raw_vcf.
    Returns a description of the caller for the VCF header and the run log.
    '''
    profile, profile_text = mpileup_profile_args()
    mpileup = (['bcftools', 'mpileup', '--threads', str(threads)] + profile
               + ['-d', str(MPILEUP_MAX_DEPTH), '-a', 'FORMAT/AD,FORMAT/DP', '-Ou',
                  '-f', reference, bam])
    call = ['bcftools', 'call', '--threads', str(threads), '-mv', '--ploidy', '2',
            '-Ov', '-o', raw_vcf]
    vsnp3_run.run_pipe(mpileup, call, log=log)
    # bcftools call exits 0 with a header and no record for a sample identical
    # to the reference; that is an answer, not a failure (freebayes is treated
    # the same way in vsnp3_alignment_vcf.py).
    vsnp3_run.require_output(raw_vcf, vcf_header=True, what='bcftools call output')
    return (f'bcftools {bcftools_version()} mpileup {profile_text} -d {MPILEUP_MAX_DEPTH} '
            f'| call -mv --ploidy 2')


def read_fasta(path):
    '''{sequence id: upper-case sequence} for the reference the reads were aligned to.'''
    seqs, name, parts = {}, None, []
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(parts).upper()
    return seqs


def _info_value(info, key):
    m = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]*)', info)
    return m.group(1) if m else None


def _dp4(info, fmt, sample):
    '''
    (ref_fwd, ref_rev, alt_fwd, alt_rev), from INFO/DP4 (bcftools: the
    high-quality bases the caller used), else from FORMAT/AD with no strand
    information, else None.
    '''
    dp4 = _info_value(info, 'DP4')
    if dp4:
        try:
            values = [int(v) for v in dp4.split(',')]
            if len(values) == 4:
                return tuple(values)
        except ValueError:
            pass
    if fmt and sample:
        keys, values = fmt.split(':'), sample.split(':')
        if 'AD' in keys and len(values) == len(keys):
            try:
                ad = [int(v) for v in values[keys.index('AD')].split(',')]
                if len(ad) >= 2:
                    return ad[0], 0, sum(ad[1:]), 0
            except ValueError:
                pass
    return None


def model_header_lines(caller_text, error_rate, error_detail,
                       prior=LONGREAD_VARIANT_PRIOR, strand_min_alt=LONGREAD_STRAND_MIN_ALT,
                       strand_min_per_strand=LONGREAD_STRAND_MIN_PER_STRAND,
                       homopolymer_min=LONGREAD_HOMOPOLYMER_MIN, limit=None):
    limit = QUAL_THRESHOLD - 1 if limit is None else limit
    weight = per_read_weight(error_rate)
    return [
        '##vsnp3_platform=nanopore',
        f'##vsnp3_longread_caller={caller_text}; records with caller QUAL <= '
        f'{LONGREAD_CALLER_MIN_QUAL} or not single-base SNPs removed',
        f'##vsnp3_longread_model=error_rate={error_rate:.4f};error_rate_source={error_detail};'
        f'prior={prior};weight_per_read={weight:.1f};floor_qual={N_THRESHOLD};genotype=evidence model;'
        f'strand_min_alt={strand_min_alt};strand_min_per_strand={strand_min_per_strand};'
        f'homopolymer_min={homopolymer_min};limit_qual={limit}',
        '##INFO=<ID=VSNP_RAWQ,Number=1,Type=Float,Description="QUAL as reported by bcftools '
        'call, before the vSNP3 long-read evidence model replaced it">',
        '##INFO=<ID=VSNP_EVQ,Number=1,Type=Float,Description="vSNP3 evidence QUAL: Phred-scaled '
        'posterior probability that the site is reference-only, from DP4 read counts, the '
        'measured per-read error rate and the variant prior, before the floor at '
        f'{N_THRESHOLD} and the strand and homopolymer limits">',
        '##INFO=<ID=VSNP_CALLER_GT,Number=1,Type=String,Description="genotype as reported by '
        'bcftools call, where the vSNP3 long-read evidence model decided the genotype '
        'differently (AC and GT carry the model\'s decision)">',
        '##INFO=<ID=VSNP_FLAG,Number=.,Type=String,Description="strand: ALT seen on one strand '
        'only although both strands are sequenced; homopolymer: in or beside a reference '
        f'homopolymer run of at least {homopolymer_min} bases. A flagged record is given a QUAL '
        f'between {N_THRESHOLD} and {limit}, so step 2 reads it as N">',
        f'##FILTER=<ID={FILTER_NAMES[FLAG_STRAND]},Description="ALT support on one strand only">',
        f'##FILTER=<ID={FILTER_NAMES[FLAG_HOMOPOLYMER]},Description="in or beside a homopolymer '
        f'run of at least {homopolymer_min} bases">',
    ]


def apply_model(raw_vcf, out_vcf, reference_fasta, error_rate, error_detail, caller_text,
                prior=LONGREAD_VARIANT_PRIOR, strand_min_alt=LONGREAD_STRAND_MIN_ALT,
                strand_min_per_strand=LONGREAD_STRAND_MIN_PER_STRAND,
                homopolymer_min=LONGREAD_HOMOPOLYMER_MIN, caller_min_qual=LONGREAD_CALLER_MIN_QUAL):
    '''
    Rewrite the bcftools VCF raw_vcf as out_vcf with the evidence QUAL.

    Kept: single-base SNP records with caller QUAL above caller_min_qual.
    Each keeps its genotype, INFO and FORMAT, gains VSNP_RAWQ and VSNP_EVQ, and
    carries VSNP_FLAG and a FILTER name when a rule applied.  QUAL is the
    evidence QUAL, never below N_THRESHOLD (the caller reported a variant, so
    the sample is at least uncertain here); a flagged record's QUAL is placed in
    the N band whatever its evidence.  Returns the counts of what was done.
    '''
    seqs = read_fasta(reference_fasta)
    limit = QUAL_THRESHOLD - 1
    counts = {'records': 0, 'not_snp': 0, 'below_caller_qual': 0, 'no_depth': 0,
              'written': 0, 'ac1': 0, 'ac2': 0, 'flag_strand': 0, 'flag_homopolymer': 0,
              'limited': 0, 'floored': 0, 'genotype_changed': 0, 'contig_missing': 0}
    header_added = model_header_lines(caller_text, error_rate, error_detail, prior,
                                      strand_min_alt, strand_min_per_strand,
                                      homopolymer_min, limit)
    # Both sides explicit: the read tolerated a non-ASCII header byte, but the
    # U+FFFD it produced could not then be written under an ASCII locale.
    with open(raw_vcf, encoding='utf-8', errors='replace') as src, \
            open(out_vcf, 'w', encoding='utf-8', errors='replace') as out:
        for line in src:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('##'):
                print(line, file=out)
                continue
            if line.startswith('#'):
                for added in header_added:
                    print(added, file=out)
                print(line, file=out)
                continue
            fields = line.split('\t')
            if len(fields) < 8:
                continue
            counts['records'] += 1
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            alts = alt.split(',')
            if len(ref) != 1 or alt == '.' or any(len(a) != 1 for a in alts):
                counts['not_snp'] += 1
                continue
            try:
                raw_qual = float(fields[5])
            except ValueError:
                raw_qual = 0.0
            if not raw_qual > caller_min_qual:
                counts['below_caller_qual'] += 1
                continue
            info = fields[7]
            fmt = fields[8] if len(fields) > 8 else ''
            sample = fields[9] if len(fields) > 9 else ''
            dp4 = _dp4(info, fmt, sample)
            if dp4 is None:
                counts['no_depth'] += 1
                continue
            ref_fwd, ref_rev, alt_fwd, alt_rev = dp4
            evq = evidence_qual(ref_fwd, ref_rev, alt_fwd, alt_rev, error_rate, prior)
            flags = []
            if one_strand_only(ref_fwd, ref_rev, alt_fwd, alt_rev, strand_min_alt, strand_min_per_strand):
                flags.append(FLAG_STRAND)
                counts['flag_strand'] += 1
            seq = seqs.get(chrom)
            if seq is None:
                counts['contig_missing'] += 1
            else:
                try:
                    index0 = int(pos) - 1
                except ValueError:
                    index0 = -1
                if homopolymer_touching(seq, index0) >= homopolymer_min:
                    flags.append(FLAG_HOMOPOLYMER)
                    counts['flag_homopolymer'] += 1
            # Never below the N band: bcftools reported a variant here, so the
            # sample is uncertain at the least, and step 2 must not show the
            # reference base.  Counted when the floor raised the value.
            qual = max(evq, float(N_THRESHOLD))
            if qual > evq:
                counts['floored'] += 1
            if flags:
                # In the N band whatever the evidence: the sample is uncertain
                # here, so neither a SNP nor a reference call is asserted.
                qual = float(min(qual, limit))
                counts['limited'] += 1
            # The genotype from the evidence model, for a record with one ALT.
            # bcftools' AC and GT are rewritten when they disagree, and its
            # genotype is kept in INFO.  A multi-allelic record is left as it is.
            ac = _info_value(info, 'AC')
            caller_gt = None
            if len(alts) == 1 and ac in ('1', '2'):
                want = '1' if mixed_by_evidence(ref_fwd, ref_rev, alt_fwd, alt_rev, error_rate) else '2'
                if want != ac:
                    counts['genotype_changed'] += 1
                    info = re.sub(r'(^|;)AC=[^;]*', lambda m: f'{m.group(1)}AC={want}', info, count=1)
                    if fmt and sample and fmt.split(':')[0] == 'GT':
                        parts = sample.split(':')
                        caller_gt = parts[0]
                        parts[0] = '0/1' if want == '1' else '1/1'
                        fields[9] = ':'.join(parts)
                    ac = want
            if ac == '1':
                counts['ac1'] += 1
            elif ac == '2':
                counts['ac2'] += 1
            extra = f'VSNP_RAWQ={raw_qual:g};VSNP_EVQ={evq:g}'
            if caller_gt:
                extra += f';VSNP_CALLER_GT={caller_gt}'
            if flags:
                extra += ';VSNP_FLAG=' + ','.join(flags)
            fields[5] = f'{qual:g}'
            fields[6] = ';'.join(FILTER_NAMES[f] for f in flags) if flags else '.'
            fields[7] = (info + ';' + extra) if info and info != '.' else extra
            print('\t'.join(fields), file=out)
            counts['written'] += 1
    return counts


def run_longread_calling(reference, bam, raw_vcf, out_vcf, threads=8, log=None,
                         error_rate=None):
    '''
    The whole nanopore calling stage: call, measure, model.  Returns a dict
    with the caller text, error rate and its source, the per-read weight and
    the counts from apply_model, and appends the system calls to log.
    '''
    caller_text = call_variants(reference, bam, raw_vcf, threads=threads, log=log)
    if error_rate is None:
        rate, detail = measure_error_rate(bam)
        if log is not None:
            log.append(f'samtools stats {bam}')
    else:
        rate = min(max(float(error_rate), LONGREAD_ERROR_RATE_MIN), LONGREAD_ERROR_RATE_MAX)
        detail = f'{rate} given on the command line'
    counts = apply_model(raw_vcf, out_vcf, reference, rate, detail, caller_text)
    return {'caller': caller_text, 'error_rate': rate, 'error_rate_detail': detail,
            'weight_per_read': per_read_weight(rate), 'counts': counts,
            'limit_qual': QUAL_THRESHOLD - 1}


def describe(result):
    '''One-line summaries of a run_longread_calling result, for the run log and report.'''
    c = result['counts']
    return [
        f"Long-read SNP caller: {result['caller']}",
        f"Long-read error rate: {result['error_rate_detail']}; evidence per clean read "
        f"{result['weight_per_read']:.1f} (Phred), variant prior {LONGREAD_VARIANT_PRIOR}",
        f"Long-read calls: {c['records']} records from bcftools; {c['not_snp']} not single-base "
        f"SNPs and {c['below_caller_qual']} with caller QUAL <= {LONGREAD_CALLER_MIN_QUAL} "
        f"removed; {c['written']} written with the evidence QUAL ({c['ac2']} AC=2, {c['ac1']} AC=1); "
        f"{c['flag_strand']} flagged one-strand, {c['flag_homopolymer']} flagged homopolymer "
        f"(>= {LONGREAD_HOMOPOLYMER_MIN} bp): {c['limited']} placed in the N band "
        f"(QUAL {N_THRESHOLD}-{result['limit_qual']}, N in step 2); {c['floored']} thin calls raised "
        f"to the floor of QUAL {N_THRESHOLD} (N in step 2, never the reference); genotype decided by the "
        f"evidence model, {c['genotype_changed']} differing from bcftools (its genotype kept as VSNP_CALLER_GT)",
    ]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='vsnp3_longread.py', formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
        ---------------------------------------------------------
        Call SNPs from nanopore reads aligned with minimap2, and write a VCF whose
        QUAL is on the same scale as freebayes' so that vSNP3's thresholds apply.

        Usage:
          vsnp3_longread.py -r reference.fasta -b sample_nodup.bam
          vsnp3_longread.py -r reference.fasta -b sample_nodup.bam --vcf raw.vcf -o out.vcf

        The first form runs bcftools mpileup/call and the model; the second applies
        the model to a VCF that bcftools call has already written.  The BAM is still
        needed for the error rate unless --error_rate is given.
        ---------------------------------------------------------
        '''), epilog='---------------------------------------------------------')
    parser.add_argument('-r', '--reference', dest='reference', required=True,
                        help='reference FASTA the reads were aligned to')
    parser.add_argument('-b', '--bam', dest='bam', required=True, help='sorted, indexed BAM')
    parser.add_argument('--vcf', dest='vcf', default=None,
                        help='apply the model to this bcftools VCF instead of calling')
    parser.add_argument('-o', '--output', dest='output', default=None,
                        help='output VCF (default: <sample>_filtered_hapall_nanopore.vcf)')
    parser.add_argument('-n', '--SAMPLE_NAME', dest='sample_name', default=None,
                        help='sample name for the output file (default: from the BAM name)')
    parser.add_argument('-t', '--threads', dest='threads', type=int, default=8)
    parser.add_argument('--error_rate', dest='error_rate', type=float, default=None,
                        help='per-read error rate to use instead of measuring it from the BAM')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug',
                        help='keep the raw bcftools VCF')
    parser.add_argument('-v', '--version', action='version',
                        version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()

    sample = args.sample_name or re.sub('[_.].*', '', os.path.basename(args.bam))
    output = args.output or f'{sample}_filtered_hapall_nanopore.vcf'
    log = []
    if args.vcf:
        if args.error_rate is None:
            rate, detail = measure_error_rate(args.bam)
        else:
            rate = min(max(args.error_rate, LONGREAD_ERROR_RATE_MIN), LONGREAD_ERROR_RATE_MAX)
            detail = f'{rate} given on the command line'
        caller_text = f'bcftools {bcftools_version()} (VCF supplied with --vcf)'
        counts = apply_model(args.vcf, output, args.reference, rate, detail, caller_text)
        result = {'caller': caller_text, 'error_rate': rate, 'error_rate_detail': detail,
                  'weight_per_read': per_read_weight(rate), 'counts': counts,
                  'limit_qual': QUAL_THRESHOLD - 1}
    else:
        raw = f'{sample}_unfiltered_hapall.vcf'
        result = run_longread_calling(args.reference, args.bam, raw, output,
                                      threads=args.threads, log=log, error_rate=args.error_rate)
        if not args.debug:
            os.remove(raw)
    for call in log:
        print(f'SYSTEM CALL: {call}')
    for line in describe(result):
        print(line)
    print(f'Wrote {output}')
    sys.exit(0)

# Created October 2026 by Tod Stuber
