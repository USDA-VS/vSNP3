#!/usr/bin/env python3

from vsnp3_version import __version__

import os
import json
import shlex
import subprocess
import re
import glob
import shutil
import locale
import argparse
import textwrap
from datetime import datetime
import pandas as pd
from Bio import SeqIO

import vsnp3_run
import vsnp3_longread
from vsnp3_version import LONGREAD_VARIANT_PRIOR, LONGREAD_HOMOPOLYMER_MIN
from vsnp3_version import LONGREAD_CALLER_MIN_QUAL
from vsnp3_file_setup import Setup
from vsnp3_file_setup import Excel_Stats
from vsnp3_file_setup import worker_count
from vsnp3_fastq_stats_seqkit import FASTQ_Stats
from vsnp3_vcf_annotation import VCF_Annotation
from vsnp3_assembly import Assemble
from vsnp3_zero_coverage import Zero_Coverage

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


# Threads for bwa, minimap2, samtools, bcftools and freebayes-parallel.  Every
# call used to say 8 whatever the machine: step 1 never read VSNP3_MAX_CPUS,
# which the vSNP GUI sets to the session's core budget, and the cluster
# wrappers start several samples at once.  Now min(8, worker_count()): still 8
# on a large machine, fewer when the cap or the core count says so.
DEFAULT_THREADS = 8

# bwa mem reads its input in batches of 10 Mbp per thread and estimates the
# insert-size distribution batch by batch, so a paired-end alignment depends on
# the thread count (measured: 160 of 600,000 read pairs placed differently at
# -t 2 than at -t 8).  vSNP3 has always run bwa with 8 threads; when the count
# is lowered the batch is pinned to the 8-thread size with -K, which gives the
# same SAM as an uncapped run.
BWA_BATCH_BASES = 10_000_000 * DEFAULT_THREADS


def step1_threads():
    '''min(DEFAULT_THREADS, worker_count()), never below 1.'''
    return max(1, min(DEFAULT_THREADS, worker_count()))


def move_replacing(src, dest_dir):
    '''
    Move src into dest_dir, replacing a file or directory of the same name.

    A re-run in the same folder used to stop here with shutil.Error
    "Destination path ... already exists" after the whole alignment had
    completed.  Returns the destination path.
    '''
    dest = os.path.join(dest_dir, os.path.basename(src.rstrip(os.sep)))
    if os.path.isdir(dest) and not os.path.islink(dest):
        shutil.rmtree(dest)
    elif os.path.lexists(dest):
        os.remove(dest)
    shutil.move(src, dest)
    return dest


def concatenated_reference_name(FASTAs, directory='.'):
    '''
    Name of the concatenated reference written into the working directory.

    <stems joined by _>.fasta, as always -- unless a file of that name is
    already there, which is the user's own FASTA given with -r from the
    working directory: concat_fasta then rewrote it in place and the run moved
    it into alignment_<ref>/ at the end, so the second sample run in that
    folder could not find its reference.  Then <stems>.vsnp3_reference.fasta
    is used instead.  The stem before the first '.' is the same either way,
    so alignment_<ref>/, the reference name in the statistics and the sample's
    file names do not change.
    '''
    stems = [re.sub(r'\..*', '', os.path.basename(each)) for each in FASTAs]
    joined = '_'.join(stems)
    plain = f'{joined}.fasta'
    if os.path.lexists(os.path.join(directory, plain)):
        return f'{joined}.vsnp3_reference.fasta'
    return plain


def contigs_with_mapped_reads(bam):
    '''Contigs with at least one mapped read, from samtools idxstats (needs the index).'''
    proc = vsnp3_run.run(['samtools', 'idxstats', bam])
    contigs = []
    for line in vsnp3_run._text(proc.stdout).splitlines():
        parts = line.split('\t')
        if len(parts) >= 3 and parts[0] != '*':
            try:
                if int(parts[2]) > 0:
                    contigs.append(parts[0])
            except ValueError:
                continue
    return contigs


def contigs_without_records(contigs, vcf_path):
    '''The contigs among `contigs` for which the VCF holds no record.'''
    seen = set()
    with open(vcf_path, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            seen.add(line.split('\t', 1)[0])
    return [contig for contig in contigs if contig not in seen]


class Alignment(Setup):
    '''
    '''

    def __init__(self, SAMPLE_NAME=None, FASTQ_R1=None, FASTQ_R2=None, reference=None, nanopore=False, gbk=None, assemble_unmap=None, debug=False):
        '''
        Start at class call
        '''
        Setup.__init__(self, SAMPLE_NAME=SAMPLE_NAME, FASTA=reference, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.print_run_time('Align and make VCF file')
        self.nanopore = nanopore
        self.gbk = gbk
        self.assemble_unmap = assemble_unmap
        self.threads = step1_threads()
        # (title, [lines]) of external tool output worth keeping -- a failed
        # SPAdes run's log tail, anything freebayes-parallel wrote to stderr.
        # vsnp3_step1.py writes it into the run log; it is not in the HTML
        # report, which carries the command lines and notes only.
        self.tool_output = []

    def filter_vcf(self, input_file, output_file, quality_threshold=LONGREAD_CALLER_MIN_QUAL):
        with open(input_file, 'r', encoding='utf-8', errors='replace') as infile, \
                open(output_file, 'w', encoding='utf-8', errors='replace') as outfile:
            for line in infile:
                if line.startswith('#'):
                    outfile.write(line)
                else:
                    fields = line.strip().split('\t')
                    if len(fields) >= 6:
                        try:
                            qual = float(fields[5])
                            if qual > quality_threshold:
                                outfile.write(line)
                        except ValueError:
                            continue
        
    def _read_stats_from_flagstat(self, bam):
        '''
        The read counts the markdup statistics file supplies on the Illumina
        path, taken from samtools flagstat for a long-read BAM, which is not
        duplicate-marked.  Primary reads only, as markdup counts them.
        '''
        argv = ['samtools', 'flagstat', '-O', 'json', bam]
        proc = vsnp3_run.run(argv)
        try:
            passed = json.loads(vsnp3_run._text(proc.stdout))['QC-passed reads']
            primary = int(passed.get('primary', passed.get('total', 0)))
            mapped = int(passed.get('primary mapped', passed.get('mapped', 0)))
        except (ValueError, KeyError, TypeError) as e:
            raise vsnp3_run.ExternalToolError(argv, detail=f'could not read the counts ({e})') from e
        self.MARKDUPLICATE_COMMAND = 'not run: long reads are not duplicate-marked'
        self.READS_READ = primary
        self.READS_WRITTEN = primary
        self.READS_EXCLUDED = max(0, primary - mapped)        # unmapped
        self.READS_EXAMINED = mapped
        self.READS_PAIRED = 0
        self.READS_SINGLE = mapped
        for name in ('DUPLICATE_PAIR', 'DUPLICATE_SINGLE', 'DUPLICATE_PAIR_OPTICAL',
                     'DUPLICATE_SINGLE_OPTICAL', 'DUPLICATE_NON_PRIMARY',
                     'DUPLICATE_NON_PRIMARY_OPTICAL', 'DUPLICATE_PRIMARY_TOTAL',
                     'DUPLICATE_TOTAL', 'ESTIMATED_LIBRARY_SIZE'):
            setattr(self, name, 0)
        self.DUPLICATION_RATIO = 0.0

    def run(self,):
        '''
        description
        '''
        # A FASTQ_Stats(...).run() used to sit here.  Its result was assigned to a
        # local that nothing else in this file referenced, and FASTQ_Stats.run()
        # makes two full passes over each read file -- `seqkit stats -a` and
        # `seqkit fx2tab | awk`.  So every sample paid four complete gzip
        # decompress-and-scan passes, ~25 ms per MB of compressed reads, for a
        # value that was thrown away.  vsnp3_step1.py already computes and uses
        # the same statistics.
        sample_name = self.sample_name
        reference = self.reference
        threads = self.threads
        samfile = f'{sample_name}.sam'
        fixmate_bamfile = f'{sample_name}_fixmate.bam'
        pos_srt_bamfile = f'{sample_name}_pos_srt.bam'
        nodup_bamfile = f'{sample_name}_nodup.bam'
        unfiltered_hapall = f'{sample_name}_unfiltered_hapall.vcf'
        mapfix_hapall = f'{sample_name}_mapfix_hapall.vcf'
        filtered_hapall = f'{sample_name}_filtered_hapall.vcf'
        unmapped_read1 = f'{sample_name}_unmapped_R1.fastq'
        unmapped_read2 = f'{sample_name}_unmapped_R2.fastq'
        unmapped_read = f'{sample_name}_unmapped.fastq'
        # Working files named for the sample.  markduplicate_stats.txt and
        # chrom_ranges.txt were fixed names, so two samples run by hand in one
        # folder read each other's duplicate statistics and called variants over
        # each other's reference regions.
        markdup_stats = f'{sample_name}_markduplicate_stats.txt'
        chrom_ranges_file = f'{sample_name}_chrom_ranges.txt'
        alignment_vcf_run_summary = []

        def stamp():
            return datetime.now().strftime("%Y-%m-%d_%H:%M:%S")

        # Every external command below is run as an argument list through
        # vsnp3_run.run, so a '(' or ';' in a file name is an ordinary
        # character and a non-zero exit is an error here, not an empty file
        # three stages later.  The command recorded in the run summary is the
        # same text the shell string used to be; stderr goes where the old
        # redirect sent it (the terminal, or /dev/null) unless it is kept.
        vsnp3_run.run(['samtools', 'faidx', reference], stderr=None)
        # os.system(f'picard CreateSequenceDictionary REFERENCE={reference} OUTPUT={reference.rsplit(".", 1)[0]}.dict 2> /dev/null')
        vsnp3_run.run(['bwa', 'index', reference], stderr=subprocess.DEVNULL) # Needed for freebayes-parallel.
        if self.nanopore:
            self.aligner = "Minimap2"
            # map-ont suits every nanopore chemistry.  --secondary=no keeps out of
            # the BAM the secondary alignments bcftools skips anyway, and the read
            # group now says what the reads are (it used to say ILLUMINA).
            read_group = f'@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ONT'
            run_set = f'minimap2 -a -x map-ont --secondary=no -R "{read_group}" -t {threads} {reference} {self.FASTQ_R1} -o {samfile}'
            vsnp3_run.run(['minimap2', '-a', '-x', 'map-ont', '--secondary=no', '-R', read_group,
                           '-t', threads, reference, self.FASTQ_R1, '-o', samfile], stderr=None)
        else:
            self.aligner = "BWA"
            read_group = f'@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA\\tPI:250'
            # See BWA_BATCH_BASES: only when the thread count is not the 8 the
            # alignment has always been made with.
            batch = [] if threads == DEFAULT_THREADS else ['-K', str(BWA_BATCH_BASES)]
            batch_text = ''.join(f' {b}' for b in batch)
            reads = [self.FASTQ_R1, self.FASTQ_R2] if self.paired else [self.FASTQ_R1]
            # run_set = f'minimap2 -t 8 -a -x sr {reference} {self.FASTQ_R1} {self.FASTQ_R2} -o {samfile}'
            run_set = f'bwa mem -M -R "{read_group}" -t {threads}{batch_text} {reference} {" ".join(reads)} > {samfile}'
            vsnp3_run.run(['bwa', 'mem', '-M', '-R', read_group, '-t', threads] + batch
                          + [reference] + reads, stdout_path=samfile, stderr=None)
        alignment_vcf_run_summary.append(f'SYSTEM CALL: {run_set} -- {stamp()}')
        # The exit status is checked by vsnp3_run.run; the artifact is checked
        # too.  An empty SAM is what a failed bwa leaves behind, and every later
        # stage happily processes it.
        vsnp3_run.require_output(samfile, min_bytes=100, what=f'{self.aligner} SAM output')

        if self.nanopore:
            # Long reads are sorted and indexed, not duplicate-marked.  Nanopore
            # libraries are not PCR-amplified, and samtools markdup on single-end
            # reads compares start positions only, so at depth it discarded reads
            # that merely began at the same base.  The file keeps the _nodup name
            # every other tool and wrapper looks for; the statistics and the run
            # log say that no marking was done.
            samsort = f'samtools sort -l 1 -@{threads} -o {nodup_bamfile} {samfile}'
            vsnp3_run.run(['samtools', 'sort', '-l', '1', f'-@{threads}', '-o', nodup_bamfile, samfile],
                          stderr=None)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samsort} -- {stamp()}')
            alignment_vcf_run_summary.append(f'NOTE: Long reads: sorted only, no duplicate marking (nanopore libraries are not PCR-amplified, and markdup on single-end reads discards reads that share a start position) -- {stamp()}')
            vsnp3_run.require_output(nodup_bamfile, bam=True, what='samtools sort output')
        else:
            # http://www.htslib.org/workflow/fastq.html
            samfixmate = f'samtools fixmate -O bam,level=1 -m {samfile} {fixmate_bamfile}'
            samsort = f'samtools sort -l 1 -@{threads} -o {pos_srt_bamfile} {fixmate_bamfile}'
            sammarkup = f'samtools markdup -f {markdup_stats} -r -O bam,level=1 {pos_srt_bamfile} {nodup_bamfile}'
            vsnp3_run.run(['samtools', 'fixmate', '-O', 'bam,level=1', '-m', samfile, fixmate_bamfile],
                          stderr=None)
            vsnp3_run.run(['samtools', 'sort', '-l', '1', f'-@{threads}', '-o', pos_srt_bamfile, fixmate_bamfile],
                          stderr=None)
            vsnp3_run.run(['samtools', 'markdup', '-f', markdup_stats, '-r', '-O', 'bam,level=1',
                           pos_srt_bamfile, nodup_bamfile], stderr=None)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samfixmate} -- {stamp()}')
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samsort} -- {stamp()}')
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {sammarkup} -- {stamp()}')
            vsnp3_run.require_output(fixmate_bamfile, bam=True, what='samtools fixmate output')
            vsnp3_run.require_output(pos_srt_bamfile, bam=True, what='samtools sort output')
            vsnp3_run.require_output(nodup_bamfile, bam=True, what='samtools markdup output')
            vsnp3_run.require_output(markdup_stats,
                                     what='samtools markdup statistics')
        mapped = vsnp3_run.mapped_read_count(nodup_bamfile)
        if mapped == 0:
            raise vsnp3_run.ExternalToolError(
                [self.aligner.lower(), reference, self.FASTQ_R1],
                detail=f'no reads aligned to {os.path.basename(reference)}. Coverage '
                       f'would be 0x everywhere, so no meaningful VCF can be produced. '
                       f'Check that the reference matches the organism sequenced.')
        alignment_vcf_run_summary.append(
            f'NOTE: {mapped} reads mapped -- {datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}')

        if self.nanopore:
            self._read_stats_from_flagstat(nodup_bamfile)
            alignment_vcf_run_summary.append(f'NOTE: Read stats gathered by samtools flagstat; no duplicate marking for long reads -- {stamp()}')
        else:
            # http://www.htslib.org/doc/samtools-markdup.html
            alignment_vcf_run_summary.append(f'NOTE: Read stats gathered by {markdup_stats} -- {stamp()}')
            mdf = pd.read_csv(markdup_stats, delimiter=':', index_col=0, header=None,
                              encoding='utf-8', encoding_errors='replace')
            mdict = mdf.to_dict(orient='dict')[1]
            self.MARKDUPLICATE_COMMAND = mdict['COMMAND']
            self.READS_READ = int(mdict['READ'])
            self.READS_WRITTEN = int(mdict['WRITTEN'])
            self.READS_EXCLUDED = int(mdict['EXCLUDED'])
            self.READS_EXAMINED = int(mdict['EXAMINED'])
            self.READS_PAIRED = int(mdict['PAIRED'])
            self.READS_SINGLE = int(mdict['SINGLE'])
            self.DUPLICATE_PAIR = int(mdict['DUPLICATE PAIR'])
            self.DUPLICATE_SINGLE = int(mdict['DUPLICATE SINGLE'])
            self.DUPLICATE_PAIR_OPTICAL = int(mdict['DUPLICATE PAIR OPTICAL'])
            self.DUPLICATE_SINGLE_OPTICAL = int(mdict['DUPLICATE SINGLE OPTICAL'])
            self.DUPLICATE_NON_PRIMARY = int(mdict['DUPLICATE NON PRIMARY'])
            self.DUPLICATE_NON_PRIMARY_OPTICAL = int(mdict['DUPLICATE NON PRIMARY OPTICAL'])
            self.DUPLICATE_PRIMARY_TOTAL = int(mdict['DUPLICATE PRIMARY TOTAL'])
            self.DUPLICATE_TOTAL = int(mdict['DUPLICATE TOTAL'])
            self.ESTIMATED_LIBRARY_SIZE = int(mdict['ESTIMATED_LIBRARY_SIZE'])
            try:
                self.DUPLICATION_RATIO = float(self.DUPLICATE_TOTAL/self.READS_EXAMINED)
            except ZeroDivisionError:
                self.DUPLICATION_RATIO = 0

        vsnp3_run.run(['samtools', 'index', nodup_bamfile], stderr=None)
        vsnp3_run.require_output(f'{nodup_bamfile}.bai', what='samtools index output')

        if self.nanopore:
            # bcftools mpileup/call with the nanopore profile, then the vSNP3
            # long-read evidence QUAL (vsnp3_longread.py).  bcftools' QUAL is
            # capped near 228, so it never meant what freebayes' QUAL means and
            # the thresholds could not read it; the +100 that used to be added
            # to every nanopore QUAL is gone, and the VCF header records the
            # model.  A header-only VCF is a sample identical to the reference,
            # accepted as on the freebayes path.
            filtered_hapall = f'{sample_name}_filtered_hapall_nanopore.vcf'
            calls = []
            longread = vsnp3_longread.run_longread_calling(
                reference, nodup_bamfile, unfiltered_hapall, filtered_hapall,
                threads=threads, log=calls)
            for call in calls:
                alignment_vcf_run_summary.append(f'SYSTEM CALL: {call} -- {stamp()}')
            for line in vsnp3_longread.describe(longread):
                alignment_vcf_run_summary.append(f'NOTE: {line} -- {stamp()}')
            if longread['counts']['written'] == 0:
                print('\n### bcftools called no SNP: the sample matches the reference '
                      'wherever it has coverage\n')
            self.longread = longread
        else:
            with open(chrom_ranges_file, 'w', encoding='utf-8', errors='replace') as chrom_ranges, \
                    open(reference, encoding='utf-8', errors='replace') as reference_handle:
                for record in SeqIO.parse(reference_handle, "fasta"):
                    chrom = record.id
                    total_len = len(record.seq)
                    min_number = 0
                    step = 100000
                    if step < total_len:
                        for chunk in range(min_number, total_len, step)[1:]:
                            print("{}:{}-{}".format(chrom, min_number, chunk), file=chrom_ranges)
                            min_number = chunk
                    print("{}:{}-{}".format(chrom, min_number, total_len), file=chrom_ranges)
            # freebayes-parallel is a bash script that hands everything after the
            # thread count to GNU parallel, which joins it into a command line and
            # runs that through a shell again.  So the reference and BAM paths are
            # shell-quoted here; for an ordinary name shlex.quote adds nothing,
            # and for 'sample(1)_nodup.bam' it is what lets freebayes run.  The
            # regions file is read by the script itself (cat "$regionsfile"), so
            # it is passed as is.
            reference_q, nodup_q = shlex.quote(reference), shlex.quote(nodup_bamfile)
            freebayes_parallel = f'freebayes-parallel {chrom_ranges_file} {threads} -E -1 -e 1 -u --strict-vcf -f {reference_q} {nodup_q} > {unfiltered_hapall}' #2>/dev/null'
            # stderr is kept: freebayes-parallel's exit status is vcfuniq's, so a
            # region whose freebayes failed (killed for memory, say) is reported
            # only on stderr, and the VCF simply has no record there -- which
            # step 2 reads as "matches the reference".
            freebayes_proc = vsnp3_run.run(
                ['freebayes-parallel', chrom_ranges_file, threads, '-E', '-1', '-e', '1', '-u',
                 '--strict-vcf', '-f', reference_q, nodup_q],
                stdout_path=unfiltered_hapall, stderr=subprocess.PIPE)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {freebayes_parallel} -- {stamp()}')
            freebayes_stderr = vsnp3_run._text(freebayes_proc.stderr).strip().splitlines()
            if freebayes_stderr:
                self.tool_output.append(
                    (f'freebayes-parallel stderr ({len(freebayes_stderr)} lines, last {min(30, len(freebayes_stderr))} shown)',
                     freebayes_stderr[-30:]))
            # An empty VCF is what a failed run leaves behind.
            # A VCF with its #CHROM header and no record is not that: the reads
            # mapped (checked above) and freebayes called no variant, which is
            # what reads identical to the reference give.  That used to stop
            # step 1 ("contains no variant records") with no VCF written, so such
            # a sample could never reach step 2.  It now goes on to a VCF of
            # zero-coverage records only, and the report's no-SNP notice.
            has_records = vsnp3_run.require_output(unfiltered_hapall, vcf_header=True,
                                                   what='freebayes output')
            if not has_records:
                print('\n### freebayes called no variant: the sample matches the reference '
                      'wherever it has coverage\n')
                alignment_vcf_run_summary.append(
                    f'NOTE: freebayes called no variant; the VCF holds zero-coverage records '
                    f'only -- {stamp()}')
            else:
                # A contig with mapped reads and no record at all is either
                # identical to the reference along its length or a region whose
                # freebayes failed silently (see above).  Said, not failed on.
                try:
                    silent = contigs_without_records(contigs_with_mapped_reads(nodup_bamfile),
                                                     unfiltered_hapall)
                except Exception as e:                           # noqa: BLE001
                    silent = []
                    print(f'Note: could not check the VCF for contigs without records: {e}')
                for contig in silent:
                    message = (f'NOTE: freebayes-parallel wrote no record for contig {contig}, '
                               f'which has mapped reads: either the sample matches the reference '
                               f'along it, or a freebayes region failed (freebayes-parallel does '
                               f'not report that; see its stderr in the run log)')
                    print(f'\n### {message}\n')
                    alignment_vcf_run_summary.append(f'{message} -- {stamp()}')
            # os.system(f'bcftools mpileup --threads 8 -Ou -f {reference} {nodup_bamfile} | bcftools call --threads 8 -mv -Ov -o {unfiltered_hapall}')
            with open(mapfix_hapall, 'w+', encoding='utf-8', errors='replace') as write_fix, \
                    open(unfiltered_hapall, 'r', encoding='utf-8', errors='replace') as unfiltered:
                for line in unfiltered:
                    line = line.strip()
                    new_line = re.sub(r';MQM=', r';MQ=', line)
                    new_line = re.sub(r'ID=MQM,', r'ID=MQ,', new_line)
                    print(new_line, file=write_fix)
            # remove clearly poor positions: the same QUAL gate as the long-read
            # caller's (LONGREAD_CALLER_MIN_QUAL, 20)
            self.filter_vcf(mapfix_hapall, filtered_hapall, quality_threshold=LONGREAD_CALLER_MIN_QUAL)
            alignment_vcf_run_summary.append(f'NOTE: Freebayes MQM= changed to MQ= to stay in the same naming as GATK -- {stamp()}')

        unmapped_dir = 'unmapped_reads'
        if not os.path.exists(unmapped_dir):
            os.makedirs(unmapped_dir)
        if self.paired and not self.nanopore:
            samtools_fastq_unmapped = f'samtools fastq -f4 -1 {unmapped_read1} -2 {unmapped_read2} --reference {reference} --threads {threads} {nodup_bamfile} 2> /dev/null'
            vsnp3_run.run(['samtools', 'fastq', '-f4', '-1', unmapped_read1, '-2', unmapped_read2,
                           '--reference', reference, '--threads', threads, nodup_bamfile],
                          stderr=subprocess.DEVNULL)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samtools_fastq_unmapped} -- {stamp()}')
            vsnp3_run.run(['pigz', unmapped_read1], stderr=None)
            vsnp3_run.run(['pigz', unmapped_read2], stderr=None)
            move_replacing(f'{unmapped_read1}.gz', unmapped_dir)
            move_replacing(f'{unmapped_read2}.gz', unmapped_dir)
            self.unmapped_read1 = f'{self.cwd}/{unmapped_dir}/{unmapped_read1}.gz'
            self.unmapped_read2 = f'{self.cwd}/{unmapped_dir}/{unmapped_read2}.gz'
            self.unmapped_read_list = [self.unmapped_read1, self.unmapped_read2]
            if self.assemble_unmap:
                assemble = Assemble(FASTQ_R1=self.unmapped_read1, FASTQ_R2=self.unmapped_read2, debug=True)
        else:
            # Single-end Illumina reads and long reads alike: the unmapped reads
            # are written out, so a nanopore sample's folder looks like any other.
            samtools_fastq_unmapped = f'samtools fastq -f4 -0 {unmapped_read} --reference {reference} --threads {threads} {nodup_bamfile} 2> /dev/null'
            vsnp3_run.run(['samtools', 'fastq', '-f4', '-0', unmapped_read,
                           '--reference', reference, '--threads', threads, nodup_bamfile],
                          stderr=subprocess.DEVNULL)
            alignment_vcf_run_summary.append(f'SYSTEM CALL: {samtools_fastq_unmapped} -- {stamp()}')
            vsnp3_run.run(['pigz', unmapped_read], stderr=None)
            move_replacing(f'{unmapped_read}.gz', unmapped_dir)
            self.unmapped_read = f'{self.cwd}/{unmapped_dir}/{unmapped_read}.gz'
            self.unmapped_read_list = [self.unmapped_read]
            if self.assemble_unmap and not self.nanopore:
                assemble = Assemble(FASTQ_R1=self.unmapped_read, debug=True)

        if not self.assemble_unmap or self.nanopore:
            assembly_message = "skipped assembly"
            alignment_vcf_run_summary.append(f'NOTE: Skipped unmapped read assembly -- {stamp()}')
            assembly_error = f'{self.cwd}/{unmapped_dir}/skipped_assembly'
            self.unmapped_assemble = None
        else:
            assemble.run()
            try:
                shutil.move(assemble.FASTA, f'{assemble.fastq_name}_unmapped.fasta')
                move_replacing(f'{assemble.fastq_name}_unmapped.fasta', unmapped_dir)
                assemble.FASTA = f'{self.cwd}/{unmapped_dir}/{assemble.fastq_name}_unmapped.fasta'
                assemble.stats(assemble.FASTA)
                self.assemble = assemble
                self.unmapped_assemble = assemble.FASTA
                if not self.debug and os.path.isdir('spades_assembly'):
                    shutil.rmtree('spades_assembly')
                assembly_message = f'{assemble.contig_count:,}'
                alignment_vcf_run_summary.append(f'NOTE: SPAdes unmapped read assembly completed -- {stamp()}')
            except TypeError:
                # assemble.FASTA is None: SPAdes produced no scaffolds, failed,
                # or is not installed at all.  Recorded; the alignment
                # goes on to the VCF.  SPAdes used to exit the whole run when
                # it was missing, and its output was discarded; the last lines
                # of it now go to the run log.
                self.unmapped_assemble = None
                assembly_error = f'{self.cwd}/{unmapped_dir}/failed_assembly'
                with open(assembly_error, "w", encoding='utf-8', errors='replace') as opened_file:
                    print("see unmapped FASTQ files for troubleshooting", file=opened_file)
                if os.path.exists('spades_assembly/spades.log'):
                    move_replacing('spades_assembly/spades.log', unmapped_dir)
                if not self.debug and os.path.isdir('spades_assembly'):
                    shutil.rmtree('spades_assembly')
                not_found = 'not found' in str(getattr(assemble, 'spades_version', '') or '')
                assembly_message = 'failed assembly (SPAdes not found)' if not_found else 'failed assembly'
                why = ' (spades.py not found on PATH)' if not_found else ''
                alignment_vcf_run_summary.append(f'NOTE: SPAdes unmapped read assembly failed{why} -- {stamp()}')
                log_tail = getattr(assemble, 'spades_log_tail', None)
                if log_tail:
                    self.tool_output.append(('SPAdes (unmapped read assembly failed), last lines', list(log_tail)))
        self.assembly_message = assembly_message

        if self.gbk:
            self.vcf_annotation = VCF_Annotation(gbk_list=self.gbk, vcf_file=filtered_hapall)
            alignment_vcf_run_summary.append(f'IMPORT: VCF_Annotation(gbk_list=self.gbk, vcf_file=filtered_hapall) -- {stamp()}')

        # The sample name is passed on: Zero_Coverage used to derive it from the
        # BAM name again, so -n with a '_' or '.' named the VCF differently from
        # the BAM, the report and the statistics.
        self.zero_coverage = Zero_Coverage(FASTA=reference, bam=nodup_bamfile, vcf=filtered_hapall,
                                           sample_name=sample_name)
        alignment_vcf_run_summary.append(f'IMPORT: Zero_Coverage(FASTA=reference, bam=nodup_bamfile, vcf=filtered_hapall,) -- {stamp()}')
        if self.nanopore:
            nanopore_zc_name = re.sub('_zc.vcf', '_nanopore_zc.vcf', self.zero_coverage.zero_coverage_vcf)
            os.replace(self.zero_coverage.zero_coverage_vcf, nanopore_zc_name)
            self.zero_coverage.zero_coverage_vcf = nanopore_zc_name

        # Deliverables into alignment_<ref>/, by name, replacing what a previous
        # run of this sample left there.  This used to glob *_nodup.bam, *_zc.vcf,
        # *_nodup.bam.bai and *_annotated.vcf from the whole working directory,
        # whichever sample they belonged to, and shutil.move raised on a file
        # already present in alignment_<ref>/.
        alignment = f'alignment_{self.reference_name}'
        os.makedirs(alignment, exist_ok=True)
        if os.path.exists(unmapped_dir):
            move_replacing(unmapped_dir, alignment)
        deliverables = [nodup_bamfile, f'{nodup_bamfile}.bai', self.zero_coverage.zero_coverage_vcf]
        if self.gbk:
            deliverables.append(os.path.basename(self.vcf_annotation.vcf))
        for each in deliverables:
            if os.path.exists(each):
                move_replacing(each, alignment)
        prefix = glob.escape(sample_name)
        for pattern in (f'{prefix}*_nodup.bam', f'{prefix}*_zc.vcf', f'{prefix}*_nodup.bam.bai',
                        f'{prefix}*_annotated.vcf'):
            for each in sorted(glob.glob(pattern)):
                move_replacing(each, alignment)
        self.zero_coverage_vcf_file_path = f'{self.cwd}/{alignment}/{os.path.basename(self.zero_coverage.zero_coverage_vcf)}'
        move_replacing(reference, alignment)
        if os.path.exists(f'{reference}.fai'):
            move_replacing(f'{reference}.fai', alignment)
        self.reference = f'{self.cwd}/{alignment}/{os.path.basename(reference)}'
        self.nodup_bamfile = f'{self.cwd}/{alignment}/{nodup_bamfile}'
        self.zc_vcf = self.zero_coverage.zero_coverage_vcf
        if self.gbk:
            self.annotated_vcf = f'{self.cwd}/{alignment}/{os.path.basename(self.vcf_annotation.vcf)}'

        # Working files, this sample's and this reference's only.  The patterns
        # used to be unprefixed ('*.sam', '*.bai', '*.fai', '*_unfiltered_hapall.vcf',
        # ...) and swept every other sample's files in the directory into temp/
        # and, without --debug, deleted them.
        temp_dir = './temp'
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        reference_base = glob.escape(os.path.basename(reference))
        reference_stem = glob.escape(os.path.basename(reference).rsplit('.', 1)[0])
        files_grab = []
        for files in (f'{prefix}_unmapped*.fastq.gz', f'{prefix}_all.bam', f'{prefix}_fixmate.bam',
                      f'{prefix}_pos_srt.bam', markdup_stats, f'{prefix}*.bai',
                      f'{prefix}_filtered_hapall.vcf', f'{prefix}_mapfix_hapall.vcf',
                      f'{prefix}_unfiltered_hapall.vcf', f'{prefix}_filtered_hapall_nanopore.vcf',
                      f'{prefix}.sam', f'{reference_base}.amb', f'{reference_base}.ann',
                      f'{reference_base}.bwt', f'{reference_base}.pac', f'{reference_base}.sa',
                      f'{prefix}_sorted.bam', f'{reference_stem}.dict', chrom_ranges_file,
                      f'{reference_base}.fai', 'dup_metrics.csv'):
            files_grab.extend(sorted(glob.glob(files)))
        for each in files_grab:
            move_replacing(each, temp_dir)

        if self.debug is False:
            shutil.rmtree(temp_dir)
            alignment_vcf_run_summary.append(f'NOTE: Files moved to temp_dir and removed: *_unmapped*.fastq.gz, *_all.bam, *_fixmate.bam, *_pos_srt.bam, {markdup_stats}, *.bai, *_filtered_hapall.vcf, *_mapfix_hapall.vcf, *_unfiltered_hapall.vcf, *_filtered_hapall_nanopore.vcf, *.sam, *.amb, *.ann, *.bwt, *.pac, *.fasta.sa, *_sorted.bam, *.dict, {chrom_ranges_file}, *.fai, dup_metrics.csv -- {stamp()}')

        #get program versions
        programs = []
        try:
            # Use more modern subprocess calls with error handling
            try:
                minimap2_version = subprocess.run(['minimap2', "--version"], capture_output=True, text=True, check=False).stdout.rstrip()
                programs.append(f'Minimap2: {minimap2_version}')
            except (subprocess.SubprocessError, FileNotFoundError):
                programs.append('Minimap2: version unknown')

            try:
                bcftools_output = subprocess.run(['bcftools', '--version'], capture_output=True, text=True, check=False).stdout
                first = bcftools_output.strip().splitlines()[0] if bcftools_output.strip() else ''
                programs.append(f'Bcftools: {first.split()[1]}' if len(first.split()) > 1 else 'Bcftools: version unknown')
            except (subprocess.SubprocessError, FileNotFoundError, IndexError):
                programs.append('Bcftools: version unknown')

            try:
                freebayes_output = subprocess.run(['freebayes'], capture_output=True, text=True, check=False).stdout
                for line in freebayes_output.splitlines():
                    if 'version' in line:
                        programs.append(f'Freebayes: {line.rstrip().replace("version:  ", "")}')
            except (subprocess.SubprocessError, FileNotFoundError):
                programs.append('Freebayes: version unknown')
                
            try:
                samtools_output = subprocess.run(['samtools', 'version'], capture_output=True, text=True, check=False).stdout.splitlines()
                if len(samtools_output) >= 2:
                    programs.append(samtools_output[0])
                    programs.append(samtools_output[1])
                else:
                    # Fallback for older samtools versions
                    samtools_stderr = subprocess.run(['samtools'], capture_output=True, text=True, check=False).stderr
                    for line in samtools_stderr.splitlines():
                        if 'Version' in line:
                            programs.append(f'Samtools: {line}')
            except (subprocess.SubprocessError, FileNotFoundError, IndexError):
                programs.append('Samtools: version unknown')
                
            # `assemble` is a local of run(), only bound when --assemble_unmap was
            # used, so referencing it here raised UnboundLocalError on every default
            # run.  The generic handler below turned that into an "Error getting
            # program versions" line, which vsnp3_step1.py then filtered out of the
            # run log by matching on the message text -- hiding the symptom rather
            # than the cause.  run() records the object on self, so ask for that.
            spades_version = getattr(getattr(self, 'assemble', None), 'spades_version', None)
            if spades_version:
                programs.append(spades_version)
        except Exception as e:
            programs.append(f'Error getting program versions: {str(e)}')
        
        self.programs = programs
        self.alignment_vcf_run_summary = alignment_vcf_run_summary
        # Which caller made the VCF, for the statistics and the report: a VCF
        # from either platform can sit in the same step 2 run, and the reader
        # should not have to open the file to tell them apart.
        if self.nanopore:
            self.snp_caller = (f"{self.longread['caller']}; QUAL = vSNP3 long-read evidence "
                               f"QUAL (vsnp3_longread.py)")
        else:
            freebayes = next((p.replace('Freebayes: ', 'freebayes ') for p in programs
                              if p.startswith('Freebayes')), 'freebayes')
            self.snp_caller = f'{freebayes} (freebayes-parallel -E -1 -e 1 -u --strict-vcf)'


    def excel(self, excel_dict):
        # Computed here rather than read off an attribute: it used to be set as a
        # side effect of building the LaTeX table, so excel() only worked if that
        # ran first.
        total_reads = self.READS_PAIRED + self.READS_SINGLE + self.READS_EXCLUDED
        try:
            self.freq_unmapped_reads = self.READS_EXCLUDED / total_reads
        except ZeroDivisionError:
            self.freq_unmapped_reads = 0
        excel_dict['Aligner'] = f'{self.aligner}'
        excel_dict['Mapped Paired Reads'] = f'{self.READS_PAIRED:,}'
        excel_dict['Mapped Single Reads'] = f'{self.READS_SINGLE:,}'
        excel_dict['Unmapped Reads'] = f'{self.READS_EXCLUDED:,}'
        excel_dict['Unmapped Percent'] = f'{(self.freq_unmapped_reads*100):,.1f}%'
        excel_dict['Unmapped Assembled Contigs'] = f'{self.assembly_message}'
        if self.nanopore:
            excel_dict['Duplicate Paired Reads'] = 'not marked (long reads)'
            excel_dict['Duplicate Single Reads'] = 'not marked (long reads)'
            excel_dict['Duplicate Percent of Mapped Reads'] = 'not marked (long reads)'
        else:
            excel_dict['Duplicate Paired Reads'] = f'{self.DUPLICATE_PAIR:,}'
            excel_dict['Duplicate Single Reads'] = f'{self.DUPLICATE_SINGLE:,}'
            excel_dict['Duplicate Percent of Mapped Reads'] = f'{(self.DUPLICATION_RATIO*100):,.1f}%'
        excel_dict['SNP Caller'] = getattr(self, 'snp_caller', '')
        longread = getattr(self, 'longread', None)
        if longread:
            c = longread['counts']
            excel_dict['Long-read Error Rate'] = (f"{longread['error_rate'] * 100:.2f}% "
                                                  f"({longread['error_rate_detail']})")
            excel_dict['Long-read QUAL Model'] = (
                f"evidence QUAL, {longread['weight_per_read']:.1f} per clean read, variant prior "
                f"{LONGREAD_VARIANT_PRIOR}; replaces bcftools QUAL, kept as INFO/VSNP_RAWQ")
            excel_dict['Long-read Calls Limited to N'] = (
                f"{c['limited']:,} of {c['written']:,} ({c['flag_strand']:,} one-strand, "
                f"{c['flag_homopolymer']:,} homopolymer >= {LONGREAD_HOMOPOLYMER_MIN} bp); "
                f"{c.get('floored', 0):,} thin calls raised to the N floor; genotype from the model, "
                f"{c.get('genotype_changed', 0):,} differing from bcftools")
        excel_dict['BAM/Reference File'] = f'{self.zero_coverage.bam} made with {self.reference_name}'
        excel_dict['Reference Length'] = f'{self.zero_coverage.reference_length:,}'
        excel_dict['Genome with Coverage'] = f'{(self.zero_coverage.genome_coverage*100):,.2f}%'
        excel_dict['Average Depth'] = f'{self.zero_coverage.ave_coverage:,.1f}'
        excel_dict['No Coverage Bases'] = f'{self.zero_coverage.total_zero_coverage:,}'
        excel_dict['Percent Ref with Zero Coverage'] = f'{self.zero_coverage.percent_ref_with_zero_coverage:,.6f}%'
        excel_dict['Ambiguous SNPs'] = f'{self.zero_coverage.ac1_count:,}'
        excel_dict['Quality SNPs'] = f'{self.zero_coverage.good_snp_count:,}'


if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    Usage:
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta
    alignment_vcf.py -r1 *fastq.gz -r *fasta
    alignment_vcf.py -r1 *_R1*fastq.gz -r2 *_R2*fastq.gz -r *fasta -g *gbk

    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-n', '--SAMPLE_NAME', action='store', dest='SAMPLE_NAME', required=False, help='Force output files to this sample name')
    parser.add_argument('-r1', '--read1', action='store', dest='FASTQ_R1', required=True, help='Required: single read, R1 when Illumina read')
    parser.add_argument('-r2', '--read2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: R2 Illumina read')
    parser.add_argument('-r', '--reference', nargs='+', dest='FASTA', required=True, help='Required: FASTA file to be used as reference.  Multiple can be specified with wildcard')
    parser.add_argument('-o', '--nanopore', action='store_true', dest='nanopore', default=False, help='nanopore (long-read) path: minimap2 map-ont, bcftools calling and the vSNP3 long-read evidence QUAL (vsnp3_longread.py)')
    parser.add_argument('-b', '--gbk', nargs='*', dest='gbk', required=False, default=None, help='Optional: gbk to annotate VCF file.  Multiple can be specified with wildcard')
    parser.add_argument('-assemble_unmap', '--assemble_unmap', action='store_true', dest='assemble_unmap', help='Optional: assemble the unmapped reads with SPAdes (not done by default).  See also vsnp3_assembly.py')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='keep temp file')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.basename(__file__)}: version {__version__}')
    args = parser.parse_args()
    
    print(f'\n{os.path.basename(__file__)} SET ARGUMENTS:')
    print(args)
    print("\n")

    def concat_fasta(FASTAs):
        fastas_used = ", ".join(os.path.basename(each) for each in FASTAs)
        #concatenate FASTA list, under a name that is not a file already here
        concatenated_FASTA = concatenated_reference_name(FASTAs)
        concatenated_TEMP = f'{concatenated_FASTA}.tmp'
        with open(concatenated_TEMP,'wb') as wfd:
            for each in FASTAs:
                with open(each,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
        shutil.move(concatenated_TEMP, concatenated_FASTA)
        return concatenated_FASTA, fastas_used
    concatenated_FASTA, fastas_used = concat_fasta(args.FASTA) # -r option for FASTA will take wildcard for multiple FASTAs

    alignment = Alignment(SAMPLE_NAME=args.SAMPLE_NAME, FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, reference=concatenated_FASTA, nanopore=args.nanopore, gbk=args.gbk, assemble_unmap=args.assemble_unmap, debug=args.debug)
    alignment.run()

    #Excel Stats
    excel_stats = Excel_Stats(alignment.sample_name)
    excel_stats.excel_dict["FASTA/s"] = fastas_used
    alignment.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    prefix = glob.escape(alignment.sample_name)
    files_grab = []
    for files in ('*.log', '*out', f'{prefix}_all.bam', f'{prefix}*.bai', f'{prefix}_filtered_hapall.vcf', f'{prefix}_mapfix_hapall.vcf', f'{prefix}_unfiltered_hapall.vcf', f'{prefix}.sam', f'{prefix}_sorted.bam', f'{prefix}_chrom_ranges.txt', 'dup_metrics.csv'):
        files_grab.extend(sorted(glob.glob(files)))
    for each in files_grab:
        move_replacing(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Created March 2021 by Tod Stuber