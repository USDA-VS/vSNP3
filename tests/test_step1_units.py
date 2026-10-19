#!/usr/bin/env python3
"""
Unit checks for the step 1 modules: the functions changed for the 2026-10-16
audit (findings H6-H11, M9-M21, L1-L5, G5), each exercised on its failure path.

    python3 tests/test_step1_units.py

Pure Python: no reads, no reference set, no bwa/samtools/freebayes/seqkit/SPAdes.
External programs are mocked where a function would call one.  pandas and
Biopython are needed, as they are by step 1 itself.
"""

import gzip
import io
import os
import shutil
import subprocess
import sys
import tempfile
import threading

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, 'bin')
sys.path.insert(0, BIN)

from vsnp3_version import __version__, MQ_THRESHOLD, QUAL_THRESHOLD  # noqa: E402

PASS, FAIL = [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f'\n        {detail}' if detail and not ok else ''))
    return ok


def skip(name, why):
    print(f'  SKIP  {name}   ({why})')


class scratch:
    '''A temporary directory that is the working directory inside the block.'''

    def __enter__(self):
        self.path = tempfile.mkdtemp(prefix='vsnp3_step1_units_')
        self.previous = os.getcwd()
        os.chdir(self.path)
        return self.path

    def __exit__(self, *exc):
        os.chdir(self.previous)
        shutil.rmtree(self.path, ignore_errors=True)
        return False


def write(path, text, mode='w'):
    with open(path, mode, encoding=None if 'b' in mode else 'utf-8') as f:
        f.write(text)
    return path


# --------------------------------------------------------------------- krona

def section_krona():
    print('\n[1] krona_lca_all: the LCA mapping is parsed, not destroyed (H10)')
    import pandas as pd
    import krona_lca_all

    check('"0:12 1763:5 |:| 0:10" resolves to taxid 1763, not a TypeError',
          krona_lca_all.best_taxid('0:12 1763:5 |:| 0:10') == 1763)
    check('a mapping of only taxid 0 stays 0', krona_lca_all.best_taxid('0:30 |:| 0:12') == 0)
    check('the largest summed count wins across both mates',
          krona_lca_all.best_taxid('9606:3 1763:4 |:| 9606:2 1763:1') == 9606)

    with scratch():
        rows = ['C\tread1\t0\t150|150\t0:12 1763:5 |:| 0:10',
                'C\tread2\t562\t150|150\t562:20 |:| 562:19',
                'U\tread3\t0\t150|150\t0:30 |:| 0:30']
        write('kraken.txt', '\n'.join(rows) + '\n')
        try:
            krona_lca_all.force_tax_number('kraken.txt')
            out = pd.read_csv('kronaInput.txt', sep='\t', header=None)
            by_read = dict(zip(out[0], out[1].astype(int)))
            check('force_tax_number writes kronaInput.txt with read1 forced to 1763',
                  by_read.get('read1') == 1763 and by_read.get('read2') == 562
                  and by_read.get('read3') == 0, f'got {by_read}')
        except Exception as e:                                       # noqa: BLE001
            check('force_tax_number runs on a classified read with taxid 0', False, repr(e))


# -------------------------------------------------------------------- seqkit

# `seqkit stats -a` default output, seqkit 2.10 layout (N50_num present); the
# values are the strings the workbook has always held.
SEQKIT_HUMAN = (
    'file                 format  type  num_seqs     sum_len  min_len  avg_len  max_len   Q1   Q2   Q3  sum_gap  N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n\n'
    '/data/s_R1.fastq.gz  FASTQ   DNA    150,000  34,088,785       35    227.3      251  250  250  251        0  251        1   93.61   89.92    30.14  65.48      0\n')

# seqkit 2.1-2.3 layout: no N50_num column.
SEQKIT_HUMAN_OLD = (
    'file               format  type  num_seqs    sum_len  min_len  avg_len  max_len     Q1     Q2     Q3  sum_gap    N50  Q20(%)  Q30(%)  AvgQual  GC(%)\n'
    '/data/long.fastq   FASTQ   DNA      1,500  2,538,885      903  1,692.6    2,499  1,306  1,699  2,092        0  1,872     100   52.34    26.38  50.06\n')

# What seqkit 2.13 prints with -T for the same reads: Q20(%) and Q30(%) rounded
# to whole numbers.  The reader must not be fed this form.
SEQKIT_T_213 = (
    'file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n\n'
    '/data/s_R1.fastq.gz\tFASTQ\tDNA\t150000\t34088785\t35\t227.3\t251\t218\t250\t251\t0\t250\t2\t94\t90\t24.52\t61.70\t9946\n')


def section_seqkit():
    print('\n[2] vsnp3_fastq_stats_seqkit: columns by header name, values as before (M10, L4)')
    try:
        import vsnp3_fastq_stats_seqkit as sk
    except ImportError as e:
        skip('seqkit parsing', f'import failed: {e}')
        return
    stats = sk.parse_seqkit_stats_table(SEQKIT_HUMAN, 's_R1.fastq.gz')
    expect = {'file_name': '/data/s_R1.fastq.gz', 'read_format': 'FASTQ', 'seq_type': 'DNA',
              'num_seqs': '150,000', 'sum_len': '34,088,785', 'min_len': '35', 'avg_len': '227.3',
              'max_len': '251', 'Q1': '250', 'Q2': '250', 'Q3': '251', 'sum_gap': '0', 'N50': '251',
              'passQ20': '93.61', 'passQ30': '89.92'}
    wrong = {k: (stats.get(k), v) for k, v in expect.items() if stats.get(k) != v}
    check('the seqkit 2.10 default output is read by header name into the stored strings', not wrong, f'{wrong}')

    old = sk.parse_seqkit_stats_table(SEQKIT_HUMAN_OLD, 'long.fastq')
    check('a pre-2.4 layout (no N50_num) gives Q20/Q30 from their own columns',
          (old['passQ20'], old['passQ30'], old['avg_len']) == ('100', '52.34', '1,692.6'),
          f"got {(old['passQ20'], old['passQ30'], old['avg_len'])}")
    check('the old index-based code read AvgQual as Q30 on that layout: now it does not',
          old['passQ30'] != '26.38')
    t213 = sk.parse_seqkit_stats_table(SEQKIT_T_213, 'x')
    check('seqkit 2.13 -T output would store rounded percentages ("94"), which is why -T is not run',
          (t213['passQ20'], t213['passQ30']) == ('94', '90'), f"{t213['passQ20']} {t213['passQ30']}")
    src = open(sk.__file__, encoding='utf-8').read()
    check('seqkit is run without -T', '"-T"' not in src.split('def parse_seqkit_stats_table')[-1].split('class FASTQ_Stats')[-1] and '"seqkit", "stats", "-a", "-o"' in src)

    try:
        sk.parse_seqkit_stats_table('file  format  num_seqs\n/x  FASTQ  5\n', 'x')
        check('a table missing a needed column raises ValueError', False, 'no exception')
    except ValueError as e:
        check('a table missing a needed column raises ValueError', 'lacks column' in str(e), str(e))
    try:
        sk.parse_seqkit_stats_table('', 'x')
        check('empty seqkit output raises ValueError', False, 'no exception')
    except ValueError:
        check('empty seqkit output raises ValueError', True)

    # avg_len with a comma still parses the way vsnp3_step1.py reads it
    check('step 1 can read the long-read avg_len back as a float',
          float(str(old['avg_len']).replace(',', '')) == 1692.6)


# ------------------------------------------------------------- zero coverage

def section_zero_coverage():
    print('\n[3] vsnp3_zero_coverage: sample name from the caller, header-aware VCF read (H7, L1, L4)')
    try:
        import vsnp3_zero_coverage as zc
    except ImportError as e:
        skip('zero_coverage', f'import failed: {e}')
        return
    check('-n B26-0290_S3 keeps the full forced name',
          zc.sample_name_for('B26-0290_S3_nodup.bam', 'B26-0290_S3') == 'B26-0290_S3')
    check('without a caller name the BAM basename up to _ or . is used',
          zc.sample_name_for('B26-0290_S3_nodup.bam') == 'B26-0290')
    check('a BAM given with a path gives a name, not a path',
          zc.sample_name_for('out/run1/13-1941-6_nodup.bam') == '13-1941-6')
    check('a dotted name stops at the dot', zc.sample_name_for('NZ_CP0123.1_nodup.bam') == 'NZ')

    with scratch():
        vcf = ('##fileformat=VCFv4.2\n'
               '##reference=file:///data/réf/NC_002945v4.fasta\n'
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample#1\n'
               'chr1\t100\t.\tA\tG\t200\t.\tAC=2;MQM=60;NOTE=#hash\tGT\t1/1\n'
               '\n'
               'chr1\t200\t.\tC\tT\t50\t.\tAC=1;MQM=58\tGT\t0/1\n')
        write('s.vcf', vcf)
        header, body = zc.read_vcf_lines('s.vcf')
        check('header lines are every # line, #CHROM included',
              len(header) == 3 and header[-1].startswith('#CHROM'), f'{header}')
        check('a non-ASCII header line is read, not fatal', 'réf' in header[1], header[1])
        frame = zc.vcf_body_frame(body)
        check('a record with "#" in INFO is read whole (comment="#" cut it at the #)',
              list(frame['INFO']) == ['AC=2;MQM=60;NOTE=#hash', 'AC=1;MQM=58'], f"{list(frame['INFO'])}")
        check('blank lines are not records', len(frame) == 2)
        check('POS is numeric', str(frame['POS'].dtype).startswith('int'))
        empty = zc.vcf_body_frame([])
        check('a header-only VCF gives an empty frame with the VCF columns',
              empty.empty and list(empty.columns) == zc.VCF_COLUMNS)
        counts = zc.step2_snp_counts(frame)
        check('step2_snp_counts still reads the frame', counts == (1, 1), f'{counts}')


# --------------------------------------------------------------- spoligotype

def section_spoligotype():
    print('\n[4] vsnp3_spoligotype: at least one dask worker, plain FASTQ accepted (M11, M20)')
    try:
        import vsnp3_spoligotype as sp
    except ImportError as e:
        skip('spoligotype', f'optional dependency missing: {e}')
        return
    check('half of 1 CPU is 1 worker, not 0 (the hang)', sp.half_cpus(1) == 1)
    check('half of 2 CPUs is 1 worker', sp.half_cpus(2) == 1)
    check('half of 3 CPUs is 1 worker', sp.half_cpus(3) == 1)
    check('half of 14 CPUs is 7 workers', sp.half_cpus(14) == 7)

    with scratch():
        record = '@r1\nACGTACGT\n+\nIIIIIIII\n'
        write('plain.fastq', record)
        with gzip.open('zipped.fastq.gz', 'wt') as f:
            f.write(record)
        # extension lies: a gzipped file called .fastq and a plain file called .gz
        with gzip.open('lying.fastq', 'wt') as f:
            f.write(record)
        write('lying.fastq.gz', record)
        got = {}
        for name in ('plain.fastq', 'zipped.fastq.gz', 'lying.fastq', 'lying.fastq.gz'):
            try:
                with sp.open_fastq_text(name) as handle:
                    got[name] = handle.read()
            except Exception as e:                                   # noqa: BLE001
                got[name] = repr(e)
        check('open_fastq_text reads plain and gzipped FASTQ whatever the extension says',
              all(v == record for v in got.values()), f'{got}')


# ------------------------------------------------------------------ run_pipe

def section_run():
    print('\n[5] vsnp3_run: a chatty first process cannot deadlock run_pipe; stderr routing (M17, M13)')
    import vsnp3_run

    chatty = [sys.executable, '-c',
              'import sys\n'
              'sys.stderr.write("w" * 300000)\n'       # well past a 64 KB pipe buffer
              'sys.stderr.flush()\n'
              'sys.stdout.write("payload\\n")\n']
    consumer = [sys.executable, '-c', 'import sys; sys.stdout.write(sys.stdin.read())']
    result = {}

    def worker():
        try:
            with scratch():
                codes = vsnp3_run.run_pipe(chatty, consumer, stdout_path='out.txt')
                with open('out.txt') as f:
                    result['out'] = f.read()
                result['codes'] = codes
        except Exception as e:                                       # noqa: BLE001
            result['error'] = repr(e)

    t = threading.Thread(target=worker, daemon=True)
    t.start()
    t.join(timeout=60)
    if t.is_alive():
        check('run_pipe completes with 300 KB on the first process stderr', False, 'still running after 60 s: deadlock')
    else:
        check('run_pipe completes with 300 KB on the first process stderr',
              result.get('codes') == (0, 0) and result.get('out') == 'payload\n', f'{result}')

    try:
        vsnp3_run.run_pipe(chatty, ['definitely-not-a-program-xyz'])
        check('a missing consumer raises ExternalToolError', False, 'no exception')
    except vsnp3_run.ExternalToolError as e:
        check('a missing consumer raises ExternalToolError', 'not found' in str(e), str(e))

    failing = [sys.executable, '-c', 'import sys; sys.stderr.write("boom\\n"); sys.exit(3)']
    try:
        vsnp3_run.run(failing)
        check('run() raises on a non-zero exit', False, 'no exception')
    except vsnp3_run.ExternalToolError as e:
        check('run() raises on a non-zero exit with the stderr tail',
              e.returncode == 3 and 'boom' in str(e), str(e))
    proc = vsnp3_run.run(failing, check=False, stderr=subprocess.DEVNULL)
    check('stderr=DEVNULL discards stderr (the old "2> /dev/null")',
          proc.returncode == 3 and proc.stderr is None)
    proc = vsnp3_run.run([sys.executable, '-c', 'print("hi")'], stderr=None)
    check('stderr=None lets stderr through to the terminal (the old os.system) and still captures stdout',
          proc.stderr is None and vsnp3_run._text(proc.stdout).strip() == 'hi')
    with scratch():
        vsnp3_run.run([sys.executable, '-c', 'print("to file")'], stdout_path='o.txt')
        check('stdout_path replaces "> file"', open('o.txt').read() == 'to file\n')


# ------------------------------------------------------------------ assembly

def section_assembly():
    print('\n[6] vsnp3_assembly: a missing or failing SPAdes is recorded, not sys.exit (H9)')
    try:
        import vsnp3_assembly
    except ImportError as e:
        skip('assembly', f'import failed: {e}')
        return
    Assemble = vsnp3_assembly.Assemble
    argv2 = Assemble.spades_argv(['a_R1.fq.gz', 'a_R2.fq.gz'])
    argv1 = Assemble.spades_argv(['a.fq.gz'])
    check('two read files use --pe1-s/--pe2-s', argv2[:5] == ['spades.py', '--pe1-s', 'a_R1.fq.gz', '--pe2-s', 'a_R2.fq.gz'])
    check('one read file uses --iontorrent -s', argv1[:4] == ['spades.py', '--iontorrent', '-s', 'a.fq.gz'])
    try:
        Assemble.spades_argv([])
        check('no read file raises ValueError, not SystemExit', False, 'no exception')
    except ValueError:
        check('no read file raises ValueError, not SystemExit', True)

    bare = Assemble.__new__(Assemble)
    real_run = vsnp3_assembly.subprocess.run

    def not_found(*a, **k):
        raise FileNotFoundError(2, 'No such file or directory', 'spades.py')

    vsnp3_assembly.subprocess.run = not_found
    try:
        ran = bare._run_spades(argv2)
    except SystemExit:
        ran = 'SystemExit'
    finally:
        vsnp3_assembly.subprocess.run = real_run
    check('spades.py absent: returns False and records "not found" instead of exiting',
          ran is False and bare.spades_version == 'SPAdes: not found on PATH'
          and bare.spades_log_tail and 'not found' in bare.spades_log_tail[0], f'{ran} {bare.__dict__}')

    def failed(*a, **k):
        return subprocess.CompletedProcess(a[0], 1, stdout=b'== Error ==  spades failed\nline 2\n', stderr=b'')

    bare2 = Assemble.__new__(Assemble)
    bare2.spades_version = None
    vsnp3_assembly.subprocess.run = failed
    try:
        ran = bare2._run_spades(argv1)
    finally:
        vsnp3_assembly.subprocess.run = real_run
    check('SPAdes exit 1: returns False and keeps the log tail (stdout when stderr is empty)',
          ran is False and bare2.spades_returncode == 1 and bare2.spades_log_tail == ['== Error ==  spades failed', 'line 2'],
          f'{ran} {bare2.__dict__}')

    with scratch() as here:
        # run(): no spades_assembly/ directory after a failed start must not raise in rmtree
        obj = Assemble.__new__(Assemble)
        obj.sample_name, obj.FASTQ_list, obj.cwd, obj.debug = 's', ['a.fq.gz'], here, False
        obj.print_run_time = lambda *a, **k: None
        vsnp3_assembly.subprocess.run = not_found
        try:
            obj.run()
            ok = obj.FASTA is None and obj.spades_version == 'SPAdes: not found on PATH'
            detail = f'{obj.__dict__}'
        except BaseException as e:                                  # noqa: BLE001
            ok, detail = False, repr(e)
        finally:
            vsnp3_assembly.subprocess.run = real_run
        check('Assemble.run() with no SPAdes leaves FASTA None and returns normally', ok, detail)


# ------------------------------------------------------------ fasta_to_fastq

def section_fasta_to_fastq():
    print('\n[7] vsnp3_fasta_to_fastq: reads per record from its own length; seeded (M9)')
    try:
        import vsnp3_fasta_to_fastq as f2f
    except ImportError as e:
        skip('fasta_to_fastq', f'import failed: {e}')
        return
    import random
    random.seed(7)
    chromosome = ''.join(random.choice('ACGT') for _ in range(3000))
    plasmid = ''.join(random.choice('ACGT') for _ in range(1000))
    with scratch():
        write('two.fasta', f'>chromosome\n{chromosome}\n>plasmid\n{plasmid}\n')
        conv = f2f.Fasta_to_Paired_Fastq('two.fasta', coverage=10, read_length=100)
        with gzip.open(conv.fastq_r1_file, 'rt') as f:
            first = f.read()
        n_reads = first.count('\n') // 4
        check('a 3,000 + 1,000 bp FASTA at 10x with 100 bp reads gives 300 + 100 = 400 pairs (was 800)',
              n_reads == 400, f'{n_reads} reads')
        headers = [line for line in first.splitlines() if line.startswith('@')]
        check('read headers keep the SIM:<n>:FC:1:1:<x>:<y> 1:N:0:<barcode> form',
              headers and headers[0].startswith('@SIM:1:FC:1:1:') and headers[0].endswith(' 1:N:0:AAACCCGGGTTT'), headers[:1])
        conv2 = f2f.Fasta_to_Paired_Fastq('two.fasta', coverage=10, read_length=100)
        with gzip.open(conv2.fastq_r1_file, 'rt') as f:
            second = f.read()
        check('the same FASTA gives identical reads on a second run (seeded)', first == second)
        conv3 = f2f.Fasta_to_Paired_Fastq('two.fasta', coverage=10, read_length=100, seed=1)
        with gzip.open(conv3.fastq_r1_file, 'rt') as f:
            third = f.read()
        check('another seed gives different reads', third != first)
    check('reads_for_length is the whole-genome formula applied per record',
          f2f.reads_for_length(4000, 10, 100) == 400 and f2f.reads_for_length(3, 10, 100) == 0)


# ------------------------------------------------------------ group reporter

def section_group_reporter():
    print('\n[8] vsnp3_group_reporter: the step 2 record rule; no workbook is not a failure (M16, M15)')
    try:
        import vsnp3_group_reporter as gr
    except ImportError as e:
        skip('group_reporter', f'import failed: {e}')
        return
    check('NoDefiningSnpWorkbook is a GroupLookupError', issubclass(gr.NoDefiningSnpWorkbook, gr.GroupLookupError))
    reporter = gr.GroupReporter.__new__(gr.GroupReporter)
    reporter.error = None
    reporter.defining_snps = {'chr1:100': 'grpA', 'chr1:200': 'grpB', 'chr1:300': 'grpC', 'chr1:400': 'grpD'}
    reporter.inverted_defining_snps = {}
    reporter.mixed_defining_snps = {}
    q = QUAL_THRESHOLD + 50
    with scratch():
        vcf = ('##fileformat=VCFv4.2\n'
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n'
               f'chr1\t100\t.\tA\tG\t{q}\t.\tAC=2;MQ={MQ_THRESHOLD}\tGT\t1/1\n'       # MQ exactly at threshold
               f'chr1\t200\t.\tA\tG\t{q}\t.\tAC=2;MQ={MQ_THRESHOLD - 1}\tGT\t1/1\n'   # below
               f'chr1\t300\t.\tA\tGT\t{q}\t.\tAC=2;MQ={MQ_THRESHOLD + 4}\tGT\t1/1\n'  # multi-base ALT
               f'chr1\t400\t.\tA\tG\t{q}\t.\tAC=1;MQ={MQ_THRESHOLD + 4}\tGT\t0/1\n')  # mixed
        write('s_zc.vcf', vcf)
        _, found, mixed = reporter.find_initial_positions('s_zc.vcf')
        check('MQ exactly at the threshold passes (>=, as in step 2; was >)', 'chr1:100' in found, f'{found}')
        check('MQ below the threshold does not pass', 'chr1:200' not in found)
        check('a multi-base ALT does not pass (as in step 2; was accepted)', 'chr1:300' not in found and 'chr1:300' not in mixed)
        check('an AC=1 record goes to the mixed set', 'chr1:400' in mixed and 'chr1:400' not in found)
        groups = reporter.bin_and_html_table('s_zc.vcf', found, mixed)
        check('groups follow: grpA (MQ 56) and grpD (mixed), not grpB or grpC', groups == ['grpA', 'grpD'], f'{groups}')


# --------------------------------------------------------- alignment helpers

def section_alignment():
    print('\n[9] vsnp3_alignment_vcf: threads, bwa batch pin, replacing moves, reference name, contig check (M18, M12, L2, M19)')
    try:
        import vsnp3_alignment_vcf as av
    except ImportError as e:
        skip('alignment_vcf', f'import failed: {e}')
        return
    from vsnp3_file_setup import worker_count
    saved = os.environ.get('VSNP3_MAX_CPUS')
    try:
        os.environ['VSNP3_MAX_CPUS'] = '2'
        check('VSNP3_MAX_CPUS=2 caps the step 1 thread count at 2', av.step1_threads() == 2, f'{av.step1_threads()}')
        os.environ['VSNP3_MAX_CPUS'] = '64'
        check('a large cap leaves the thread count at 8 (the historical value)', av.step1_threads() == 8, f'{av.step1_threads()}')
        os.environ.pop('VSNP3_MAX_CPUS')
        check('unset: min(8, worker_count())', av.step1_threads() == max(1, min(8, worker_count())))
    finally:
        if saved is None:
            os.environ.pop('VSNP3_MAX_CPUS', None)
        else:
            os.environ['VSNP3_MAX_CPUS'] = saved
    check('the bwa batch pin is 10 Mbp x 8 threads, bwa\'s own default at -t 8', av.BWA_BATCH_BASES == 80_000_000)

    with scratch():
        os.makedirs('dest')
        write('a.txt', 'new')
        write(os.path.join('dest', 'a.txt'), 'old')
        av.move_replacing('a.txt', 'dest')
        check('move_replacing replaces an existing file (shutil.move raised)',
              not os.path.exists('a.txt') and open('dest/a.txt').read() == 'new')
        os.makedirs('unmapped_reads')
        write('unmapped_reads/r.fq', 'x')
        os.makedirs('dest/unmapped_reads')
        write('dest/unmapped_reads/stale.fq', 'y')
        av.move_replacing('unmapped_reads', 'dest')
        check('move_replacing replaces an existing directory',
              os.path.exists('dest/unmapped_reads/r.fq') and not os.path.exists('dest/unmapped_reads/stale.fq'))

        check('concatenated reference keeps its usual name when nothing is in the way',
              av.concatenated_reference_name(['/ref/NC_002945v4.fasta']) == 'NC_002945v4.fasta')
        check('two FASTAs join their stems as before',
              av.concatenated_reference_name(['/r/a.fasta', '/r/b.fna']) == 'a_b.fasta')
        write('NC_002945v4.fasta', '>x\nACGT\n')
        distinct = av.concatenated_reference_name(['NC_002945v4.fasta'])
        check('with the user\'s own FASTA in the working directory a distinct name is used',
              distinct == 'NC_002945v4.vsnp3_reference.fasta', distinct)
        import re
        check('...whose stem before the first "." is unchanged, so alignment_<ref> is unchanged',
              re.sub('[.].*', '', distinct) == 'NC_002945v4')

        write('calls.vcf', '##x\n#CHROM\tPOS\n chrA\t1\nchrA\t5\n'.replace(' chrA', 'chrA'))
        missing = av.contigs_without_records(['chrA', 'chrB', 'chrC'], 'calls.vcf')
        check('contigs with mapped reads but no VCF record are named', missing == ['chrB', 'chrC'], f'{missing}')


# ---------------------------------------------------------------- step 1 CLI

def section_step1_cli():
    print('\n[10] vsnp3_step1.py command line: fast --version, help text, -f with -b refused (G5, H11, L14, M14)')
    step1 = os.path.join(BIN, 'vsnp3_step1.py')
    for flag in ('--version', '-v'):
        proc = subprocess.run([sys.executable, '-X', 'importtime', step1, flag], capture_output=True, text=True, timeout=120)
        check(f'{flag} prints exactly what argparse printed',
              proc.returncode == 0 and proc.stdout == f'vsnp3_step1.py: version {__version__}\n', repr(proc.stdout))
        heavy = [m for m in ('pandas', 'Bio', 'dask', 'numpy') if f' {m}\n' in proc.stderr or f'| {m}\n' in proc.stderr or f' {m}.' in proc.stderr]
        check(f'{flag} is answered before pandas, Biopython, numpy or dask are imported', not heavy, f'imported: {heavy}')

    proc = subprocess.run([sys.executable, step1, '--help'], capture_output=True, text=True, timeout=120)
    text = proc.stdout
    check('--help exits 0', proc.returncode == 0, proc.stderr[-300:])
    check(f'--help names the current version in the conda line', f'conda install vsnp3={__version__}' in text)
    check('--help examples give the reference FASTA with -r', '-r ../myreference/*fasta -b ../myreference/*gbk' in text
          and '-r *fasta -b *gbk' in text)
    check('--help no longer shows "-f *fasta -b *gbk"', '-f *fasta -b' not in text and '-f ../myreference' not in text)
    check('-assemble_unmap help says the flag enables assembly',
          'assemble the unmapped reads' in text and 'skip assembly' not in text)

    with scratch():
        write('g.fasta', '>x\nACGT\n')
        write('g.gbk', 'LOCUS x\n')
        proc = subprocess.run([sys.executable, step1, '-f', 'g.fasta', '-b', 'g.gbk'], capture_output=True, text=True, timeout=120)
        check('-f with -b and no -r/-t is refused by argparse (exit 2) with the reason',
              proc.returncode == 2 and 'not the reference' in proc.stderr, f'exit {proc.returncode}: {proc.stderr[-300:]}')


# ----------------------------------------------------------------- bruc mlst

def section_bruc_mlst():
    print('\n[11] vsnp3_bruc_mlst: missing picard is a clear message; an unreadable VCF raises (M21)')
    try:
        import vsnp3_bruc_mlst as bm
    except ImportError as e:
        skip('bruc_mlst', f'import failed: {e}')
        return
    check('missing_tools names an absent program', bm.missing_tools(('no-such-tool-xyz', 'python3' if shutil.which('python3') else 'sh')) == ['no-such-tool-xyz'])
    real_which = bm.shutil.which
    bm.shutil.which = lambda name, *a, **k: None if name == 'picard' else '/usr/bin/' + name
    try:
        with scratch():
            try:
                bm.Bruc_MLST('x_R1.fastq.gz', None)
                check('Bruc_MLST without picard raises RuntimeError naming picard', False, 'no exception')
            except RuntimeError as e:
                check('Bruc_MLST without picard raises RuntimeError naming picard',
                      'picard' in str(e) and 'install' in str(e), str(e))
            check('...before writing ST1-MLST.fasta', not os.path.exists('ST1-MLST.fasta'))
    finally:
        bm.shutil.which = real_which
    obj = bm.Bruc_MLST.__new__(bm.Bruc_MLST)
    with scratch():
        try:
            obj.read_vcf('does_not_exist.vcf')
            check('read_vcf on a missing file raises RuntimeError (was an empty frame, KeyError later)', False, 'no exception')
        except RuntimeError as e:
            check('read_vcf on a missing file raises RuntimeError (was an empty frame, KeyError later)', 'could not read' in str(e), str(e))
        write('empty.vcf', '')
        try:
            obj.read_vcf('empty.vcf')
            check('read_vcf on a zero-byte file raises RuntimeError', False, 'no exception')
        except RuntimeError as e:
            check('read_vcf on a zero-byte file raises RuntimeError', 'no #CHROM' in str(e), str(e))
        write('header_only.vcf', '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n')
        frame = obj.read_vcf('header_only.vcf')
        check('a header-only VCF (no variant) is an empty frame with columns, a valid result',
              frame.empty and 'POS' in frame.columns and 'ALT' in frame.columns)


# ------------------------------------------------------------ kraken2 wrapper

def section_kraken2():
    print('\n[12] vsnp3_kraken2_wrapper: threads never below 1 (M18)')
    try:
        import vsnp3_kraken2_wrapper as kw
    except ImportError as e:
        skip('kraken2_wrapper', f'import failed: {e}')
        return
    real = kw.multiprocessing.cpu_count
    kw.multiprocessing.cpu_count = lambda: 1
    try:
        k = kw.Kraken2_Identification(FASTQ_R1='x_R1.fastq.gz', database='db')
    finally:
        kw.multiprocessing.cpu_count = real
    check('a one-core machine gives --threads 1, not -1', k.threads == 1, f'{k.threads}')


# ------------------------------------------------------------ download GCA

def section_download():
    print('\n[13] vsnp3_download_GCA_fasta_get_metadata: exhausted retries are reported (L5)')
    try:
        import vsnp3_download_GCA_fasta_get_metadata as dl
    except ImportError as e:
        skip('download_GCA', f'import failed: {e}')
        return
    import contextlib
    real_system, real_sleep = dl.os.system, dl.time.sleep
    dl.os.system = lambda cmd: 1          # the datasets tool "runs" and leaves no zip
    dl.time.sleep = lambda s: None
    try:
        with scratch(), contextlib.redirect_stdout(io.StringIO()):    # the retry countdown
            try:
                dl.Downloader(gca_accession='GCA_000000000')
                outcome = (False, 'no exception')
            except RuntimeError as e:
                outcome = ('download failed after 4 attempts' in str(e), str(e))
        check('four failed downloads raise RuntimeError instead of falling through', *outcome)
    finally:
        dl.os.system, dl.time.sleep = real_system, real_sleep


# -------------------------------------------------------- validator, longread

def section_encodings():
    print('\n[14] non-ASCII bytes in text files are tolerated (H6)')
    try:
        import vsnp3_input_validator as iv
        import vsnp3_longread as lr
    except ImportError as e:
        skip('encodings', f'import failed: {e}')
        return
    with scratch():
        write('reads.fastq', '@read-é\nACGT\n+\nIIII\n')
        ok, msg = iv.InputValidator().validate_fastq_format('reads.fastq')
        check('a plain FASTQ with an accented read name validates', ok, msg)
        write('ref.fasta', '>contig description Müller\nACGT\n')
        ok, msg = iv.InputValidator().validate_fasta_format('ref.fasta')
        check('a FASTA with an accented description validates', ok, msg)
        seqs = lr.read_fasta('ref.fasta')
        check('vsnp3_longread.read_fasta reads it', seqs == {'contig': 'ACGT'}, f'{seqs}')
        # a Latin-1 byte, which is not valid UTF-8, must not be fatal either
        with open('latin.fasta', 'wb') as f:
            f.write(b'>contig M\xfcller\nACGT\n')
        ok, msg = iv.InputValidator().validate_fasta_format('latin.fasta')
        check('a FASTA with a Latin-1 byte validates (errors="replace")', ok, msg)
        check('vsnp3_longread.read_fasta reads it too', lr.read_fasta('latin.fasta') == {'contig': 'ACGT'})


# --------------------------------------------------------------- vcf_search

def section_vcf_search():
    print('\n[15] vcf_search.sh: portable splitting (L5)')
    script = os.path.join(BIN, 'vcf_search.sh')
    if not shutil.which('bash'):
        skip('vcf_search.sh', 'bash not found')
        return
    with scratch():
        head = '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        write('both.vcf', head + 'chr1\t100\t.\tA\tG\t200\t.\tAC=2\nchr1\t200\t.\tC\tT\t200\t.\tAC=2\n')
        write('first_only.vcf', head + 'chr1\t100\t.\tA\tG\t200\t.\tAC=2\nchr1\t1000\t.\tC\tT\t200\t.\tAC=2\n')
        write('neither.vcf', head + 'chr1\t300\t.\tA\tG\t200\t.\tAC=2\n')
        proc = subprocess.run(['bash', script, 'chr1:100, chr1:200!'], capture_output=True, text=True, timeout=60)
        check('required chr1:100 and excluded chr1:200! select first_only.vcf alone',
              proc.returncode == 0 and proc.stdout.split() == ['first_only.vcf'], f'exit {proc.returncode}: {proc.stdout!r} {proc.stderr!r}')
        proc = subprocess.run(['bash', script, 'chr1:100,chr1:200'], capture_output=True, text=True, timeout=60)
        check('two required positions without a space select both.vcf alone',
              proc.stdout.split() == ['both.vcf'], f'{proc.stdout!r} {proc.stderr!r}')


def main():
    print('=' * 78)
    print('vSNP3 step 1 unit tests')
    print('=' * 78)
    section_krona()
    section_seqkit()
    section_zero_coverage()
    section_spoligotype()
    section_run()
    section_assembly()
    section_fasta_to_fastq()
    section_group_reporter()
    section_alignment()
    section_step1_cli()
    section_bruc_mlst()
    section_kraken2()
    section_download()
    section_encodings()
    section_vcf_search()
    print('\n' + '=' * 78)
    print(f'PASS {len(PASS)}   FAIL {len(FAIL)}')
    if FAIL:
        for name in FAIL:
            print(f'  - {name}')
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
