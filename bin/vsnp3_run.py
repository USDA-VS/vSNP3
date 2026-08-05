#!/usr/bin/env python3

'''
Run external tools so that a failure is reported instead of inferred.

vsnp3_alignment_vcf.py called bwa, samtools, bcftools and freebayes through
os.system with f-string interpolated paths.  Nothing checked a return code and
most calls sent stderr to /dev/null, so a failed aligner produced an empty SAM,
an empty BAM, an empty VCF, and a run that looked successful at 0x coverage.
The paths were also unquoted, so a working directory containing a space broke the
command line silently.

Two rules here:

  - argv is a list.  Quoting is then structural rather than textual, which is what
    actually makes a path with a space safe.  Adding quote characters to a shell
    string does not, because the next interpolation reintroduces the problem.
  - a tool that exits 0 having written nothing useful is still a failure, so
    require_output() checks the post-condition separately from the return code.
'''

import os
import shlex
import subprocess

from vsnp3_version import __version__          # noqa: F401  re-exported for callers


class ExternalToolError(RuntimeError):
    '''An external tool failed, or produced no usable output.'''

    def __init__(self, argv, returncode=None, stderr=None, detail=None):
        self.argv = list(argv)
        self.returncode = returncode
        self.stderr = stderr or ''
        printable = ' '.join(shlex.quote(str(a)) for a in self.argv)
        message = f'{self.argv[0]} failed'
        if returncode is not None:
            message += f' (exit {returncode})'
        if detail:
            message += f': {detail}'
        message += f'\n  command: {printable}'
        tail = self.stderr.strip().splitlines()[-15:]
        if tail:
            message += '\n  stderr:\n' + '\n'.join(f'    {line}' for line in tail)
        super().__init__(message)


def run(argv, cwd=None, stdout_path=None, check=True, log=None):
    '''
    Run argv as a list.  stdout_path replaces a `> file` shell redirect without
    needing shell=True.  Raises ExternalToolError on a non-zero exit when check.
    '''
    argv = [str(a) for a in argv]
    if log is not None:
        log.append(' '.join(shlex.quote(a) for a in argv))
    stdout = open(stdout_path, 'wb') if stdout_path else subprocess.PIPE
    try:
        proc = subprocess.run(argv, cwd=cwd, stdout=stdout,
                              stderr=subprocess.PIPE, check=False)
    except FileNotFoundError as e:
        raise ExternalToolError(argv, detail=f'not found on PATH ({e.strerror})') from e
    finally:
        if stdout_path:
            stdout.close()
    if check and proc.returncode != 0:
        raise ExternalToolError(argv, proc.returncode,
                                _text(proc.stderr))
    return proc


def run_pipe(argv_a, argv_b, cwd=None, stdout_path=None, check=True, log=None):
    '''
    argv_a | argv_b, for the two bcftools mpileup-into-call pipelines.  A shell
    pipeline reports only the last exit status, so both are checked here.
    '''
    argv_a = [str(a) for a in argv_a]
    argv_b = [str(a) for a in argv_b]
    if log is not None:
        log.append(' '.join(shlex.quote(a) for a in argv_a) + ' | '
                   + ' '.join(shlex.quote(a) for a in argv_b))
    stdout = open(stdout_path, 'wb') if stdout_path else subprocess.PIPE
    try:
        first = subprocess.Popen(argv_a, cwd=cwd, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
        second = subprocess.Popen(argv_b, cwd=cwd, stdin=first.stdout,
                                  stdout=stdout, stderr=subprocess.PIPE)
        first.stdout.close()          # let the first process see a broken pipe
        err_b = second.communicate()[1]
        err_a = first.stderr.read()
        first.stderr.close()
        first.wait()
    except FileNotFoundError as e:
        raise ExternalToolError(argv_a + ['|'] + argv_b,
                                detail=f'not found on PATH ({e.strerror})') from e
    finally:
        if stdout_path:
            stdout.close()
    if check and first.returncode != 0:
        raise ExternalToolError(argv_a, first.returncode, _text(err_a))
    if check and second.returncode != 0:
        raise ExternalToolError(argv_b, second.returncode, _text(err_b))
    return first.returncode, second.returncode


def require_output(path, min_bytes=1, bam=False, vcf_records=False, what=None):
    '''
    Assert a stage actually produced something usable.

    Checked separately from the return code because several of these tools exit 0
    after writing an empty or truncated file -- which is exactly how a failed
    alignment used to become a complete-looking sample at 0x coverage.
    '''
    label = what or os.path.basename(path)
    if not os.path.exists(path):
        raise ExternalToolError(['require_output', path],
                                detail=f'{label} was not created')
    size = os.path.getsize(path)
    if size < min_bytes:
        raise ExternalToolError(['require_output', path],
                                detail=f'{label} is {size} bytes, expected at least {min_bytes}')
    if bam:
        proc = subprocess.run(['samtools', 'quickcheck', '-v', path],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              check=False)
        if proc.returncode != 0:
            raise ExternalToolError(['samtools', 'quickcheck', path],
                                    proc.returncode, _text(proc.stderr),
                                    detail=f'{label} is not a valid BAM')
    if vcf_records:
        with open(path, 'r', errors='replace') as f:
            for line in f:
                if not line.startswith('#') and line.strip():
                    return
        raise ExternalToolError(['require_output', path],
                                detail=f'{label} contains no variant records')


def mapped_read_count(bam):
    '''Mapped reads in a BAM. Zero means the alignment produced nothing.'''
    proc = run(['samtools', 'view', '-c', '-F', '4', bam])
    try:
        return int(_text(proc.stdout).strip())
    except (AttributeError, ValueError) as e:
        raise ExternalToolError(['samtools', 'view', '-c', '-F', '4', bam],
                                detail=f'could not read a count from samtools ({e})') from e


def _text(raw):
    if raw is None:
        return ''
    if isinstance(raw, bytes):
        return raw.decode('utf-8', errors='replace')
    return raw
