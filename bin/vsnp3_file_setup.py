#!/usr/bin/env python

from vsnp3_version import __version__
from vsnp3_version import FASTQ_USABILITY, REFERENCE_USABILITY

import os
import shutil
import re
import time
import locale
import logging
import pandas as pd
import multiprocessing
from datetime import datetime
import yaml

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class bcolors:
    PURPLE = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE='\033[37m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    ENDC = '\033[0m'

class Setup:
    ''' 
    Standarize setup
    '''

    def __init__(self, SAMPLE_NAME=None, FASTA=None, FASTQ_R1=None, FASTQ_R2=None, reference=None, gbk=None, debug=False):
        self.cwd = os.getcwd()
        self._refuse_whitespace_paths(FASTA, FASTQ_R1, FASTQ_R2, reference, gbk)
        self.FASTA = FASTA
        self.FASTQ_R1 = FASTQ_R1
        self.FASTQ_R2 = FASTQ_R2
        self.reference = reference
        self.gbk = gbk
        self.debug = debug
        self.remove_copied_gbk = False
        if FASTQ_R2:
            self.paired = True
        else:
            self.paired = False
        if FASTA:
            try: #IF FASTA provided as path variable copy local
                shutil.copy(FASTA, self.cwd)
            except shutil.SameFileError:
                pass
            FASTA = os.path.basename(FASTA)
            self.FASTA = os.path.join(self.cwd, FASTA)
            self.sample_name = re.sub('[_.].*', '', FASTA)
            self.fasta_name = re.sub('[.].*', '', FASTA) #explict FASTA name if needed
        if gbk:
            try:
                updated_gbk_list=[]
                for each in gbk:
                    shutil.copy(each, self.cwd)
                    updated_gbk_list.append(os.path.join(self.cwd, os.path.basename(each)))
                self.gbk = updated_gbk_list
                self.remove_copied_gbk = True
            except shutil.SameFileError:
                self.gbk = gbk
                self.remove_copied_gbk = False
        if reference:
            self.reference_name = re.sub('[.].*', '', reference)
            self.reference = reference
            try: #IF FASTA provided as path variable copy local
                shutil.copy(reference, self.cwd)
            except shutil.SameFileError:
                pass
            reference = os.path.basename(reference)
            self.reference = os.path.join(self.cwd, reference)
        elif FASTA: # if no explict reference make the FASTA a "reference"
            self.reference_name = re.sub('[.].*', '', FASTA)
            self.reference = self.FASTA
        self.FASTQ_list=[]
        self.FASTQ_dict={}
        try:
            try: #IF FASTQ provided as path variable copy local
                shutil.copy(FASTQ_R1, self.cwd)
            except shutil.SameFileError:
                pass
            FASTQ_R1 = os.path.basename(FASTQ_R1)
            self.FASTQ_R1 = os.path.join(self.cwd, FASTQ_R1)
            self.sample_name = re.sub('[_.].*', '', FASTQ_R1) # default to FASTQ_R1 as sample name if also FASTA
            self.fastq_name = re.sub('[_.].*', '', FASTQ_R1) #explict FASTQ name if needed
            self.FASTQ_list.append(self.FASTQ_R1)
            self.FASTQ_dict = {'FASTQ_R1': self.FASTQ_R1}
        except TypeError:
            self.FASTQ_R1 = None
        if FASTQ_R2:
            try: #IF FASTQ provided as path variable copy local
                shutil.copy(FASTQ_R2, self.cwd)
            except shutil.SameFileError:
                pass
            FASTQ_R2 = os.path.basename(FASTQ_R2)
            self.FASTQ_R2 = os.path.join(self.cwd, FASTQ_R2)
            self.FASTQ_list.append(self.FASTQ_R2)
            self.FASTQ_dict['FASTQ_R2'] = self.FASTQ_R2
        else:
            self.FASTQ_R2 = None

        self.startTime = datetime.now()
        # At least one CPU. cpu_count() - 2 goes to zero or negative on a 1-2 core
        # machine, and a worker count of zero is not a usable value downstream.
        self.cpus = max(1, multiprocessing.cpu_count() - 2)
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        if SAMPLE_NAME:
            self.sample_name = SAMPLE_NAME

    def _refuse_whitespace_paths(self, *paths):
        '''
        Refuse up front if any path, or the working directory, contains whitespace.

        vsnp3_alignment_vcf.py builds its bwa, samtools, bcftools and freebayes
        command lines by interpolating these paths into a shell string, unquoted.
        A single space splits one argument into two, so the aligner reads a file
        that does not exist, writes an empty SAM, and every later stage produces an
        empty-but-present output: the run completes, reports success, and shows 0x
        coverage.  Nothing downstream can distinguish that from a genuinely
        unalignable sample.

        Checked here rather than fixed at each of the 15 call sites because this
        catches the whole class at the one point every path passes through.
        '''
        offenders = []
        if any(ch.isspace() for ch in self.cwd):
            offenders.append(f'working directory: {self.cwd}')
        for path in paths:
            for each in (path if isinstance(path, (list, tuple)) else [path]):
                if each and any(ch.isspace() for ch in str(each)):
                    offenders.append(str(each))
        if offenders:
            raise ValueError(
                'vSNP3 cannot run with whitespace in a file or directory path, '
                'because the external aligner and variant caller command lines '
                'would break silently and produce an empty result that looks like '
                'a successful run at 0x coverage.\n  '
                + '\n  '.join(offenders)
                + '\nRename or move so no path contains a space, tab or newline.')

    def print_time(self,):
        '''
        description
        '''
        self.run_time = datetime.now() - self.startTime
        runtime_str = str(datetime.now() - self.startTime)
        print('\n\nruntime: {}\n'.format(runtime_str))

    def print_run_time(self, tool):
        print('{}\n{} {}'.format(bcolors.RED, tool, bcolors.ENDC))
        now = datetime.now()
        print('{}{}{}'.format(bcolors.WHITE, now.strftime("%Y-%m-%d %H:%M:%S"), bcolors.ENDC))


def format_elapsed(seconds):
    '''"2 hours 5 minutes 9 seconds", with zero hours and minutes left out.'''
    total = int(round(seconds))
    hours, rest = divmod(total, 3600)
    minutes, secs = divmod(rest, 60)
    parts = []
    if hours:
        parts.append(f'{hours} hour{"s" if hours != 1 else ""}')
    if minutes:
        parts.append(f'{minutes} minute{"s" if minutes != 1 else ""}')
    parts.append(f'{secs} second{"s" if secs != 1 else ""}')
    return ' '.join(parts)


class PhaseTimer:
    '''
    Wall-clock time per named phase of a run, in the order the phases ran.

    Monotonic, so a clock adjustment during a 13-hour run cannot move the
    numbers.  started_at is the wall-clock start for the report; elapsed is
    read at the moment the report is written, not before, so it covers
    everything that ran.
    '''

    def __init__(self):
        self.started_at = datetime.now()
        self._start = time.monotonic()
        self._last = self._start
        self.phases = []            # [(name, seconds)]

    def mark(self, name):
        '''Close the phase that began at the previous mark (or at the start).'''
        now = time.monotonic()
        self.phases.append((name, now - self._last))
        self._last = now

    def absorb(self, phases):
        '''Take over phases timed elsewhere, then continue from now.'''
        self.phases.extend(phases)
        self._last = time.monotonic()

    @property
    def elapsed(self):
        return time.monotonic() - self._start

    def print(self):
        width = max((len(name) for name, _ in self.phases), default=10)
        print('\nTime by phase:')
        for name, seconds in self.phases:
            print(f'  {name:<{width}}  {format_elapsed(seconds)}')
        print(f'  {"total":<{width}}  {format_elapsed(self.elapsed)}')

# Sections, in report order.  (title, [(label, excel_dict key), ...])
# Keys absent from excel_dict are skipped, and a section with no present keys is
# omitted entirely -- so a run without spoligotyping or assembly simply has no
# such section rather than a table of blanks.
SECTIONS = [
    ('Reference and alignment', [
        ('Reference', 'Reference'),
        ('Reference length', 'Reference Length'),
        ('Aligner', 'Aligner'),
        ('Sequencing platform', 'Platform'),
        ('BAM file', 'BAM File'),
        ('BAM / reference', 'BAM/Reference File'),
        ('Mapped paired reads', 'Mapped Paired Reads'),
        ('Mapped single reads', 'Mapped Single Reads'),
        ('Duplicate paired reads', 'Duplicate Paired Reads'),
        ('Duplicate single reads', 'Duplicate Single Reads'),
        ('Duplicates, percent of mapped', 'Duplicate Percent of Mapped Reads'),
        ('Unmapped reads', 'Unmapped Reads'),
        ('Unmapped percent', 'Unmapped Percent'),
    ]),
    ('Coverage', [
        ('Average depth', 'Average Depth'),
        ('Genome with coverage', 'Genome with Coverage'),
        ('Bases with no coverage', 'No Coverage Bases'),
        ('Percent of reference with zero coverage', 'Percent Ref with Zero Coverage'),
    ]),
    ('SNPs', [
        ('SNP caller', 'SNP Caller'),
        ('Quality SNPs', 'Quality SNPs'),
        ('Ambiguous SNPs', 'Ambiguous SNPs'),
        ('Defining SNP groups', 'Groups'),
        # Long-read (nanopore) samples only; absent keys are not shown.
        ('Long-read error rate', 'Long-read Error Rate'),
        ('Long-read evidence QUAL', 'Long-read QUAL Model'),
        ('Long-read calls limited to N', 'Long-read Calls Limited to N'),
    ]),
    ('Spoligotype', [
        ('SB number', 'Spoligotype SB Number'),
        ('Octal code', 'Spoligotype Octal Code'),
        ('Binary code', 'Spoligotype Binary Code'),
        ('Spacer counts', 'Spoligotype Spacer Counts'),
    ]),
    ('Assembly of unmapped reads', [
        ('Contigs', 'Contig count'),
        ('Total length', 'Total length'),
        ('Longest contig', 'Longest contig'),
        ('N50', 'N50'),
        ('Contig length distribution', 'Contig length counts <|301-999bp|>'),
        ('Unmapped assembled contigs', 'Unmapped Assembled Contigs'),
    ]),
]

# Read-level metrics, rendered as one row per metric with an R1 and an R2 column.
READ_METRICS = [
    ('File', 'FASTQ_R1', 'FASTQ_R2'),
    ('File size', 'R1 File Size', 'R2 File Size'),
    ('Read count', 'R1 Read Count', 'R2 Read Count'),
    ('Total bases', 'R1 Length Sum', 'R2 Length Sum'),
    ('Minimum length', 'R1 Min Length', 'R2 Min Length'),
    ('Average length', 'R1 Ave Length', 'R2 Ave Length'),
    ('Maximum length', 'R1 Max Length', 'R2 Max Length'),
    ('Passing Q20', 'R1 Passing Q20', 'R2 Passing Q20'),
    ('Passing Q30', 'R1 Passing Q30', 'R2 Passing Q30'),
    ('Mean read quality', 'R1 Read Quality Ave', 'R2 Read Quality Ave'),
]

# Keys shown in the summary tiles, and therefore not repeated in the tables below.
TILE_KEYS = ('Average Depth', 'Genome with Coverage', 'Quality SNPs')

VERDICT_CLASS = {'acceptable': 'good', 'questionable': 'warn', 'poor': 'bad'}


# An absolute POSIX path: a leading slash, at least one more slash, and no
# whitespace or quotes.  The lookbehind keeps URLs intact -- the "//" in
# "http://www.htslib.org/" is rejected because the first slash follows a colon
# and the second follows a slash, and tool version banners do print URLs.
_ABSOLUTE_PATH = re.compile(r'(?<![\w:/])/(?:[^\s/"\']+/)+[^\s/"\']*')


def redact_paths(value):
    '''
    Replace absolute paths in report text with just the filename.

    The run log records each command line verbatim, and those carry the full
    path of the reference and of the reads -- which on a shared system names the
    user, the project and the directory layout.  A report is meant to be
    sendable to whoever needs the result, so it should not have to be read for
    disclosures first.  The basename is what the rest of the report already
    shows for a file, and it leaves the command lines readable.

    Applied at the report boundary rather than at the ~20 places that build
    these strings: that is one place to audit, and it also covers paths that
    arrive inside an exception message, which no call site controls.
    '''
    if not isinstance(value, str):
        return value

    def shorten(match):
        path = match.group(0)
        if path.startswith('/dev/'):
            return path                  # /dev/null discloses nothing, and
                                         # rewriting it makes a command look wrong
        return os.path.basename(path.rstrip('/')) or path

    return _ABSOLUTE_PATH.sub(shorten, value)


# Rendered lazily via _template(): see the note on the import below.
_TEMPLATE_SOURCE = '''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ sample }} &middot; vSNP3 step 1 report</title>
<style>
  /*
    One light palette, deliberately not following prefers-color-scheme.  A
    dark-mode browser rendered this report on a near-black background, which is
    wrong for a document that is read next to its own printout and next to the
    PDF rendered from it -- the PDF is always light, because no print pipeline
    honours a dark scheme.
  */
  :root {
    --bg:#f6f7f9; --card:#ffffff; --ink:#1a1d21; --muted:#5b6572;
    --line:#dfe3e8; --accent:#1f5c8b; --accent-soft:#eaf1f7;
    --good:#1a7f4b; --good-bg:#e8f5ee;
    --warn:#8a5a00; --warn-bg:#fdf3e0;
    --bad:#a52121;  --bad-bg:#fbeaea;
  }
  * { box-sizing:border-box; }
  body {
    margin:0; padding:2rem 1.25rem 4rem; background:var(--bg); color:var(--ink);
    font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  }
  .wrap { max-width:1060px; margin:0 auto; }
  header {
    display:flex; flex-wrap:wrap; align-items:center; gap:1.25rem;
    padding-bottom:1.25rem; border-bottom:3px solid var(--accent); margin-bottom:1.75rem;
  }
  header img { max-height:56px; max-width:220px; }
  .title { flex:1 1 260px; }
  .title h1 { margin:0; font-size:1.6rem; letter-spacing:-.01em; }
  .title .sub { color:var(--muted); font-size:.9rem; margin-top:.2rem; }
  .stamp { text-align:right; color:var(--muted); font-size:.85rem; }
  @media (max-width:640px) { .stamp { text-align:left; } }

  /*
    Flex rather than grid: WeasyPrint renders the PDF from this same file and
    does not implement CSS grid before 63, where a grid container collapses to
    one tile per row.  Flex degrades to the same layout everywhere.
  */
  .tiles { display:flex; flex-wrap:wrap; gap:.85rem; margin-bottom:1.75rem; }
  .tile { flex:1 1 170px; min-width:170px; background:var(--card);
          border:1px solid var(--line); border-radius:10px; padding:.9rem 1rem; }
  .tile .k { color:var(--muted); font-size:.74rem; text-transform:uppercase; letter-spacing:.06em; }
  .tile .v { font-size:1.45rem; font-weight:600; margin-top:.3rem; word-break:break-word; }
  .tile .n { color:var(--muted); font-size:.74rem; margin-top:.25rem; }
  .tile.good { background:var(--good-bg); border-color:var(--good); }
  .tile.warn { background:var(--warn-bg); border-color:var(--warn); }
  .tile.bad  { background:var(--bad-bg);  border-color:var(--bad); }
  .tile.good .v { color:var(--good); } .tile.warn .v { color:var(--warn); } .tile.bad .v { color:var(--bad); }

  section { background:var(--card); border:1px solid var(--line); border-radius:10px; margin-bottom:1.1rem; overflow:hidden; }
  section > h2 {
    margin:0; padding:.7rem 1rem; font-size:.86rem; font-weight:600;
    text-transform:uppercase; letter-spacing:.07em; color:var(--accent);
    background:var(--accent-soft); border-bottom:1px solid var(--line);
  }
  .tbl { width:100%; border-collapse:collapse; }
  .tbl th, .tbl td { padding:.5rem 1rem; text-align:left; border-bottom:1px solid var(--line); vertical-align:top; }
  .tbl tr:last-child th, .tbl tr:last-child td { border-bottom:none; }
  .tbl th { font-weight:500; color:var(--muted); width:40%; }
  .tbl td { font-variant-numeric:tabular-nums; word-break:break-word; }
  .tbl thead th { color:var(--ink); font-weight:600; width:auto; }
  .mono { font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace; font-size:.86em; }
  .scroll { overflow-x:auto; }

  .note { background:var(--warn-bg); border:1px solid var(--warn); color:var(--warn);
          border-radius:8px; padding:.7rem 1rem; margin-bottom:1.1rem; font-size:.9rem; }
  footer { margin-top:2rem; padding-top:1rem; border-top:1px solid var(--line);
           color:var(--muted); font-size:.82rem; }
  footer .row { display:flex; flex-wrap:wrap; gap:.4rem 1.5rem; }
  /*
    The run log is a plain section rather than a <details>/<summary> disclosure.
    WeasyPrint gave a hidden-summary <details> the full remaining page height, so
    the PDF carried a page-tall empty grey box, and a widget that only opens on a
    click is of no use in a printed report anyway.
  */
  pre { background:var(--bg); border:1px solid var(--line); border-radius:6px;
        padding:.7rem; overflow-x:auto; font-size:.8rem; margin:0; }
  .log { padding:.8rem 1rem; }

  /*
    Applies to the browser's Print and to the WeasyPrint-rendered PDF, which
    selects print media.  Tiles become inline-blocks here rather than staying
    flex items: inline-block is the one layout every print engine gets right.
  */
  @page { size:letter; margin:0.6in 0.55in; }
  @media print {
    body { background:#fff; padding:0; font-size:10pt; }
    .wrap { max-width:none; }
    /*
      Rows and tiles stay whole; sections are allowed to split.  Holding a whole
      section together pushed any section that would not fit onto the next page
      and left a third of the preceding page blank.
    */
    .tile, tr { break-inside:avoid; }
    section > h2 { break-after:avoid; }
    section, .tile { border-color:#bbb; }
    header { border-bottom-color:#333; }
    .tiles { display:block; margin-bottom:1rem; }
    .tile { display:inline-block; width:31.5%; margin:0 1% .55rem 0;
            vertical-align:top; min-width:0; }
    .scroll { overflow-x:visible; }
    pre { white-space:pre-wrap; word-break:break-word; }
    a { color:inherit; text-decoration:none; }
  }
</style>
</head>
<body>
<div class="wrap">

  <header>
    {% if logo_data_uri %}<img src="{{ logo_data_uri }}" alt="">{% endif %}
    <div class="title">
      <h1>{{ sample }}</h1>
      <div class="sub">{{ description or 'vSNP3 step 1 — alignment and SNP calling' }}</div>
    </div>
    <div class="stamp">
      {{ generated }}<br>vSNP3 {{ version }}
    </div>
  </header>

  {% if tiles %}
  <div class="tiles">
    {% for t in tiles %}
    <div class="tile {{ t.cls }}">
      <div class="k">{{ t.label }}</div>
      <div class="v">{{ t.value }}</div>
      {% if t.note %}<div class="n">{{ t.note }}</div>{% endif %}
    </div>
    {% endfor %}
  </div>
  {% endif %}

  {% for message in notices %}
  <div class="note">{{ message }}</div>
  {% endfor %}

  {% if read_rows %}
  <section>
    <h2>Sequencing quality</h2>
    <div class="scroll">
    <table class="tbl">
      <thead><tr><th>Metric</th><th>Read 1</th><th>{{ 'Read 2' if paired else 'Read 2 (none)' }}</th></tr></thead>
      <tbody>
        {% for label, r1, r2 in read_rows %}
        <tr><th>{{ label }}</th>
            <td class="{{ 'mono' if label == 'File' }}">{{ r1 }}</td>
            <td class="{{ 'mono' if label == 'File' }}">{{ r2 }}</td></tr>
        {% endfor %}
      </tbody>
    </table>
    </div>
  </section>
  {% endif %}

  {% for title, rows in sections %}
  <section>
    <h2>{{ title }}</h2>
    <div class="scroll">
    <table class="tbl">
      <tbody>
        {% for label, value, mono in rows %}
        <tr><th>{{ label }}</th><td class="{{ 'mono' if mono }}">{{ value }}</td></tr>
        {% endfor %}
      </tbody>
    </table>
    </div>
  </section>
  {% endfor %}

  {% for t in tables %}
  <section>
    <h2>{{ t.title }}</h2>
    <div class="scroll">
    <table class="tbl">
      <thead><tr>{% for h in t.headers %}<th>{{ h }}</th>{% endfor %}</tr></thead>
      <tbody>
        {% for row in t.rows %}
        <tr>{% for cell in row %}<td>{{ cell }}</td>{% endfor %}</tr>
        {% endfor %}
      </tbody>
    </table>
    </div>
  </section>
  {% endfor %}

  {% if other_rows %}
  <section>
    <h2>Other metrics</h2>
    <div class="scroll">
    <table class="tbl"><tbody>
      {% for label, value in other_rows %}
      <tr><th>{{ label }}</th><td>{{ value }}</td></tr>
      {% endfor %}
    </tbody></table>
    </div>
  </section>
  {% endif %}

  {% if program_versions %}
  <section>
    <h2>Program versions and run log</h2>
    <div class="log"><pre>{{ program_versions }}</pre></div>
  </section>
  {% endif %}

  <footer>
    <div class="row">
      <span><strong>Sample:</strong> {{ sample }}</span>
      <span><strong>Run:</strong> {{ run_stamp }}</span>
      <span><strong>vSNP3:</strong> {{ version }}</span>
      {% if runtime %}<span><strong>Runtime:</strong> {{ runtime }}</span>{% endif %}
    </div>
  </footer>

</div>
</body>
</html>
'''


_TEMPLATE = None


def _template():
    '''
    Compile the report template on first use.

    jinja2 is imported here rather than at module scope because ten modules
    import vsnp3_file_setup, among them vsnp3_step2.py, which has no need of a
    step 1 report.  A module-scope import would make jinja2 a hard import-time
    requirement of all of them and turn a missing optional dependency into a
    step 2 failure.
    '''
    global _TEMPLATE
    if _TEMPLATE is None:
        import jinja2
        _TEMPLATE = jinja2.Template(_TEMPLATE_SOURCE, autoescape=True)
    return _TEMPLATE


class HtmlReport:
    '''
    Collects the step 1 report and writes it as HTML, then as a PDF rendered
    from that same HTML.

    Everything in the report comes from the same stats dict that produces the
    workbook, so a metric only has to be recorded once.  The PDF is rendered
    from the finished HTML rather than composed separately, so the two cannot
    disagree -- which is what went wrong with the LaTeX report this replaced:
    it was a second, independently maintained description of the same numbers.
    '''

    def __init__(self, sample_name, config_file=None, report_description=None):
        self.sample_name = sample_name
        self.description = report_description
        self.date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.html_file = f'{sample_name}_{self.date_stamp}_report.html'
        self.pdf_file = f'{sample_name}_{self.date_stamp}_report.pdf'
        self.notices = []
        self.tables = []
        self.program_versions = ''
        self.runtime = ''
        self.config = self._load_config(config_file)

    # ------------------------------------------------------------------ config

    @staticmethod
    def _load_config(config_file):
        if config_file is None:
            here = os.path.dirname(os.path.realpath(__file__))
            candidate = os.path.join(here, '..', 'dependencies', 'report_config.yaml')
            config_file = os.path.abspath(candidate)
        config = {'logo_path': None}
        if config_file and os.path.exists(config_file):
            try:
                with open(config_file) as f:
                    loaded = yaml.safe_load(f) or {}
                config.update(loaded)
            except Exception as e:                       # noqa: BLE001
                print(f'Note: could not read report config {config_file}: {e}')
        return config

    def _logo_data_uri(self):
        '''Embed the logo so the report is a single portable file.'''
        import base64
        import mimetypes
        path = self.config.get('logo_path')
        if not path or not os.path.exists(path):
            return None
        mime = mimetypes.guess_type(path)[0] or 'image/png'
        try:
            with open(path, 'rb') as f:
                return f'data:{mime};base64,' + base64.b64encode(f.read()).decode()
        except OSError as e:
            print(f'Note: could not embed logo {path}: {e}')
            return None

    # ------------------------------------------------------------------- build

    def add_notice(self, message):
        '''Surface something the reader needs to know, e.g. a failed group lookup.'''
        if message:
            self.notices.append(message)

    def add_table(self, title, headers, rows):
        '''
        Attach a multi-row table to the report.

        For results that are a list rather than one value per metric, and so
        have no representation in the stats dict -- the sourmash similarity
        ranking is the one case.  It appeared only in the old LaTeX report,
        which is how it came to be missing from the HTML one.
        '''
        rows = [list(r) for r in (rows or [])]
        if rows:
            self.tables.append({'title': title,
                                'headers': list(headers),
                                'rows': rows})

    def _tiles(self, stats):
        tiles = []
        for verdict_key, label in (('FASTQ Usability', 'FASTQ usability'),
                                   ('Reference Usability', 'Reference usability')):
            value = stats.get(verdict_key)
            if value:
                # "Questionable (average depth 35.2x < 40.0x)": the verdict is
                # the tile value, the measurement that set it is the note.
                word, _, note = str(value).strip().partition(' (')
                tiles.append({'label': label, 'value': word,
                              'note': note.rstrip(')') if note else '',
                              'cls': VERDICT_CLASS.get(word.strip().lower(), '')})
        for key, label in (('Average Depth', 'Average depth'),
                           ('Genome with Coverage', 'Genome with coverage'),
                           ('Quality SNPs', 'Quality SNPs')):
            if key in stats and stats[key] not in (None, ''):
                tiles.append({'label': label, 'value': stats[key], 'cls': ''})
        return tiles

    def write(self, stats):
        '''
        Render the report from the stats dict that also produces the workbook.

        Returns the path written.  Any key not claimed by a section appears under
        "Other metrics", so a metric added to excel_dict later shows up in the
        report without anyone having to remember to add it here.

        Every value is passed through redact_paths on the way in, so a metric
        added later cannot reintroduce an absolute path into a shareable report
        without anyone noticing.
        '''
        stats = {k: redact_paths(v) for k, v in (stats or {}).items()}
        claimed = set(TILE_KEYS) | {'sample', 'date', 'FASTQ Usability',
                                    'Reference Usability'}

        paired = bool(str(stats.get('FASTQ_R2') or '').strip())
        read_rows = []
        for label, k1, k2 in READ_METRICS:
            claimed.update((k1, k2))
            v1, v2 = stats.get(k1, ''), stats.get(k2, '')
            if v1 in (None, '') and v2 in (None, ''):
                continue
            read_rows.append((label,
                              self._fmt(label, v1),
                              self._fmt(label, v2) if paired else '—'))

        sections = []
        for title, entries in SECTIONS:
            rows = []
            for label, key in entries:
                claimed.add(key)
                value = stats.get(key)
                if value in (None, ''):
                    continue
                mono = key in ('BAM File', 'BAM/Reference File',
                               'Spoligotype Binary Code', 'Spoligotype Spacer Counts')
                rows.append((label, value, mono))
            if rows:
                sections.append((title, rows))

        other = [(k, v) for k, v in sorted(stats.items())
                 if k not in claimed and v not in (None, '')]

        html = _template().render(
            sample=self.sample_name,
            description=self.description,
            version=__version__,
            generated=datetime.now().strftime('%B %d, %Y at %H:%M'),
            run_stamp=self.date_stamp,
            logo_data_uri=self._logo_data_uri(),
            tiles=self._tiles(stats),
            notices=[redact_paths(n) for n in self.notices],
            paired=paired,
            read_rows=read_rows,
            sections=sections,
            tables=[{'title': t['title'],
                     'headers': t['headers'],
                     'rows': [[redact_paths(c) for c in row] for row in t['rows']]}
                    for t in self.tables],
            other_rows=other,
            program_versions=redact_paths(self.program_versions),
            runtime=self.runtime,
        )
        with open(self.html_file, 'w', encoding='utf-8') as out:
            out.write(html)
        return self.html_file

    # --------------------------------------------------------------------- pdf

    def to_pdf(self, html_file=None):
        '''
        Render the HTML report to PDF with WeasyPrint.

        Returns the PDF path, or None after printing why not.  A missing or
        failing converter is reported and does not raise: the HTML report is the
        primary artifact and is already on disk by this point, so aborting the
        run over the second copy of it would lose the alignment and the VCF for
        no reason.  Unlike the pdflatex call this replaced, the result is
        checked -- that one discarded both its exit status and its output, so a
        machine with no LaTeX produced no report and said nothing.
        '''
        html_file = html_file or self.html_file
        if not os.path.exists(html_file):
            print(f'Note: no PDF written, {html_file} is not on disk.')
            return None
        try:
            from weasyprint import HTML
        except ImportError as e:
            print(f'Note: no PDF written ({e}).  The HTML report is complete.\n'
                  f'      Install the renderer with:  '
                  f'conda install -c conda-forge weasyprint')
            return None
        # WeasyPrint logs every CSS property it does not implement at WARNING.
        # None of them change this report, and on a 500-sample run the noise
        # buries the messages that matter.
        logging.getLogger('weasyprint').setLevel(logging.ERROR)
        try:
            HTML(filename=html_file).write_pdf(self.pdf_file)
        except Exception as e:                            # noqa: BLE001
            print(f'Warning: HTML report written but PDF conversion failed: {e}')
            return None
        if not os.path.exists(self.pdf_file):
            print('Warning: WeasyPrint reported success but wrote no PDF.')
            return None
        return self.pdf_file

    @staticmethod
    def _fmt(label, value):
        if value in (None, ''):
            return '—'
        if label == 'File':
            return os.path.basename(str(value))
        return value


# Workers a pool may start when VSNP3_MAX_CPUS is unset.  The vSNP GUI's
# install step has capped vsnp3's pools at this value since 2026-05; the
# number comes from a 128-core shared server where cpu_count / 1.2 started
# 106 processes.  A pool larger than this gains nothing here: reading is I/O
# bound, and RAxML concurrency is limited separately (vsnp3_fasta_to_snps_table).
DEFAULT_MAX_WORKERS = 32


def worker_count(tasks=None):
    '''
    Size of a worker pool: cpu_count / 1.2, capped by VSNP3_MAX_CPUS when it is
    set to a positive whole number and by DEFAULT_MAX_WORKERS otherwise, and
    never more than tasks when given.

    cpu_count() reports the machine, not the allocation.  Inside an 8-core job
    on a 64-core node, or on a 128-core shared server, cpu_count / 1.2 starts
    far more processes than the run was given, and each one holds its share of
    the VCF frames -- a 1,951-sample step 2 in a 32 GB session was killed with
    no traceback.  The vSNP GUI's suite launcher passes the session's core
    budget as VSNP3_MAX_CPUS; vsnp3_fasta_to_snps_table.py reads the same
    variable to limit concurrent RAxML runs.  Until 3.37 only that file read
    it, and the GUI's install step rewrote the pool-size expression in
    vsnp3_step2.py and vsnp3_group_on_defining_snps.py to
    min(cpu_count / 1.2, VSNP3_MAX_CPUS or 32); this is that expression, so a
    3.37 install behaves as the patched 3.36 installs did.
    '''
    count = max(1, int(multiprocessing.cpu_count() / 1.2))
    try:
        cap = int(str(os.environ.get('VSNP3_MAX_CPUS') or '').strip() or 0)
    except ValueError:
        cap = 0
    count = min(count, cap if cap > 0 else DEFAULT_MAX_WORKERS)
    if tasks is not None:
        count = min(count, max(1, int(tasks)))
    return max(1, count)


# Shown as verdict tiles in the step 1 HTML report, never written to the stats
# workbook.  The workbooks are merged across samples and sorted on, and a one
# word judgement next to the measurements it was derived from is not a
# measurement.  The report shows the verdict with the threshold that set it.
HTML_ONLY_KEYS = ('FASTQ Usability', 'Reference Usability')


def _usability(limits, checks):
    '''
    (verdict, reason) from a limits table and the measurements.

    checks: key -> (value, crosses_when_below, label, format).  The first
    limit crossed decides, Poor before Questionable, in the order the limits
    are written in vsnp3_version.py.
    '''
    for verdict in ('poor', 'questionable'):
        for key, limit in limits[verdict].items():
            if key not in checks:
                continue
            value, below, label, fmt = checks[key]
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if (value < limit) if below else (value > limit):
                op = '<' if below else '>'
                return verdict.capitalize(), f'{label} {fmt.format(value)} {op} {fmt.format(limit)}'
    return 'Acceptable', 'all limits met'


def fastq_usability_verdict(depth, duplicate_ratio, q20_percent, q30_percent):
    '''FASTQ usability as (verdict, reason), limits from FASTQ_USABILITY.'''
    checks = {
        'depth': (depth, True, 'average depth', '{:.1f}x'),
        'duplicate_ratio': (duplicate_ratio, False, 'duplicate reads', '{:.0%}'),
        'q20_percent': (q20_percent, True, 'reads passing Q20', '{:.1f}%'),
        'q30_percent': (q30_percent, True, 'reads passing Q30', '{:.1f}%'),
    }
    return _usability(FASTQ_USABILITY, checks)


def reference_usability_verdict(genome_coverage, zero_coverage_percent, good_snp_percent):
    '''Reference usability as (verdict, reason), limits from REFERENCE_USABILITY.'''
    checks = {
        'genome_coverage': (genome_coverage, True, 'reference with coverage', '{:.1%}'),
        'zero_coverage_percent': (zero_coverage_percent, False, 'reference with zero coverage', '{:.2f}%'),
        'good_snp_percent': (good_snp_percent, False, 'quality SNPs, percent of reference', '{:.3f}%'),
    }
    return _usability(REFERENCE_USABILITY, checks)


class Excel_Stats:

    def __init__(self, sample_name):
        self.sample_name = sample_name
        date_stamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.excel_filename = '{}_{}_stats.xlsx'.format(sample_name, date_stamp)
        excel_dict = {}
        excel_dict['sample'] = sample_name
        excel_dict['date'] = date_stamp
        self.excel_dict = excel_dict

    def post_excel(self,):
        rows = {k: v for k, v in self.excel_dict.items() if k not in HTML_ONLY_KEYS}
        df = pd.DataFrame.from_dict(rows, orient='index').T
        df = df.set_index('sample')
        df.to_excel(self.excel_filename)

# Created 2021 by Tod Stuber