#!/usr/bin/env python3
"""
Checks that must hold before a release is tagged.

    python3 tests/test_release_consistency.py
    VSNP3_VERIFY_SHA256=1 python3 tests/test_release_consistency.py   # after tagging

The version number lives in bin/vsnp3_version.py, but a dozen other files carry
it too and cannot import Python: the conda recipe, the README, the container
definition and its README, and the internal/*.sh wrappers.  Before this check they
drifted badly -- meta.yaml said 3.23 while the code said 3.35.

Runs in a release folder with no .git as well as in a checkout; what cannot be
checked outside a checkout is reported as WARN, never silently skipped.
"""

import fnmatch
import os
import re
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PASS, FAIL, WARN = [], [], []


def check(name, ok, detail=''):
    (PASS if ok else FAIL).append(name)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f'\n        {detail}' if detail and not ok else ''))
    return ok


def warn(name, detail):
    WARN.append(name)
    print(f'  WARN  {name}\n        {detail}')


def read(*parts):
    path = os.path.join(REPO, *parts)
    if not os.path.exists(path):
        return None
    with open(path, errors='replace') as f:
        return f.read()


def gitignored(rel, patterns):
    """
    Whether a .gitignore line in `patterns` matches the relative path `rel`.

    Only for a folder with no .git, where `git check-ignore` cannot answer.  It
    covers the forms this repository's .gitignore uses (names, globs, `dir/`,
    `**/name`), not the full gitignore grammar.
    """
    parts = rel.split(os.sep)
    for pattern in patterns:
        p = pattern.rstrip('/')
        if p.startswith('**/'):
            p = p[3:]
        if '/' in p.lstrip('/'):
            if fnmatch.fnmatch(rel, p.lstrip('/')):
                return True
        elif any(fnmatch.fnmatch(part, p) for part in parts):
            return True
    return False


def canonical_version():
    body = read('bin', 'vsnp3_version.py')
    return re.search(r'__version__ = "([^"]+)"', body).group(1)


def specs(block):
    """
    Parse a conda requirements block into {package: version_spec}.

    The spec is the part after the name, normalized only for whitespace, so
    'numpy >=1.26.4,<3' -> {'numpy': '>=1.26.4,<3'} and a bare '- pigz' ->
    {'pigz': ''}.  Keeping the spec rather than discarding it is the point:
    package names alone cannot see a pin that moved in one file and not another.
    """
    out = {}
    for line in block.splitlines():
        line = line.strip()
        if not line.startswith('- '):
            continue
        item = line[2:].split('#', 1)[0].strip()
        if not item:
            continue
        m = re.match(r'([A-Za-z0-9._-]+)\s*(.*)$', item)
        if m:
            out[m.group(1)] = re.sub(r'\s+', '', m.group(2))
    return out


def tarball_sha256(url, timeout=120):
    """sha256 of what `url` serves, streamed; raises urllib errors as they come."""
    import hashlib
    import urllib.request
    digest = hashlib.sha256()
    with urllib.request.urlopen(url, timeout=timeout) as response:
        for chunk in iter(lambda: response.read(1 << 16), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    print('=' * 78)
    version = canonical_version()
    print(f'vSNP3 release consistency -- canonical version {version}')
    print('=' * 78)

    print('\n[1] version agrees everywhere it is written')
    # (file, regex capturing the version, must the file exist).  Every file here is
    # one internal/bump_version.sh rewrites, so a mismatch is never legitimate and
    # is a FAIL; `required` only says whether the file's absence is one too.
    targets = [
        ('conda_build/meta.yaml', r'{%\s*set version\s*=\s*"([^"]+)"\s*%}', True),
        ('README.md', r'vsnp3=([0-9.]+)', True),
        ('container/Singularity.def', r'vsnp3=([0-9.]+)', False),
        ('container/README.md', r'vsnp3(?:=|:v)([0-9.]+)', False),
    ]
    for rel, pattern, required in targets:
        body = read(*rel.split('/'))
        if body is None:
            if required:
                check(f'{rel} present', False, f'{rel} is missing')
            else:
                warn(f'{rel} not present', 'not release-blocking; it is documentation')
            continue
        found = sorted(set(re.findall(pattern, body)))
        check(f'{rel} says {version}', found == [version],
              f'found {found}; internal/bump_version.sh rewrites this file, so run it '
              f'(or fix the line by hand)')

    # The shell wrappers carry version="X.YZ" (some "vsnp3_X.YZ") and cannot import
    # Python.  bump_version.sh rewrites exactly that line in every internal/*.sh, so
    # a wrapper that disagrees was added or edited after the bump: a FAIL, with no
    # legitimate exception.  internal/attic/ is not scanned; nothing there runs.
    for directory in ('internal',):
        d = os.path.join(REPO, directory)
        if not os.path.isdir(d):
            continue
        stale, literals, wrappers = [], [], 0
        for name in sorted(os.listdir(d)):
            if not name.endswith('.sh'):
                continue
            wrappers += 1
            body = read(directory, name)
            for found in re.findall(r'^version="(?:vsnp3_)?([0-9.]+)"', body, re.M):
                if found != version:
                    stale.append(f'{directory}/{name}: version="{found}"')
            # Hint lines used to hard-code env names such as vsnp3_3.26 and
            # vsnp3=3.31; they now derive from ${version}.  A number that differs
            # from the canonical one is worth a look but may be a legitimate
            # mention of an older release in a comment, so WARN, not FAIL.
            for found in re.findall(r'vsnp3(?:[=_]|:v)([0-9]+\.[0-9]+)\b', body):
                if found != version:
                    literals.append(f'{directory}/{name}: vsnp3 {found}')
        check(f'every version= line in {directory}/*.sh says {version} ({wrappers} wrappers)',
              not stale, '\n        '.join(stale) +
              '\n        internal/bump_version.sh rewrites these lines; run it')
        if literals:
            warn(f'{len(literals)} other vsnp3 version literal(s) in {directory}/*.sh differ '
                 f'from {version}', '\n        '.join(literals))

    print('\n[2] conda recipe')
    meta = read('conda_build', 'meta.yaml')
    if meta:
        license_file = re.search(r'license_file:\s*(\S+)', meta)
        if license_file:
            name = license_file.group(1)
            check(f'declared license_file "{name}" exists',
                  os.path.exists(os.path.join(REPO, name)),
                  f'meta.yaml declares it; bioconda resolves it from the source '
                  f'tarball, so this repo and the tarball differ')
        # Anything imported at module scope by an entry point must be declared, and
        # so must the programs the tools shell out to.
        run_block = meta.split('run:', 1)[-1].split('test:', 1)[0] if 'run:' in meta else ''
        for package, why in [('pyyaml', 'imported at module scope by bin/vsnp3_file_setup.py, '
                                        'which both entry points import'),
                             ('jinja2', 'imported lazily by HtmlReport in bin/vsnp3_file_setup.py '
                                        'for the step 1 report'),
                             ('weasyprint', 'imported lazily by HtmlReport.to_pdf in '
                                            'bin/vsnp3_file_setup.py to render the report PDF'),
                             ('matplotlib', 'imported by bin/vsnp3_kernel_plots.py'),
                             ('seaborn', 'imported by bin/vsnp3_kernel_plots.py'),
                             ('snp-dists', 'run with subprocess by bin/vsnp3_kernel_plots.py; '
                                           'it is the program that computes the distances plotted')]:
            check(f'{package} declared in meta.yaml run requirements',
                  package in run_block.lower(), why)

        # The sha256 is of the GitHub tag tarball and so cannot exist before the tag.
        # internal/bump_version.sh writes a placeholder here on purpose, so that a
        # release cannot go out carrying the previous release's hash: conda-build
        # would then fail with a hash mismatch that names nothing useful, and the
        # 3.37 recipe did sit with the 3.36 hash until an audit noticed.
        sha = re.search(r'\{%\s*set sha256\s*=\s*"([^"]+)"', meta)
        if sha:
            value = sha.group(1)
            url = f'https://github.com/USDA-VS/vsnp3/archive/{version}.tar.gz'
            regenerate = f'curl -L {url} | shasum -a 256'
            is_hex = bool(re.fullmatch(r'[0-9a-f]{64}', value))
            check('meta.yaml sha256 is a 64-hex hash, not the bump placeholder', is_hex,
                  f'sha256 is "{value}".  After tagging {version}, regenerate it and write '
                  f'it into meta.yaml:\n        {regenerate}')
            if is_hex:
                if os.environ.get('VSNP3_VERIFY_SHA256'):
                    import urllib.error
                    try:
                        actual = tarball_sha256(url)
                    except urllib.error.HTTPError as e:
                        warn(f'meta.yaml sha256 could not be verified: {url} returned {e.code}',
                             f'the {version} tag is not published yet (404) or GitHub is '
                             f'unavailable; verify after tagging with\n        {regenerate}')
                    except (urllib.error.URLError, OSError, TimeoutError) as e:
                        warn('meta.yaml sha256 could not be verified: no network',
                             f'{e}; verify after tagging with\n        {regenerate}')
                    else:
                        check(f'meta.yaml sha256 matches the {version} tag tarball',
                              actual == value,
                              f'meta.yaml has {value}\n        tarball is   {actual}\n'
                              f'        write the tarball value into meta.yaml')
                else:
                    warn(f'meta.yaml sha256 not verified against the {version} tag tarball',
                         f'a 64-hex value passes the format check even when it is the '
                         f'previous release\'s hash.  After tagging, run\n'
                         f'        VSNP3_VERIFY_SHA256=1 {os.path.basename(sys.argv[0])}\n'
                         f'        (downloads the tarball and compares) or\n        {regenerate}')

        # conda_build/vsnp3.yaml is the development environment for running a
        # checkout without installing the package.  It had drifted from meta.yaml by
        # five dependencies, so an environment built from it was missing jinja2,
        # pyyaml, matplotlib-base, seaborn and vcflib.  Specs are compared, not just
        # names: the versions are the half that drifted silently, because comparing
        # names alone cannot see a pin that moved in one file and not the other.
        meta_specs = specs(run_block)
        env = read('conda_build', 'vsnp3.yaml')
        if env:
            env_specs = specs(env.split('dependencies:', 1)[-1])
            missing = sorted(set(meta_specs) - set(env_specs))
            extra = sorted(set(env_specs) - set(meta_specs))
            check('conda_build/vsnp3.yaml declares the same packages as meta.yaml',
                  not missing and not extra,
                  f'missing from vsnp3.yaml: {missing}; not in meta.yaml: {extra}')
            differing = [f'{p}: meta.yaml "{meta_specs[p]}" vs vsnp3.yaml "{env_specs[p]}"'
                         for p in sorted(set(meta_specs) & set(env_specs))
                         if meta_specs[p] != env_specs[p]]
            check('conda_build/vsnp3.yaml pins the same versions as meta.yaml',
                  not differing, '\n        '.join(differing))

        # `<=` does not mean what it looks like.  Conda pads the shorter version
        # with zeros before comparing, so 1.26.4 > 1.26 and `numpy <=1.26` admits
        # only 1.26.0 exactly.  Three pins were written this way and nobody noticed
        # until an environment turned out to be held on numpy 1.26.0 and Python
        # 3.12.0.  Cap at the next major instead.
        for rel, block in (('conda_build/meta.yaml', run_block),
                           ('conda_build/vsnp3.yaml',
                            (env or '').split('dependencies:', 1)[-1])):
            offenders = [f'{p} {s}' for p, s in specs(block).items() if '<=' in s]
            check(f'{rel} uses no <= version caps', not offenders,
                  f'{", ".join(offenders)}\n        '
                  f'conda zero-pads, so "<=1.26" excludes 1.26.4 and means exactly '
                  f'1.26.0. Use "<2" or "<1.27".')

        # An exact pin on a package reached only through another one breaks the
        # solve the day that other package moves; vcflib=1.0.12 did exactly that.
        # The recipe's policy is floors (and major-version ceilings) only.
        exact = [f'{p} {s}' for p, s in meta_specs.items()
                 if re.fullmatch(r'={1,2}[0-9][^,]*', s)]
        check('meta.yaml has no exact (=) pins', not exact,
              f'{", ".join(exact)}\n        use a floor (>=) per the policy comment in meta.yaml')

        # A bare package name lets the solver pick anything, which for a compiled
        # numpy consumer means it may pair an old ABI with a new numpy.  These are
        # the ones where that matters.
        needs_floor = ['python', 'numpy', 'pandas', 'scipy', 'biopython',
                       'matplotlib-base', 'seaborn', 'dask']
        bare = [p for p in needs_floor
                if p in meta_specs and '>=' not in meta_specs[p]]
        check('every version-sensitive package declares a floor in meta.yaml',
              not bare,
              f'no lower bound on: {", ".join(bare)}\n        '
              f'these constrain or are constrained by the numpy ABI')

        spdx = re.search(r'^\s*license:\s*(\S+)', meta, re.M)
        check('meta.yaml license is an SPDX identifier',
              bool(spdx) and re.fullmatch(r'[A-Za-z0-9.+-]+(-or-later|-only)?', spdx.group(1))
              is not None and spdx.group(1) != 'GPL3',
              f'license is {spdx.group(1) if spdx else "missing"}; bioconda wants an SPDX id '
              f'such as GPL-3.0-or-later, with license_family for the old name')

        # The two CI install lines are the fourth and fifth copies of these pins.
        # The conda line had already lost the `pandas >=1.3` floor the pip line
        # carried, with nothing to notice.
        ci = read('.github', 'workflows', 'ci.yml')
        if ci:
            declared = {p: meta_specs[p] for p in ('numpy', 'pandas', 'biopython')
                        if p in meta_specs}
            for package, spec in declared.items():
                # Matrix legs deliberately pin numpy and pandas to one end of the
                # declared range, so require only that CI names the package.  That
                # each leg resolves *inside* the range is checked at runtime by
                # tests/test_numpy_contract.py, which is the only place that can
                # see what actually got installed.
                check(f'ci.yml constrains {package}',
                      re.search(rf'[\'"]?{re.escape(package)}[\'"]?\s*[><=]', ci) is not None,
                      f'meta.yaml says "{spec}" but ci.yml does not pin {package}')
        else:
            # A release folder copied out of the repository has no .github/.  The
            # pin checks cannot run here, and meta.yaml's "floors are the oldest
            # versions CI exercises" rests on that file.
            warn('.github/workflows/ci.yml not present: CI pin checks skipped',
                 'run this in the repository checkout as well, where ci.yml exists; '
                 'tests/make_test_env.sh falls back to conda_build/vsnp3.yaml meanwhile')

    print('\n[3] no assistant or AI-tooling references')
    # These must not reach GitHub or GitLab.  Checked here rather than trusted,
    # because the first audit of this release missed nine commits that were already
    # pushed: only the new commits had been looked at.
    #
    # The words are assembled at run time rather than written out.  This file is
    # scanned by its own check, and the PASS/FAIL lines below are archived in
    # validation records, so neither the source nor the output may carry them.
    patterns = tuple(''.join(p) for p in (
        ('cl', 'aude'), ('anthr', 'opic'), ('cop', 'ilot'), ('chat', 'gpt'),
        ('open', 'ai'), ('co-auth', 'ored-by'), ('generat', 'ed with')))
    # The local working-notes file: the one place such notes may live, and it must
    # never be committed.  .gitignore covers it with a glob, so its name is not
    # written there either.
    notes_file = ''.join(('cl', 'aude', '-handoff.txt'))

    # A release folder copied out of the repository has no .git.  git then answers
    # every question below with "not a git repository", which used to read as an
    # empty file list (so the reference scan passed without reading anything) and
    # as "not ignored".  Outside a checkout, scan every file .gitignore would let
    # through and read .gitignore itself; commit history cannot be checked there.
    in_git = subprocess.run(['git', '-C', REPO, 'rev-parse', '--is-inside-work-tree'],
                            capture_output=True, text=True,
                            check=False).stdout.strip() == 'true'
    ignore_patterns = [line.strip() for line in (read('.gitignore') or '').splitlines()
                       if line.strip() and not line.lstrip().startswith(('#', '!'))]
    if in_git:
        tracked = subprocess.run(['git', '-C', REPO, 'ls-files'],
                                 capture_output=True, text=True, check=False).stdout.split()
    else:
        tracked = []
        for root, dirs, files in os.walk(REPO):
            dirs[:] = [d for d in dirs if d != '.git' and not gitignored(
                os.path.relpath(os.path.join(root, d), REPO), ignore_patterns)]
            for name in files:
                rel = os.path.relpath(os.path.join(root, name), REPO)
                if not gitignored(rel, ignore_patterns):
                    tracked.append(rel)
    offenders = []
    for rel in tracked:
        path = os.path.join(REPO, rel)
        try:
            with open(path, errors='replace') as f:
                body = f.read().lower()
        except (OSError, IsADirectoryError):
            continue
        hits = sorted({p for p in patterns if p in body})
        if hits:
            offenders.append(f'{rel}: {", ".join(hits)}')
    check('no assistant references in tracked files', not offenders,
          '\n        '.join(offenders))

    if in_git:
        log = subprocess.run(['git', '-C', REPO, 'log', '--format=%B', '-n', '400'],
                             capture_output=True, text=True, check=False)
        body = log.stdout.lower()
        in_messages = sorted({p for p in patterns if p in body})
        check('no assistant references in commit messages', not in_messages,
              f'found {in_messages} in the last 400 commit messages; strip the trailer '
              f'lines and rewrite, or the references reach the remote')

        people = subprocess.run(['git', '-C', REPO, 'log', '--format=%an|%ae|%cn|%ce',
                                 '-n', '400'], capture_output=True, text=True, check=False)
        in_people = sorted({p for p in patterns if p in people.stdout.lower()})
        check('no assistant in author or committer fields', not in_people, f'{in_people}')
    else:
        warn('not a git checkout: commit messages and authors not checked',
             'run this again in the release repository before tagging')

    if in_git:
        ignored = subprocess.run(['git', '-C', REPO, 'check-ignore', '-q', notes_file],
                                 capture_output=True, check=False).returncode == 0
    else:
        ignored = gitignored(notes_file, ignore_patterns)
    is_tracked = notes_file in tracked
    check('working-notes file is gitignored', ignored,
          'add the `*-handoff.txt` rule back to .gitignore; that file is the only place '
          'such notes may live')
    check('working-notes file is not tracked', not is_tracked,
          'the working-notes file is committed and would reach the remote')

    print('\n[4] repository hygiene')
    big = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__', '.release')]
        for name in files:
            path = os.path.join(root, name)
            try:
                size = os.path.getsize(path)
            except OSError:
                continue
            if size > 10 * 1024 * 1024:
                big.append(f'{os.path.relpath(path, REPO)}: {size / 1e6:.0f} MB')
    check('no file over 10 MB in the tree', not big,
          'bioconda builds from the GitHub source tarball, so every user '
          'downloads these:\n        ' + '\n        '.join(big))

    bin_dir = os.path.join(REPO, 'bin')
    tests_in_bin = [f for f in os.listdir(bin_dir) if f.startswith('test_')]
    check('no test files in bin/', not tests_in_bin,
          f'conda_build/build.sh does `mv bin/*.py`, so these ship to users: {tests_in_bin}')

    # Gitignored, so harmless in a checkout; but a release folder tarred directly
    # ships them, and Finder recreates them whenever the folder is browsed.
    strays = []
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d != '.git']
        strays += [os.path.relpath(os.path.join(root, f), REPO) for f in files
                   if f == '.DS_Store']
        strays += [os.path.relpath(os.path.join(root, d), REPO) for d in dirs
                   if d == '__pycache__']
    if strays:
        warn(f'{len(strays)} stray file(s) or folder(s) in the tree',
             ', '.join(sorted(strays)) + '\n        delete before tarring a release folder')

    print('\n' + '=' * 78)
    print(f'PASS {len(PASS)}   FAIL {len(FAIL)}   WARN {len(WARN)}')
    for name in FAIL:
        print(f'  FAIL {name}')
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
