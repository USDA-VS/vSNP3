#!/usr/bin/env python3
"""
Checks that must hold before a release is tagged.

    python3 tests/test_release_consistency.py

The version number lives in bin/vsnp3_version.py, but a dozen other files carry
it too and cannot import Python: the conda recipe, the README, the container
definition and the internal/*.sh wrappers.  Before this check they drifted badly
-- meta.yaml said 3.23 while the code said 3.35.
"""

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


def main():
    print('=' * 78)
    version = canonical_version()
    print(f'vSNP3 release consistency -- canonical version {version}')
    print('=' * 78)

    print('\n[1] version agrees everywhere it is written')
    # (file, regex capturing the version, required)
    targets = [
        ('conda_build/meta.yaml', r'{%\s*set version\s*=\s*"([^"]+)"\s*%}', True),
        ('README.md', r'vsnp3=([0-9.]+)', True),
        ('container/Singularity.def', r'vsnp3=([0-9.]+)', False),
    ]
    for rel, pattern, required in targets:
        body = read(*rel.split('/'))
        if body is None:
            check(f'{rel} present', not required, f'{rel} is missing')
            continue
        found = sorted(set(re.findall(pattern, body)))
        ok = found == [version]
        if required:
            check(f'{rel} says {version}', ok, f'found {found}')
        elif not ok:
            warn(f'{rel} says {found}, canonical is {version}',
                 'not release-blocking, but it is user-facing documentation')

    # The shell wrappers carry version="X.YZ" and cannot import Python.
    for directory in ('internal',):
        d = os.path.join(REPO, directory)
        if not os.path.isdir(d):
            continue
        stale = []
        for name in sorted(os.listdir(d)):
            if not name.endswith('.sh'):
                continue
            body = read(directory, name)
            for found in re.findall(r'^version="(?:vsnp3_)?([0-9.]+)"', body, re.M):
                if found != version:
                    stale.append(f'{directory}/{name}: {found}')
        if stale:
            warn(f'{len(stale)} wrapper(s) in {directory}/ carry a different version',
                 '\n        '.join(stale))

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
        # Anything imported at module scope by an entry point must be declared.
        run_block = meta.split('run:', 1)[-1].split('test:', 1)[0] if 'run:' in meta else ''
        for package, why in [('pyyaml', 'imported at module scope by bin/vsnp3_file_setup.py, '
                                        'which both entry points import'),
                             ('jinja2', 'imported lazily by HtmlReport in bin/vsnp3_file_setup.py '
                                        'for the step 1 report'),
                             ('weasyprint', 'imported lazily by HtmlReport.to_pdf in '
                                            'bin/vsnp3_file_setup.py to render the report PDF'),
                             ('matplotlib', 'imported by bin/vsnp3_kernel_plots.py'),
                             ('seaborn', 'imported by bin/vsnp3_kernel_plots.py')]:
            check(f'{package} declared in meta.yaml run requirements',
                  package in run_block.lower(), why)

        # A placeholder sha256 must not reach a tag: conda-build would fail with a
        # hash mismatch rather than anything that names the real problem.
        sha = re.search(r'\{%\s*set sha256\s*=\s*"([^"]+)"', meta)
        if sha:
            value = sha.group(1)
            check('meta.yaml sha256 is a real hash',
                  bool(re.fullmatch(r'[0-9a-f]{64}', value)),
                  f'sha256 is "{value}". Regenerate against the tag tarball:\n'
                  f'        curl -sL https://github.com/USDA-VS/vsnp3/archive/'
                  f'{version}.tar.gz | shasum -a 256')

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

    print('\n[3] no assistant or AI-tooling references')
    # These must not reach GitHub or GitLab.  Checked here rather than trusted,
    # because the first audit of this release missed nine commits that were already
    # pushed: only the new commits had been looked at.
    patterns = ('claude', 'anthropic', 'copilot', 'chatgpt', 'openai',
                'co-authored-by', 'generated with')
    handoff = 'claude-handoff.txt'

    tracked = subprocess.run(['git', '-C', REPO, 'ls-files'],
                             capture_output=True, text=True, check=False)
    offenders = []
    for rel in tracked.stdout.split():
        if rel == os.path.relpath(__file__, REPO):
            continue                      # this file names the patterns on purpose
        path = os.path.join(REPO, rel)
        try:
            with open(path, errors='replace') as f:
                body = f.read().lower()
        except (OSError, IsADirectoryError):
            continue
        # .gitignore has to name the handoff file in order to ignore it, and this
        # check has to name it in order to check for it.  Only that exact filename
        # is exempt -- any other use of the words still fails.
        body = body.replace(handoff, '')
        hits = sorted({p for p in patterns if p in body})
        if hits:
            offenders.append(f'{rel}: {", ".join(hits)}')
    check('no assistant references in tracked files', not offenders,
          '\n        '.join(offenders))

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

    # The one permitted location, and it must never be committed.
    ignored = subprocess.run(['git', '-C', REPO, 'check-ignore', '-q', handoff],
                             capture_output=True, check=False).returncode == 0
    is_tracked = handoff in tracked.stdout.split()
    check(f'{handoff} is gitignored', ignored,
          f'add {handoff} to .gitignore; it is the only place such notes may live')
    check(f'{handoff} is not tracked', not is_tracked,
          f'{handoff} is committed and would reach the remote')

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

    print('\n' + '=' * 78)
    print(f'PASS {len(PASS)}   FAIL {len(FAIL)}   WARN {len(WARN)}')
    for name in FAIL:
        print(f'  FAIL {name}')
    print('=' * 78)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
