#!/usr/bin/env python3
"""
Copy conda_build/ to a temporary recipe that builds from this checkout.

    python3 tests/local_source_recipe.py <output_dir>

conda_build/meta.yaml fetches
`https://github.com/USDA-VS/vsnp3/archive/{{ version }}.tar.gz` and verifies a
sha256.  On a branch that tarball does not exist and the hash is a placeholder, so
conda-build cannot run against the recipe as written -- which is why the build was
never exercised before a bioconda PR.

Only the `source:` block is rewritten, to `path:` pointing at the checkout.
Everything else -- build.sh, the run requirements, the test commands -- is used
exactly as the maintainer wrote it, so this checks the things that actually break:
whether the environment solves, whether build.sh installs everything, and whether
the recipe's own tests pass against the installed package.
"""

import os
import re
import shutil
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The source: block runs to the next top-level key.
SOURCE_BLOCK = re.compile(r'^source:.*?(?=^\S)', re.M | re.S)


def main():
    if len(sys.argv) != 2:
        print(__doc__.strip())
        return 1
    out = os.path.abspath(sys.argv[1])

    src_recipe = os.path.join(REPO, 'conda_build')
    if not os.path.isdir(src_recipe):
        print(f'No recipe directory at {src_recipe}', file=sys.stderr)
        return 1

    if os.path.exists(out):
        shutil.rmtree(out)
    shutil.copytree(src_recipe, out)

    meta_path = os.path.join(out, 'meta.yaml')
    with open(meta_path) as f:
        meta = f.read()

    if not SOURCE_BLOCK.search(meta):
        print('Could not find a source: block in meta.yaml', file=sys.stderr)
        return 1
    meta = SOURCE_BLOCK.sub(f'source:\n  path: {REPO}\n\n', meta)

    # The sha256 placeholder is only meaningful for the url form.
    meta = re.sub(r'^\{%\s*set sha256.*$\n', '', meta, flags=re.M)

    with open(meta_path, 'w') as f:
        f.write(meta)

    print(f'recipe: {out}')
    print(f'source: {REPO}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
