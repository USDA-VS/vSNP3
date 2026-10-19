#!/usr/bin/env bash
#
# Build one leg of the CI version matrix as a local conda environment.
#
#   tests/make_test_env.sh <leg> [env_name] [-- extra conda args]
#   tests/make_test_env.sh --print-specs <leg>
#   tests/make_test_env.sh --list
#
# The legs are read from .github/workflows/ci.yml and the dependency list from
# conda_build/vsnp3.yaml, so an environment built here is the same environment
# CI tests -- there is no seventh copy of the pins to drift.  That mattered:
# before 3.36 the pins existed in five places and the CI conda line had already
# lost a floor the pip line still carried.
#
# The leg `recipe` is always available and needs no ci.yml: it is
# conda_build/vsnp3.yaml exactly as declared, with no extra pins.  A release
# folder copied out of the repository has no .github/, so that is the only leg
# there (and the script says so rather than failing).
#
# Typical use, testing both ends of the declared numpy range:
#
#   tests/make_test_env.sh floor   vsnp3_floor
#   tests/make_test_env.sh current vsnp3_current
#   tests/make_test_env.sh recipe  vsnp3_recipe
#
# Only conda-forge and bioconda are used, with `defaults` explicitly excluded --
# it is not needed, it slows the solve, and on institutional machines its terms
# of service are a question nobody wants to answer during a release.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/.." && pwd)

CI="$repo/.github/workflows/ci.yml"
ENV_YAML="$repo/conda_build/vsnp3.yaml"

[ -r "$ENV_YAML" ] || { echo "cannot read $ENV_YAML" >&2; exit 1; }
if [ ! -r "$CI" ]; then
    echo "note: $CI is not here (a release folder has no .github/);" >&2
    echo "      only the 'recipe' leg -- conda_build/vsnp3.yaml as declared -- is available" >&2
    CI=""
fi

# Parsed with the standard library only: this has to run before any environment
# exists, so pyyaml cannot be assumed.  Both files are plain enough for it.
read_specs() {
    python3 - "$CI" "$ENV_YAML" "$1" <<'PY'
import os, re, sys

ci_path, env_path, leg = sys.argv[1:4]

# The recipe leg pins nothing beyond the environment file itself.  CI legs
# follow, from the integration matrix: `- leg: floor` then python/numpy/pandas.
legs, cur = {'recipe': {}}, None
if ci_path and os.path.exists(ci_path):
    inside = False
    for line in open(ci_path):
        if re.match(r'^  \w[\w-]*:\s*$', line):
            inside = line.strip() == 'integration:'
        if not inside:
            continue
        m = re.match(r'\s*-\s*leg:\s*(\S+)\s*$', line)
        if m:
            cur = m.group(1)
            legs[cur] = {}
            continue
        m = re.match(r'\s*(python|numpy|pandas):\s*(.+?)\s*$', line)
        if m and cur:
            legs[cur][m.group(1)] = m.group(2).strip().strip("'\"")

if leg == '--list':
    for name, pins in legs.items():
        detail = '  '.join(f'{k}{v}' if v.startswith(('=', '>', '<')) else f'{k} {v}'
                           for k, v in pins.items())
        print(f"{name:9s} " + (detail or 'conda_build/vsnp3.yaml as declared'))
    sys.exit(0)

if leg not in legs:
    sys.exit(f"unknown leg '{leg}'. Known: {', '.join(legs)}")

# Everything the leg does not pin itself, taken from the dev environment file,
# which test_release_consistency.py holds identical to meta.yaml.
pinned = set(legs[leg])
base, deps = [], False
for line in open(env_path):
    if line.startswith('dependencies:'):
        deps = True
        continue
    if not deps:
        continue
    s = line.strip()
    if not s.startswith('- '):
        continue
    s = s[2:].split('#')[0].strip()
    if re.split(r'[ >=<]', s)[0] in pinned:
        continue
    base.append(s)

for package, pin in legs[leg].items():
    # ci.yml writes conda-style pins ('=1.26.4') and pip-style ('==1.26.4').
    # Strip the operator from those, then append '.*'.
    #
    # The '.*' is not cosmetic.  In a conda spec file a bare `python 3.12` means
    # the version 3.12 exactly -- conda zero-pads, so it excludes 3.12.13 and
    # resolves to 3.12.0.  That is the same trap that pinned this project to
    # numpy 1.26.0 for two releases, and without the wildcard it would quietly
    # reappear here.  `python 3.12.*` is what `python=3.12` means on the CLI.
    if re.match(r'^=?=[\d]', pin):
        pin = pin.lstrip('=')
    if pin and not re.search(r'[<>,*]', pin):
        pin += '.*'
    base.append(f'{package} {pin}' if pin else package)

print('\n'.join(base))
PY
}

case "${1:-}" in
    --list)
        read_specs --list
        exit 0
        ;;
    --print-specs)
        [ $# -ge 2 ] || { echo "usage: $0 --print-specs <leg>" >&2; exit 1; }
        read_specs "$2"
        exit 0
        ;;
    '' | -h | --help)
        sed -n '3,24p' "$0" | sed 's/^# \{0,1\}//'
        echo
        echo "Legs:"
        read_specs --list | sed 's/^/  /'
        exit 0
        ;;
esac

leg=$1
name=${2:-vsnp3_$leg}
shift $(( $# > 1 ? 2 : 1 ))
[ "${1:-}" = '--' ] && shift

specs=$(mktemp "${TMPDIR:-/tmp}/vsnp3_specs.XXXXXX")
trap 'rm -f "$specs"' EXIT
read_specs "$leg" > "$specs"

solver=$(command -v mamba || command -v conda) \
    || { echo 'neither mamba nor conda is on PATH' >&2; exit 1; }

echo "leg     $leg"
echo "env     $name"
echo "solver  $solver"
echo "specs   $(wc -l < "$specs" | tr -d ' ')"
grep -E '^(python|numpy|pandas) ' "$specs" | sed 's/^/          /'
echo

"$solver" create -y -n "$name" \
    --override-channels -c conda-forge -c bioconda \
    --file "$specs" "$@"

cat <<EOF

Built. Run the suites in it with:

    conda activate $name
    python3 tests/test_smoke.py
    python3 tests/test_numpy_contract.py
    python3 tests/test_annotation.py
EOF
