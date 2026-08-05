#!/usr/bin/env bash
#
# Run step 2 on a throwaway copy of a VCF set and snapshot the result.
#
#   tests/regression/run_step2.sh <label> <src_vcf_dir> <reference_type> [step2 args...]
#
# Writes to $VSNP3_REGRESSION_DIR (default: a temp directory), leaving
#   <workdir>/runs/<label>/       the run itself, including the copied inputs
#   <workdir>/snapshots/<label>/  the normalized snapshot for comparison
#
# Compare two snapshots with:
#   tests/regression/compare.py <snapshots>/a <snapshots>/b
#
# The copy is not a convenience.  Step 2 deletes the VCF files it reads once they
# are archived, and part of what this harness checks is that it only ever does that
# to files it owns -- so it must never be pointed at a pristine dataset.
set -uo pipefail

if [ $# -lt 3 ]; then
    sed -n '3,12p' "$0" | sed 's/^# \{0,1\}//'
    exit 1
fi

label=$1
src=$2
reference=$3
shift 3

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
workdir=${VSNP3_REGRESSION_DIR:-$(mktemp -d "${TMPDIR:-/tmp}/vsnp3_regression.XXXXXX")}

run="$workdir/runs/$label"
rm -rf "$run"
mkdir -p "$run/vcfs" "$workdir/snapshots"

shopt -s nullglob
vcfs=("$src"/*.vcf)
if [ ${#vcfs[@]} -eq 0 ]; then
    echo "No *.vcf found in $src" >&2
    exit 1
fi
cp "${vcfs[@]}" "$run/vcfs/"

# Record the inputs so a later check can prove step 2 left them alone.
( cd "$run/vcfs" && find . -name '*.vcf' -exec shasum -a 256 {} + | sort -k2 \
    > "$run/inputs_before.sha256" )

start=$(date +%s)
( cd "$run/vcfs" && python3 "$repo/bin/vsnp3_step2.py" -wd . -t "$reference" "$@" ) \
    > "$run/stdout.txt" 2>&1
rc=$?
end=$(date +%s)

# find, not `ls *.vcf`: with nullglob an unmatched glob leaves ls with no operands,
# so it lists the current directory and reports a plausible but meaningless count.
remaining=$(find "$run/vcfs" -maxdepth 1 -name '*.vcf' | wc -l | tr -d ' ')

{
    echo "label=$label"
    echo "exit=$rc"
    echo "elapsed=$((end - start))s"
    echo "source_vcfs_remaining=$remaining"
} | tee "$run/result.txt"

python3 "$here/snapshot.py" "$run/vcfs" "$workdir/snapshots/$label"
echo "workdir=$workdir"
