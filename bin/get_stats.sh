#!/usr/bin/env bash
#
# Merge the per-sample *stats.xlsx files below the current directory into one
# combined workbook, left in the current directory.
#
# This was semicolon-chained:
#   mkdir stats; cp ./*/*stats.xlsx stats; cd stats; vsnp3_excel_merge_files.py
#   mv com*xlsx ..; cd ..; rm -r ./stats/
# so every step ran whether or not the previous one succeeded.  If `cd stats`
# failed the script carried on, `cd ..` then moved ABOVE the intended directory,
# and `rm -r ./stats/` deleted a sibling run's stats directory.
set -euo pipefail

shopt -s nullglob
stats_files=(./*/*stats.xlsx)
if [ ${#stats_files[@]} -eq 0 ]; then
    echo "No */*stats.xlsx found below $(pwd)" >&2
    exit 1
fi

# mktemp -d rather than a fixed name, so a leftover directory from an interrupted
# run is not picked up and two concurrent runs cannot collide.
workdir=$(mktemp -d "${PWD}/.vsnp3_stats.XXXXXX")
trap 'rm -rf "$workdir"' EXIT

cp "${stats_files[@]}" "$workdir"/
( cd "$workdir" && vsnp3_excel_merge_files.py )

merged=("$workdir"/com*xlsx)
if [ ${#merged[@]} -eq 0 ]; then
    echo "vsnp3_excel_merge_files.py produced no combined workbook" >&2
    exit 1
fi
mv "${merged[@]}" "$PWD"/

printf '\nMerged %d stats file(s) into:\n' "${#stats_files[@]}"
for f in "${merged[@]}"; do
    printf '\t%s\n' "$(basename "$f")"
done
printf '\n\t%s\n\n' "$PWD"
