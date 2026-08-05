#!/usr/bin/env bash
#
# Flatten the dir*/ subdirectories a cluster run leaves behind, then merge the
# per-sample stats into one workbook.
#
# The previous version was:
#   mv ./dir*/* .; rm -r ./dir*/
#   mkdir stats; cp ./*/*stats.xlsx stats; cd stats; vsnp3_excel_merge_files.py
#   mv com*xlsx ..; cd ..; rm -r ./stats/; rm slurm-*
# Nothing was conditional, so a partly failed `mv` was followed unconditionally by
# `rm -r ./dir*/`, destroying whatever had not been moved.  The `cd ..` then
# `rm -r ./stats/` had the same wrong-directory hazard as get_stats.sh.
set -euo pipefail

shopt -s nullglob

dirs=(./dir*/)
if [ ${#dirs[@]} -eq 0 ]; then
    echo "No dir*/ subdirectories in $(pwd); nothing to flatten." >&2
else
    for d in "${dirs[@]}"; do
        contents=("$d"*)
        if [ ${#contents[@]} -gt 0 ]; then
            # Only remove the directory once its contents have actually moved.
            mv "${contents[@]}" .
        fi
        rmdir "$d"
    done
    printf 'Flattened %d directory/ies.\n' "${#dirs[@]}"
fi

# Same merge as get_stats.sh; call it rather than repeating it.
if command -v get_stats.sh > /dev/null 2>&1; then
    get_stats.sh
else
    "$(dirname "$(readlink -f "$0" 2>/dev/null || echo "$0")")/get_stats.sh"
fi

slurm_logs=(slurm-*)
if [ ${#slurm_logs[@]} -gt 0 ]; then
    rm -- "${slurm_logs[@]}"
    printf 'Removed %d slurm log(s).\n' "${#slurm_logs[@]}"
fi
