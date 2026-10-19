#!/bin/bash

# SEARCH ALL VCF FILES IN WORKING DIRECTORY.  PRINT NAME IF ALL ABSOLUTE POSITIONS FOUND IN VCF FILE.
# Positions suffixed with "!" must NOT be present in the VCF file.

# Check if search pairs are provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 'chrom1:pos1, chrom2:pos2, chrom3:pos3!, ...'"
    exit 1
fi

# Store the search pairs
SEARCH_PAIRS="$1"

# One position per line.  tr and awk rather than sed: "\n" and "\t" in a sed
# replacement are GNU extensions that BSD sed (macOS, FreeBSD) does not honour.
POSITIONS=$(printf '%s\n' "$SEARCH_PAIRS" | tr ',' '\n' | awk '{gsub(/^[ \t]+|[ \t]+$/, ""); if ($0 != "") print}')

# Split into required (no !) and excluded (with !) positions
REQUIRED=$(printf '%s\n' "$POSITIONS" | grep -v '!$')
EXCLUDED=$(printf '%s\n' "$POSITIONS" | grep '!$' | sed 's/!$//')

# Build grep patterns for required positions: "chrom<TAB>pos" followed by whitespace
REQUIRE_PATTERNS=$(printf '%s\n' "$REQUIRED" | awk -F: 'NF >= 2 {print $1"\t"$2"[[:space:]]"}')
REQUIRE_COUNT=$(printf '%s\n' "$REQUIRED" | grep -c .)

# Build grep patterns for excluded positions
if [ -n "$EXCLUDED" ]; then
    EXCLUDE_PATTERNS=$(printf '%s\n' "$EXCLUDED" | awk -F: 'NF >= 2 {print $1"\t"$2"[[:space:]]"}')
fi

# Iterate through all VCF files in the current directory
for vcf_file in *.vcf; do
    # Check if file exists (in case no .vcf files are found)
    [ -e "$vcf_file" ] || continue

    # Check all required positions are present
    if [ "$REQUIRE_COUNT" -gt 0 ]; then
        MATCH_COUNT=$(grep -c -f <(printf '%s\n' "$REQUIRE_PATTERNS") "$vcf_file")
        if [ "$MATCH_COUNT" -ne "$REQUIRE_COUNT" ]; then
            continue
        fi
    fi

    # Check no excluded positions are present
    if [ -n "$EXCLUDED" ]; then
        EXCLUDE_MATCH=$(grep -c -f <(printf '%s\n' "$EXCLUDE_PATTERNS") "$vcf_file")
        if [ "$EXCLUDE_MATCH" -gt 0 ]; then
            continue
        fi
    fi

    echo "$vcf_file"
done
