#!/bin/bash

# SEARCH ALL VCF FILES IN WORKING DIRECTORY.  PRINT NAME IF ALL ABSOLUTE POSITIONS FOUND IN VCF FILE.

# Check if search pairs are provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 'chrom1:pos1, chrom2:pos2, ...'"
    exit 1
fi

# Store the search pairs
SEARCH_PAIRS="$1"

# Convert the search pairs to a format suitable for grep
GREP_PATTERNS=$(echo "$SEARCH_PAIRS" | sed 's/, */\n/g' | sed 's/:/\t/' | awk '{print $1"\t"$2"[[:space:]]"}')

# Count the number of search pairs
PAIR_COUNT=$(echo "$SEARCH_PAIRS" | sed 's/, */\n/g' | wc -l)

# Iterate through all VCF files in the current directory
for vcf_file in *.vcf; do
    # Check if file exists (in case no .vcf files are found)
    [ -e "$vcf_file" ] || continue
    
    # Search for all patterns in the current file
    MATCH_COUNT=$(grep -c -f <(echo "$GREP_PATTERNS") "$vcf_file")
    
    # If all patterns are found, print the filename
    if [ "$MATCH_COUNT" -eq "$PAIR_COUNT" ]; then
        echo "$vcf_file"
    fi
done