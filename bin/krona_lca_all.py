#!/usr/bin/env python3

import os
import sys
import pandas as pd
import argparse
import textwrap


def best_taxid(lca_mapping):
    '''
    The taxid with the most k-mer hits in one kraken2 LCA mapping string, or 0.

    The mapping looks like "0:12 1763:5 |:| 0:10": "taxid:count" pairs, with
    " |:| " separating the two mates of a pair.  Taxid 0 (unclassified k-mers)
    never wins.  Returned as an int, the type of kraken2's taxid column.
    '''
    accumulate = {}
    for lca in str(lca_mapping).replace(' |:|', ' ').split():
        if ':' in lca:
            taxid, count = lca.rsplit(':', 1)
            try:
                accumulate[taxid] = accumulate.get(taxid, 0) + int(count)
            except ValueError:
                continue
    highest_taxid, highest_value = 0, 0
    for taxid, value_sum in accumulate.items():
        if value_sum > highest_value and taxid != '0':
            try:
                highest_taxid = int(taxid)
            except ValueError:
                continue
            highest_value = value_sum
    return highest_taxid


def force_tax_number(kraken_output):
    df = pd.read_csv(kraken_output, sep='\t', header=None)
    target_rows = df[(df[2] == 0) & (df[0] == 'C')]
    mapping = target_rows.iloc[:,-1]
    mapping = mapping.to_frame()
    # The mapping strings are parsed as they are.  A DataFrame.replace({' |:|': ''},
    # regex=True) used to sit here; as a regex ' |:|' means "a space, OR a colon,
    # OR the empty string", so it deleted every space and colon from every mapping
    # ("0:12 1763:5 |:| 0:10" became "01217635||010") before anything was read.

    total_read_count = len(df)
    classified_no_taxid_count = len(target_rows)
    unclassified_read_count = len(df[df[0] == 'U'])
    print("--> Total read count: {:,}, Classified reads without taxid: {:,}, {:.1%}, Unclassified reads: {:,}, {:.1%}\n\n" .format(total_read_count, classified_no_taxid_count, classified_no_taxid_count/total_read_count, unclassified_read_count, unclassified_read_count/total_read_count))

    corrected_dict = {}
    for ir, row in mapping.itertuples():
        # str.replace takes no regex= keyword: the two calls that used to be here
        # raised TypeError on the first classified read with taxid 0.
        if isinstance(row, str):
            corrected_dict[ir] = best_taxid(row)

    corrected_df = pd.DataFrame.from_dict(corrected_dict, orient='index')
    if not corrected_df.empty:  # Only rename if DataFrame is not empty
        corrected_df = corrected_df.rename(columns={0: 2})
        df.update(corrected_df)

    krona_input = df[[1, 2]]
    krona_input.to_csv("kronaInput.txt", sep='\t', index=False, header=False)

    updated_df_count = len(df[(df[2] == 0) & (df[0] == 'C')])
    print("Classified reads without taxid after update: {}\n" .format(updated_df_count))

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        2018-08-17 (Updated for Python 3.12)
        For use with Kraken2
        
        Input: Kraken ouput file
        Output: Tab delimited two column text file to be used by ktImportTaxonomy

        Purpose:
        Force a tax number for all classified reads, "C"

        Usage: krona_lca_all.py -f kraken2_output_file.txt
        Output will be: kronaInput.txt
        Two column file will contain read header and taxid

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--file_in', action='store', dest='kraken_output',  required=True, help='Required: Kraken output file (the big one)')

    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)
    print("")
    force_tax_number(args.kraken_output)