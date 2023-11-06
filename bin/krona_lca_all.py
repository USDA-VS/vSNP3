#!/usr/bin/env python

import os
import sys
import pandas as pd
import argparse
import textwrap

import locale
# Set the locale to United States
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

def force_tax_number(kraken_output):
    df = pd.read_csv(kraken_output, sep='\t', header=None)
    target_rows = df[(df[2] == 0) & (df[0] == 'C')]
    mapping = target_rows.iloc[:,-1]
    mapping = mapping.to_frame()
    mapping = mapping.replace(' |:|', '')

    total_read_count = df.count()[0]
    classified_no_taxid_count = target_rows.count()[0]
    unclassified_read_count = df[df[0] == 'U'].count()[0]
    print("--> Total read count: {:,}, Classified reads without taxid: {:,}, {:.1%}, Unclassified reads: {:,}, {:.1%}\n\n" .format(total_read_count, classified_no_taxid_count, classified_no_taxid_count/total_read_count, unclassified_read_count, unclassified_read_count/total_read_count))

    corrected_dict = {}
    for ir, row in mapping.itertuples():
        accumulate = {}
        row = row.replace(' |:|', '')
        split_row = row.replace(' |:|', '').split()
        for lca in split_row:
            taxid, count = lca.split(":")
            accumulate.setdefault(taxid, [])
            accumulate[taxid].append(count)
            highest_taxid = '0'
            highest_value = 0
            for taxid, value_list in accumulate.items():
                value_sum = sum(list(map(int, value_list)))
                if value_sum > highest_value and taxid != '0':
                    highest_taxid = taxid
                    highest_value = value_sum
        corrected_dict[ir] = highest_taxid

    corrected_df = pd.DataFrame.from_dict(corrected_dict, orient='index')
    corrected_df = corrected_df.rename(columns={0: 2})
    df.update(corrected_df)

    krona_input = df[[1, 2]]
    krona_input.to_csv("kronaInput.txt", sep='\t', index=False, header=False)

    updated_df_count = df[(df[2] == 0) & (df[0] == 'C')].count()[0]
    print("Classified reads without taxid after update: {}\n" .format(updated_df_count))

if __name__ == "__main__": # execute if directly access by the interpreter
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

        ---------------------------------------------------------
        2018-08-17
        For use with Kraken2
        
        Input: Kraken ouput file
        Output: Tab delimited two column text file to be used by ktImportTaxonomy

        Purpose:
        Force a tax number for all classified reads, "C"

        Usage: krona_lca_all.py -f kraken2_output_file.tex
        Output will be: kronaInput.txt
        Two column file will contain read header and taxid

        '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--file_in', action='store', dest='kraken_output',  required=True, help='Required: Kraken output file (the big one)')

    args = parser.parse_args()
    print ("\nSET ARGUMENTS: ")
    print (args)
    print("")
    force_tax_number(args.kraken_output)
