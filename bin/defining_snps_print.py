#!/usr/bin/env python3

import pandas as pd
import argparse

def main(file_path):
    # Load the Excel file without assuming a header row
    df = pd.read_excel(file_path, header=None)

    # Skip the first column
    df = df.iloc[:, 1:]

    # Extract row 1 and row 2
    row1 = df.iloc[0]
    row2 = df.iloc[1]

    # Print row 1 and row 2 as columns with row 1 as column 2 and row 2 as column 1
    for r1, r2 in zip(row1, row2):
        print(f"{r2}\t{r1}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Print the first two rows of an Excel file as columns with row 1 as column 2 and row 2 as column 1, skipping the first column.")
    parser.add_argument("-f", "--file", required=True, help="Path to the Excel file")
    
    args = parser.parse_args()
    
    main(args.file)
