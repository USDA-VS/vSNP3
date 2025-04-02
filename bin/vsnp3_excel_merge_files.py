#!/usr/bin/env python3

"""
Excel File Merger

This script combines multiple Excel files in a directory into a single consolidated Excel file.
The merged file will be placed in the same directory as the input files.

Usage:
    python excel_merger.py -i <input_dir> [options]
"""

__version__ = "3.28"

import glob
import os
import sys
import logging
import argparse
from datetime import datetime
from typing import List, Optional

import pandas as pd
from pathlib import Path


class ExcelMerger:
    def __init__(
        self,
        input_dir: str = ".",
        file_pattern: str = "*.xlsx",
        index_col: Optional[str] = "sample",
        log_level: str = "INFO"
    ):
        """
        Initialize the Excel merger with configuration parameters.

        Args:
            input_dir: Directory containing input Excel files
            file_pattern: Glob pattern for matching Excel files
            index_col: Column to use as index (None if no index column)
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        """
        self.input_dir = Path(input_dir)
        self.file_pattern = file_pattern
        self.index_col = index_col

        # Set up logging
        logging.basicConfig(
            level=getattr(logging, log_level.upper()),
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger(__name__)

    def get_excel_files(self) -> List[Path]:
        """Get list of Excel files matching the pattern."""
        files = list(self.input_dir.glob(self.file_pattern))
        if not files:
            raise FileNotFoundError(
                f"No Excel files found matching pattern '{self.file_pattern}' "
                f"in directory '{self.input_dir}'"
            )
        return files

    def read_excel_file(self, file_path: Path) -> pd.DataFrame:
        """
        Read a single Excel file with error handling.

        Args:
            file_path: Path to Excel file

        Returns:
            DataFrame containing the Excel file contents
        """
        try:
            return pd.read_excel(file_path, index_col=self.index_col)
        except Exception as e:
            self.logger.error(f"Error reading file {file_path}: {str(e)}")
            raise

    def merge_files(self) -> pd.DataFrame:
        """
        Merge all Excel files in the input directory maintaining column order.
        """
        files = self.get_excel_files()
        self.logger.info(f"Found {len(files)} Excel files to merge")

        df_list = []
        column_orders = []  # Store original column orders
        all_columns = set()

        # First pass: collect columns and their orders
        for file_path in files:
            try:
                df = self.read_excel_file(file_path)
                column_orders.append(list(df.columns))
                all_columns.update(df.columns)
                df_list.append(df)
                self.logger.info(f"Successfully read {file_path}")
            except Exception as e:
                self.logger.error(
                    f"Skipping {file_path} due to error: {str(e)}")
                continue

        if not df_list:
            raise ValueError("No valid Excel files were read")

        # Find reference order (from file with most columns)
        reference_order = max(column_orders, key=len)
        final_order = list(reference_order)

        # Place missing columns near related columns
        for col in all_columns:
            if col not in final_order:
                # Find first file that has this column
                for order in column_orders:
                    if col in order:
                        # Get index of column in original file
                        idx = order.index(col)
                        # Try to find a neighboring column that exists in final_order
                        for neighbor_idx in range(idx, -1, -1):
                            if neighbor_idx < len(order) and order[neighbor_idx] in final_order:
                                insert_idx = final_order.index(
                                    order[neighbor_idx]) + 1
                                final_order.insert(insert_idx, col)
                                break
                        else:
                            # If no neighbor found, append to end
                            final_order.append(col)
                        break

        # Reindex all dataframes with complete column list
        for i, df in enumerate(df_list):
            df_list[i] = df.reindex(columns=final_order)

        return pd.concat(df_list, sort=False)

    def save_merged_file(self, df: pd.DataFrame) -> Path:
        """
        Save the merged DataFrame to a new Excel file in the input directory.

        Args:
            df: DataFrame to save

        Returns:
            Path to the saved file
        """
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        output_file = self.input_dir / f'combined_excel_{timestamp}.xlsx'

        try:
            df.to_excel(output_file)
            self.logger.info(
                f"Successfully saved merged file to {output_file}")
            return output_file
        except Exception as e:
            self.logger.error(f"Error saving merged file: {str(e)}")
            raise

    def run(self) -> Path:
        """
        Execute the full merge process.

        Returns:
            Path to the output file
        """
        try:
            merged_df = self.merge_files()
            return self.save_merged_file(merged_df)
        except Exception as e:
            self.logger.error(f"Merge process failed: {str(e)}")
            raise


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge multiple Excel files into a single file in the same directory",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-i", "--input-dir",
        default=".",
        help="Directory containing input Excel files"
    )

    parser.add_argument(
        "-p", "--pattern",
        default="*.xlsx",
        help="File pattern to match Excel files (e.g., '*.xlsx', 'data_*.xls')"
    )

    parser.add_argument(
        "--index-col",
        default="sample",
        help="Column to use as index (use 'None' for no index)"
    )

    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the logging level"
    )

    return parser.parse_args()


def main():
    """Main entry point for the script."""
    try:
        # Parse command line arguments
        args = parse_arguments()

        # Convert string 'None' to actual None
        index_col = None if args.index_col.lower() == 'none' else args.index_col

        # Initialize and run merger
        merger = ExcelMerger(
            input_dir=args.input_dir,
            file_pattern=args.pattern,
            index_col=index_col,
            log_level=args.log_level
        )

        output_file = merger.run()
        print(
            f"Successfully merged Excel files. Output saved to: {output_file}")
        sys.exit(0)

    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
