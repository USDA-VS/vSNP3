#!/usr/bin/env python

__version__ = "3.35"

import os
import re
import argparse
import locale
import textwrap
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class Kernel_Plot():
    def __init__(self, FASTA=None, sample_name=None, xlim_low=None, xlim_high=None, bin_max=None, bin_increment=None, color=None, histogram=False, debug=False):
        # Validate required FASTA parameter
        if FASTA is None:
            raise ValueError("FASTA parameter is required")
        
        if not os.path.exists(FASTA):
            raise FileNotFoundError(f"FASTA file not found: {FASTA}")
        
        FASTA = os.path.basename(FASTA)
        
        # Robust sample name handling
        if sample_name:
            sample_name = str(sample_name)
        else:
            sample_name = re.sub(r'[_.].*', '', FASTA)
            if not sample_name:  # If regex resulted in empty string
                sample_name = "sample"

        # Use subprocess instead of os.system for better security and error handling
        run_set = ["snp-dists", FASTA]
        try:
            result = subprocess.run(run_set, capture_output=True, text=True, check=True)
            with open("distances.tab", "w") as f:
                f.write(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running snp-dists: {e}")
            return
        except FileNotFoundError:
            print("Error: snp-dists command not found. Please ensure it's installed and in PATH.")
            return

        try:
            # Use more modern pandas approaches
            df = pd.read_csv("distances.tab", delimiter="\t")
            if not debug:
                os.remove("distances.tab")
        except Exception as e:
            print(f"Error reading distances file: {e}")
            return
        
        # Validate dataframe has content
        if df.empty:
            print("Error: Empty distance matrix")
            return
            
        # Drop the sample names from the first column
        df = df.iloc[:, 1:] 
        # Remove last row, which is root
        df = df.iloc[:-1] 
        # Remove last column, which is root
        df = df.drop(columns=df.columns[-1]) 
        
        # Validate we still have data after processing
        if df.empty:
            print("Error: No data remaining after processing distance matrix")
            return
        
        a = np.array(df)
        n = a.shape[0]
        
        # Validate we have sufficient data
        if n < 2:
            print("Error: Insufficient samples for analysis (need at least 2)")
            return
            
        b = a[np.tril_indices(n)] 
        c = np.array(b)
        
        # Validate we have distance data
        if len(c) == 0:
            print("Error: No distance data available")
            return
        
        # Robust bin_max handling
        try:
            if bin_max is not None:
                bin_max = int(bin_max)
            else:
                bin_max = int(c.max())
        except (ValueError, TypeError):
            bin_max = int(c.max()) if len(c) > 0 else 100
            
        # Robust bin_increment handling
        try:
            bin_increment = int(bin_increment) if bin_increment is not None else 100
        except (ValueError, TypeError):
            bin_increment = 100
            
        bins = list(range(0, bin_max + bin_increment, bin_increment))
        
        # Robust xlim_low handling
        try:
            xlim_low = int(xlim_low) if xlim_low is not None else 0
        except (ValueError, TypeError):
            xlim_low = 0
        
        if histogram:
            plt.figure(figsize=(10, 6))
            plt.hist(c, bins=bins)
            plt.title("histogram")
            output_file = f"{sample_name}_histogram.png"
            plt.savefig(output_file)
            plt.close()
        else:
            try:
                plt.figure(figsize=(12, 8))
                plt.xlim(xlim_low, (bin_max + 150))
                if xlim_high is not None:
                    try:
                        xlim_high = int(xlim_high)
                        plt.xlim(xlim_low, xlim_high)
                    except (ValueError, TypeError):
                        pass  # Use default xlim if conversion fails

                ax = sns.kdeplot(c, fill=False, color='black', linewidth=1)
                ax.set(xlabel='SNP distance', ylabel='Density')
                kdeline = ax.lines[0]
                
                # Robust statistical calculations
                try:
                    median = np.median(c)
                    min_snp = np.min(c)
                    max_snp = np.max(c)
                except Exception as e:
                    print(f"Error calculating statistics: {e}")
                    return
                
                xs = kdeline.get_xdata()
                ys = kdeline.get_ydata()
                
                if len(ys) > 0:
                    yheight = np.max(ys)
                else:
                    print("Error: No KDE data generated")
                    return
                
                # Use color parameter if provided, with fallback
                fill_color = color if color else 'blue'
                
                ax.vlines(median, 0, yheight, color='red', ls=':', label=f'Median: {median:.0f}')
                ax.fill_between(xs, 0, ys, facecolor=fill_color, alpha=0.3)
                ax.legend()
                
                sns.set_style("whitegrid")
                ax.vlines(median, 0, yheight, color='red', ls='-', linewidth=2)
                legend_text = f"Min: {min_snp:.0f}\nMax: {max_snp:.0f}\nMedian: {median:.0f}"
                ax.legend([legend_text], loc="upper right", frameon=True, fancybox=True, shadow=True, handlelength=0, handletextpad=0)
                ax.set_xlabel('SNP distance', fontsize=12)
                ax.set_ylabel('Density', fontsize=12)
                ax.set_title(f"{sample_name} SNP density plot (n={n})", fontsize=14)
                ax.grid(which='both', axis='both', linewidth=0)

                output_file = f"{sample_name}_kernel.png"
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()
                
            except Exception as e:
                print(f"Error creating kernel density plot: {e}")
                return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------
    Use alignment files in FASTA format output from vSNP step2.
    kernel_plots.py *.fasta
    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', required=True, help='Alignment FASTA (required)')
    parser.add_argument('-n', '--sample_name', action='store', dest='sample_name', help='Force a sample name')
    parser.add_argument('-xh', '--xlim_low', action='store', dest='xlim_low', default=0, help='Lower limit of the x-axis')
    parser.add_argument('-yh', '--xlim_high', action='store', dest='xlim_high', default=None, help='Upper limit of the x-axis')
    parser.add_argument('-bm', '--bin_max', action='store', dest='bin_max', default=None, help='Max y length')
    parser.add_argument('-bi', '--bin_increment', action='store', dest='bin_increment', default=100, help='Bin increment for histogram')
    parser.add_argument('-c', '--color', action='store', dest='color', default='blue', help='Color of the kernel density plot')
    parser.add_argument('-g', '--histogram', action='store_true', dest='histogram', default=False, help='Generate histogram instead of kernel density plot')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Keep temp files for debugging')

    args = parser.parse_args()

    try:
        kernel_plot = Kernel_Plot(
            FASTA=args.fasta, 
            sample_name=args.sample_name, 
            xlim_low=args.xlim_low, 
            xlim_high=args.xlim_high, 
            bin_max=args.bin_max, 
            bin_increment=args.bin_increment, 
            color=args.color, 
            histogram=args.histogram, 
            debug=args.debug
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

# Created 2021 by Tod Stuber