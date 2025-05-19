#!/usr/bin/env python

__version__ = "3.29"

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

# Force 'C' locale for consistent decimal point handling
os.environ["LC_ALL"] = "C"
locale.setlocale(locale.LC_ALL, "C")


class Kernel_Plot():
    def __init__(self, FASTA=None, sample_name=None, xlim_low=None, xlim_high=None, bin_max=None, bin_increment=None, color=None, histogram=False, debug=False):
        FASTA = os.path.basename(FASTA)
        if sample_name:
            sample_name = sample_name
        else:
            sample_name = re.sub('[_.].*', '', FASTA)

        # Use subprocess instead of os.system for better security and error handling
        run_set = ["snp-dists", FASTA]
        try:
            result = subprocess.run(run_set, capture_output=True, text=True, check=True)
            with open("distances.tab", "w") as f:
                f.write(result.stdout)
        except subprocess.CalledProcessError as e:
            print(f"Error running snp-dists: {e}")
            return

        # Use more modern pandas approaches
        df = pd.read_csv("distances.tab", delimiter="\t")
        if not debug:
            os.remove("distances.tab")
        
        # Drop the sample names from the first column
        df = df.iloc[:, 1:] 
        # Remove last row, which is root
        df = df.iloc[:-1] 
        # Remove last column, which is root
        df = df.drop(columns=df.columns[-1]) 
        
        a = np.array(df)
        n = a.shape[0]
        b = a[np.tril_indices(n)] 
        c = np.array(b)
        
        if bin_max:
            bin_max = bin_max
        else:
            bin_max = int(c.max())
            
        bins = list(range(0, bin_max, int(bin_increment)))
        
        if histogram:
            plt.hist(c, bins=bins)
            plt.title("histogram")
            plt.savefig(f'{sample_name}_histogram.png')
        else:
            plt.xlim(int(xlim_low), (bin_max + 150))
            if xlim_high:
                plt.xlim(int(xlim_low), int(xlim_high))

            ax = sns.kdeplot(c, fill=False, color='black', linewidth=1)
            ax.set(xlabel='SNP distance', ylabel='Density')
            kdeline = ax.lines[0]
            median = np.median(c)
            min_snp = np.min(c)
            max_snp = np.max(c)
            xs = kdeline.get_xdata()
            ys = kdeline.get_ydata()
            yheight = np.max(ys)
            
            # Use color parameter if provided
            fill_color = color if color else 'blue'
            
            ax.vlines(median, 0, yheight, color='red', ls=':', label=f'Median: {median}')
            ax.fill_between(xs, 0, ys, facecolor=fill_color, alpha=0.3)
            ax.legend()
            
            sns.set(style="whitegrid", palette="pastel")
            ax.vlines(median, 0, yheight, color='red', ls='-', linewidth=2)
            legend_text = f"Min: {c.min():.0f}\nMax: {c.max():.0f}\nMedian: {median:.0f}"
            ax.legend([legend_text], loc="upper right", frameon=True, fancybox=True, shadow=True, handlelength=0, handletextpad=0)
            ax.set_xlabel('SNP distance', fontsize=12)
            ax.set_ylabel('Density', fontsize=12)
            ax.set_title(f'{sample_name} SNP density plot (n={n})', fontsize=14)
            ax.grid(which='both', axis='both', linewidth=0)

            plt.savefig(f'{sample_name}_kernel.png')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\
    ---------------------------------------------------------
    Use alignment files in FASTA format output from vSNP step2.
    kernel_plots.py *.fasta
    '''), epilog='''---------------------------------------------------------''')

    parser.add_argument('-f', '--fasta', action='store', dest='fasta', help='Alignment FASTA')
    parser.add_argument('-n', '--sample_name', action='store', dest='sample_name', help='Force a sample name')
    parser.add_argument('-xh', '--xlim_low', action='store', dest='xlim_low', default=0, help='Lower limit of the x-axis')
    parser.add_argument('-yh', '--xlim_high', action='store', dest='xlim_high', default=None, help='Upper limit of the x-axis')
    parser.add_argument('-bm', '--bin_max', action='store', dest='bin_max', default=None, help='Max y length')
    parser.add_argument('-bi', '--bin_increment', action='store', dest='bin_increment', default=100, help='Bin increment for histogram')
    parser.add_argument('-c', '--color', action='store', dest='color', default='blue', help='Color of the kernel density plot')
    parser.add_argument('-g', '--histogram', action='store_true', dest='histogram', default=False, help='Generate histogram instead of kernel density plot')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Keep temp files for debugging')

    args = parser.parse_args()

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

# Created 2021 by Tod Stuber