#!/usr/bin/env python

__version__ = "3.25"

import os
import re
import argparse
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class Kernel_Plot():
    def __init__(self, FASTA=None, sample_name=None, xlim_low=None, xlim_high=None, bin_max=None, bin_increment=None, color=None, histogram=False, debug=False):
        FASTA = os.path.basename(FASTA)
        if sample_name:
            sample_name = sample_name
        else:
            sample_name = re.sub('[_.].*', '', FASTA)

        run_set = f'snp-dists {FASTA} > distances.tab'
        os.system(run_set)

        df = pd.read_csv("distances.tab", delimiter="\t")
        os.remove("distances.tab")
        df = df.iloc[:, 1:] # drop the sample names from the first column
        df = df.iloc[:-1] #remove last row, which is root.
        df = df.drop(columns=df.columns[-1]) #remove last column, which is root.
        a = np.array(df)
        n = a.shape[0]
        b = a[np.tril_indices(n)] 
        c = np.array(b)
        if bin_max:
            bin_max = bin_max
        else:
            bin_max = c.max()
        bins = list(range(0, bin_max, bin_increment))
        if histogram:
            plt.hist(c, bins=bins)
            plt.title("histogram")
            plt.savefig(f'{sample_name}_histogram.png')
        else:
            plt.xlim(xlim_low, (bin_max + 150))
            if xlim_high:
                plt.xlim(xlim_low, xlim_high)

            ax = sns.kdeplot(c, fill=False, color='black', linewidth=1)
            ax.set(xlabel='SNP distance', ylabel='Density')
            kdeline = ax.lines[0]
            median = np.median(c)
            min_snp = np.min(c)
            max_snp = np.max(c)
            xs = kdeline.get_xdata()
            ys = kdeline.get_ydata()
            yheight = np.max(ys)
            ax.vlines(median, 0, yheight, color='red', ls=':', label=f'Median: {median}')
            ax.fill_between(xs, 0, ys, facecolor='blue', alpha=0.3)
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
    parser.add_argument('-c', '--color', action='store', dest='color', default='orange', help='Color of the kernel density plot')
    parser.add_argument('-g', '--histogram', action='store_true', dest='histogram', default=False, help='Generate histogram instead of kernel density plot')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='Keep temp files for debugging')

    args = parser.parse_args()

    kernel_plot = Kernel_Plot(FASTA=args.fasta, sample_name=args.sample_name, xlim_low=args.xlim_low, xlim_high=args.xlim_high, bin_max=args.bin_max, bin_increment=args.bin_increment, color=args.color, histogram=args.histogram, debug=args.debug)

    # #Latex report
    # latex_report = Latex_Report(myclass.sample_name)
    # myclass.latex(latex_report.tex)
    # latex_report.latex_ending()

    # #Excel Stats
    # excel_stats = Excel_Stats(myclass.sample_name)
    # myclass.excel(excel_stats.excel_dict)
    # excel_stats.post_excel()

    # temp_dir = './temp'
    # if not os.path.exists(temp_dir):
    #     os.makedirs(temp_dir)
    # files_grab = []
    # for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
    #     files_grab.extend(glob.glob(files))
    # for each in files_grab:
    #     shutil.move(each, temp_dir)

    # if args.debug is False:
    #     shutil.rmtree(temp_dir)

# Created 2021 by Tod Stuber