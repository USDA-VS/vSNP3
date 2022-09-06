#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task 12 # Number of CPUs to allow, 48 will use full resources of one node.
#SBATCH --job-name="vsnp3_#1"
#SBATCH --export=NONE


#This script is called by vsnp3_step1_cluster.sh.  See vsnp3_step1_cluster.sh for usage.

module load python-3.9.6-gcc-9.2.0-3m34vgo
module load bwa-0.7.17-gcc-9.2.0-75y5prf
module load py-pandas-1.3.1-gcc-9.2.0-gw265cc
module load py-pytz-2021.1-gcc-9.2.0-xb774mj
module load py-python-dateutil-2.8.2-gcc-9.2.0-fx34jn5
module load py-pysam-0.16.0.1-gcc-9.2.0-rdtuano
module load py-pyvcf-0.6.8-gcc-9.2.0-3ffhiio
module load py-humanize-3.11.0-gcc-9.2.0-dhrqxdh
module load py-biopython-1.79-gcc-9.2.0-3enl3jk
module load py-regex-2021.8.3-gcc-9.2.0-57mzr63
module load py-dask-2021.7.2-gcc-9.2.0-lse6ple
module load py-toolz-0.11.1-gcc-9.2.0-ygba23k
module load py-openpyxl-3.0.7-gcc-9.2.0-f6qf34a
module load py-jdcal-1.4.1-gcc-9.2.0-35xeiwj
module load py-et-xmlfile-1.1.0-gcc-9.2.0-xswduch
module load py-pycpuinfo-8.0.0-gcc-9.2.0-svpueek
module load py-pyparsing-2.4.7-gcc-9.2.0-qo5y2g5
module load py-xlrd-2.0.1-gcc-9.2.0-buye5dr
module load py-xlsxwriter-3.0.1-gcc-9.2.0-vzp2bql
module load py-scikit-allel-1.3.5-gcc-9.2.0-3s5cu
module load py-svgwrite-1.4.1-gcc-9.2.0-k7fwxvu
module load py-cairosvg-2.5.2-gcc-9.2.0-g3wmkcd
module load py-cairocffi-1.2.0-gcc-9.2.0-7tw6glw
module load py-pathlib2-2.3.6-gcc-9.2.0-l7tuqbf
module load py-cssselect2-0.4.1-gcc-9.2.0-ewrk5wm
module load py-webencodings-0.5.1-gcc-9.2.0-wbw7zaq
module load py-tinycss2-1.1.0-gcc-9.2.0-adilvq3
module load py-defusedxml-0.7.1-gcc-9.2.0-v3feloi
module load py-pyyaml-5.4.1-gcc-9.2.0-v5wmagc
module load texinfo-6.5-gcc-9.2.0-z6nbuiy
module load texlive-live-gcc-9.2.0-pnsn7pt

#Sourmash
module load py-sourmash-4.2.1-gcc-9.2.0-tupts2w
module load py-bz2file-0.98-gcc-9.2.0-qbcqu3f
module load py-cachetools-4.2.2-gcc-9.2.0-6wwjsi3
module load py-cffi-1.14.6-gcc-9.2.0-tej6pbg
module load py-cycler-0.10.0-gcc-9.2.0-wcvk3yj
module load py-deprecation-2.1.0-gcc-9.2.0-dmtbwga
module load py-kiwisolver-1.3.1-gcc-9.2.0-m5vzi3a
module load py-matplotlib-3.4.2-gcc-9.2.0-k2s4vt5
module load py-numpy-1.21.1-gcc-9.2.0-w5ejrxe
module load py-packaging-21.0-gcc-9.2.0-o7mjnay
module load py-pillow-8.3.1-gcc-9.2.0-nbvtthz
module load py-pycparser-2.20-gcc-9.2.0-x2tmq7u
module load py-pyparsing-2.4.7-gcc-9.2.0-qo5y2g5
module load py-python-dateutil-2.8.2-gcc-9.2.0-fx34jn5
module load py-scipy-1.7.1-gcc-9.2.0-5hw7kez
module load py-screed-1.0.5-gcc-9.2.0-pp2r3ea
module load py-setuptools-57.4.0-gcc-9.2.0-sxsblui
module load py-six-1.16.0-gcc-9.2.0-qqrbhzs

 # module load py-fsspec-2021.7.0-gcc-9.2.0-rq6fs5y #step 2

#System calls
module load abyss-2.2.3-gcc-9.2.0-waljb6i
module load bcftools-1.10.2-gcc-9.2.0-caq2tcp
module load bzip2-1.0.8-gcc-9.2.0-4q7b53i
module load freebayes-1.3.2-gcc-9.2.0-z36gxwk
module load gcc-9.2.0-gcc-9.1.0-gltu7dh
module load minimap2-2.17-gcc-9.2.0-rtdj3zn
module load parallel-20190222-gcc-9.2.0-r4rpmin
module load pigz-2.4-gcc-9.2.0-ybomdyx
module load raxml-8.2.11-gcc-9.2.0-e4o4fpc
module load samtools-1.10-gcc-9.2.0-fqappr7
module load seqkit-0.16.0-gcc-9.2.0-jawexb7
module load spades-3.13.0-gcc-9.2.0-2dgoqsd
module load vcffilter-1.0-gcc-9.2.0-zubzwdu
module load vcflib-1.0.1-gcc-9.2.0-3g4sbyw
module load vcftools-0.1.14-gcc-9.2.0-wp6mi4x
module load vsnp-3.09-gcc-10.3.0-bhtp7s3
python --version

# export PATH="${HOME}/git/gitlab/vsnp3/bin:$PATH" #test from gitlab and not loaded module
printf "\nNode:\n\t$SLURMD_NODENAME\\nn"
printf "Scripts call from:\n\t"
which vsnp3_path_adder.py
printf "\nCWD:\n\t"
pwd
printf "\n\nOptions:\n\t"
echo $@
echo ""

file_count=$(ls *gz | wc -l)
if [[ $file_count == 1 ]]; then
    echo "single read"
    vsnp3_step1.py -r1 *.fastq.gz $@ &
else        
    echo "paired reads"
    vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz $@ &
fi
wait