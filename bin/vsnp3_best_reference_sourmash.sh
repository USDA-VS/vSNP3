#!/bin/bash

#SBATCH --ntasks=48 # Number of CPUs to allow, 48 will use full resources of one node.
#SBATCH --job-name="vSNP3s1"
#SBATCH --export=NONE

#usage: vsnp3_best_reference_sourmash.sh -r1 *.fastq.gz

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


/home/tstuber/git/gitlab/vsnp3/bin/vsnp3_best_reference_sourmash.py $@

