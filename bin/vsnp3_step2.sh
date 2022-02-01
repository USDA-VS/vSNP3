#!/bin/bash

#SBATCH --job-name="vsnp3_#2"
#SBATCH --ntasks=48 # Number of CPUs to allow, 48 will use full resources of one node.
#SBATCH --export=NONE

# Provide entire option string WITHOUT quotes
# sbatch vsnp3_step2.sh -wd ../vcf_source --remove

module load python-3.9.6-gcc-9.2.0-3m34vgo
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
module load py-scikit-allel-1.3.5-gcc-9.2.0-wmnzhy6
module load py-svgwrite-1.4.1-gcc-9.2.0-k7fwxvu
module load py-cairosvg-2.5.2-gcc-9.2.0-g3wmkcd
module load py-cairocffi-1.2.0-gcc-9.2.0-7tw6glw
module load py-pathlib2-2.3.6-gcc-9.2.0-l7tuqbf
module load py-cssselect2-0.4.1-gcc-9.2.0-ewrk5wm
module load py-webencodings-0.5.1-gcc-9.2.0-wbw7zaq
module load py-tinycss2-1.1.0-gcc-9.2.0-adilvq3
module load py-defusedxml-0.7.1-gcc-9.2.0-v3feloi
module load py-pyyaml-5.4.1-gcc-9.2.0-v5wmagc
module load py-fsspec-2021.7.0-gcc-9.2.0-rq6fs5y
module load py-numpy-1.21.1-gcc-9.2.0-w5ejrxe
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
module load raxml-8.2.12-gcc-9.2.0-ckartoq
python --version

# provide options as:
# sbatch vsnp3_step2.sh "-a -t ASFV_Georgia_2007" # if vsnp3/bin is in PATH
# sbatch vsnp3_step2.sh "--remove"
# sbatch vsnp3_step2.sh "--wd ../vcf_source"
# sbatch ${HOME}/git/gitlab/vsnp3/bin/vsnp3_step2.sh -t NC_045512_wuhan-hu-1 -remove # when working directory contains vcf files
# sbatch ${HOME}/git/gitlab/vsnp3/bin/vsnp3_step2.sh -t NC_045512_wuhan-hu-1 -wd ../vcf_source

echo "Running vsnp3_step2.py ${@}"
${HOME}/git/gitlab/vsnp3/bin/vsnp3_step2.py $@