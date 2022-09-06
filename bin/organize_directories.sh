#!/usr/bin/env bash

mv ./dir*/* .; rm -r ./dir*/
mkdir stats; cp ./*/*stats.xlsx stats; cd stats; vsnp3_excel_merge_files.py
mv com*xlsx ..;
cd ..
rm -r ./stats/
rm slurm-*
ls
pwd