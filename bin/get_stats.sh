#!/usr/bin/env bash

mkdir stats; cp ./*/*stats.xlsx stats; cd stats; vsnp3_excel_merge_files.py
mv com*xlsx ..;
cd ..
rm -r ./stats/

ls
printf "\n\t"
pwd
echo ""