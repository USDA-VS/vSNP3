#!/bin/bash

# cat << EOF
# This is not a batch file, but calls batch files.

# Working directory must contain only FASTQ files.  

# Usage:
# If samples need to be forced provide options WITHOUT quotes
# vsnp3_step1_cluster.sh -t Mycobacterium_AF2122
# vsnp3_step1_cluster.sh -t NC_045512_wuhan-hu-1

# If using file options must give full file paths because vsnp3_step1_cluster.sh places files into new subdirector which causes relative paths to break.
#  vsnp3_step1_cluster.sh -f /home/tstuber/analysis/zone/vsnp3_noref_testing/dependencies/*fasta -b /home/tstuber/analysis/zone/vsnp3_noref_testing/dependencies/*gbk
# EOF

printf "\nOptions provided:\n\t" 
echo $@
echo ""

for i in *.fastq.gz; do 
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    # echo "n is : $n"
    mkdir -p $n; mv $i $n/
done

#run batch script on each sample directory, break runs into groups of 5000 so HPC doesn't overload its nodes.
COUNTER=0
root_dir=`pwd`
for folder in ./*/; do
    COUNTER=$[$COUNTER +1]
    printf "`basename $folder`\t\t"
    cd $folder
    sbatch ~/git/gitlab/vsnp3/bin/vsnp3_step1_single_sample.sh $@
    cd $root_dir
done

printf "\n$COUNTER jobs submitted\tntasks=8, 48/8=6\n"