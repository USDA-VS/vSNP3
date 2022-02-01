#!/bin/bash

'''
This is not a batch file, but calls batch files.

Working directory must contain only FASTQ files.  FASTQ files must all be expected to be found with sourmash and have a corresponding reference type that has been made available with vsnp3_path_adder.py.  Paired reads will be grouped on uniqueness of name based on everything left of [_.].  Therefore samples named sample_001 and sample_002 would need to be renamed sample-001 and sample-002 for correct automatic grouping.

At completion collect stats with `mkdir stats; cp **/*stats.xlsx stats; cd stats; vsnp3_excel_merge_files.py`

Usage:

If samples need to be forced provide options WITHOUT quotes
vsnp3_step1_cluster.sh -t Mycobacterium_AF2122
vsnp3_step1_cluster.sh -t NC_045512_wuhan-hu-1 --skip_assembly

vsnp3_step1_cluster.sh # if reference type has been added and select is available via sourmash

# need to give full file paths because vsnp3_step1_cluster.sh places files into new subdirector which causes relative paths to break.
 vsnp3_step1_cluster.sh -f /home/tstuber/analysis/zone/vsnp3_noref_testing/dependencies/*fasta -b /home/tstuber/analysis/zone/vsnp3_noref_testing/dependencies/*gbk --skip_assembly
'''

cpus_alloted_per_sample=8

reference=$rflag
printf "\nSpecified reference: $reference\n"

cpu_available=$(grep -c ^processor /proc/cpuinfo)
samples_per_node=$(awk -v x=$cpu_available -v y=$cpus_alloted_per_sample 'BEGIN { rounded = sprintf("%.0f", x/y); print rounded }')
# files_per_node=$(( samples_per_node*2 )) #multiply by 2 for paired reads
files_per_node=$(( samples_per_node )) #accounting for single reads, ie. decrease samples per node

printf "\nNumber of cpus:  $cpu_available\n"
printf "cpus begin used per sample:  $cpus_alloted_per_sample\n"
printf "Samples per node  $samples_per_node\n"

printf "\n\nparameters given: ${@}\n\n"

#set up dir with multiple FASTQ pairs
counter=0
file_count=$(ls *gz | wc -l)
while [ $file_count -gt 0 ]; do
    counter=$[$counter+1]
    counter_formated=$(printf "%02d" ${counter})
    mkdir dir${counter_formated}
    mv $(ls *gz | head -${files_per_node}) dir${counter_formated}
    echo "Moving $(ls dir${counter_formated}/*gz | wc -l) files to dir${counter_formated}"
    file_count=$(ls *gz | wc -l)
done
printf "\nNumber of nodes being used:  $counter\n"

#run batch script on each sample directory, break runs into groups of 30 so HPC doesn't overload its nodes.
COUNTER=0
SPLIT=30
root_dir=`pwd`
for folder in ./*/; do
    COUNTER=$[$COUNTER +1]
    echo $folder
    cd ./$folder
    if (( SPLIT > COUNTER )); then
        sbatch ~/git/gitlab/vsnp3/bin/vsnp3_step1.sh $@ #pass wild cards
    else
        sbatch ~/git/gitlab/vsnp3/bin/vsnp3_step1.sh $@
        SPLIT=$[$SPLIT +30]
        printf "\n"
        printf "At $COUNTER of ${counter_formated}\n"
        pwd
        printf "PAUSING 5 minutes, RUNNING NODES IN GROUPS OF 30, "
        date
        printf "\n"
        sleep 5m
    fi
    cd $root_dir
done
