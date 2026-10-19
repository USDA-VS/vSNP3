#!/bin/bash
#
# Run kSNP3 on every FASTA in the working directory and draw the trees.
#
# kSNP3 is no longer distributed.  docs/instructions/additional_tools.md describes
# kSNP4, whose programs are named kSNP4, MakeKSNP4infile and Kchooser4 and take
# different arguments, so this script is kept for sites that still have kSNP3 on
# PATH and is not a kSNP4 driver.  None of what it calls is part of the vsnp3 conda
# package; install separately: kSNP3, MakeKSNP3infile, MakeFasta, Kchooser,
# nw_display (newick-utils) and, for the PDF of each tree, inkscape.

starttime=`date +%s`
dir_name=${PWD##*/}
alias pause='read -p "$LINENO Enter"'

help () { 
if [ $# -eq 0 ]; then
    printf "\t\t\n!!! ERROR !!!\n"
    printf "\tIncorrect or no arguments given\n"
    printf "\tksnp3_run.sh -h for help\n"
    printf "\nUsage:\n"
    printf "\t-h = help\n"
    printf "\t-k = kmer size\n"
    printf "\t-m = minimum fraction of genomes with locus\n"
    printf "\t-a = annotate genomes\n"
    printf "\t-v = make vcfs\n\n"
    printf "ksnp3_run.sh with no arguments will run:\n"
    printf "\tksnp3_run.sh -k Kchooser_optimum -m 0.8\n\n"
    printf "may provide -k and -m arguments such as:\n"
    printf "\tksnp3_run.sh -k 21 -m 0.95\n\n"
fi
exit 1
}

# Cores for kSNP3: leave five for the node, never fewer than one (the bare
# `nproc - 5` went to 0 or negative on small machines).  nproc is GNU coreutils,
# so fall back to getconf where it is absent.
ncpu=$(nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)
NR_CPUS=$(( ncpu > 5 ? ncpu - 5 : 1 ))

# Newick tree -> SVG -> PDF.  Inkscape 1.0 removed the 0.92 `-f in -A out`
# options this script was written with; use whichever syntax the installed
# version speaks, and leave the SVG behind with a note when inkscape is absent.
svg_to_pdf() {  # svg_to_pdf <in.svg> <out.pdf>
    if ! command -v inkscape >/dev/null 2>&1; then
        echo "inkscape not found: $1 written, no PDF made"
        return 0
    fi
    if inkscape --version 2>/dev/null | grep -q 'Inkscape 0\.'; then
        inkscape -f "$1" -A "$2"
    else
        inkscape --export-type=pdf --export-filename="$2" "$1"
    fi
}

# Set flags
hflag=
kmer=
majority=
annotate=
vcf=
while getopts ':hk:m:av' OPTION; do
     case $OPTION in
        h) hflag=1
        ;;
        k) kmer=$OPTARG
        ;;
        m) majority=$OPTARG
        ;;
        a) annotate=1
        ;;
        v) vcf=1
        ;;
        ?) help
        ;;
    esac
done
shift $(($OPTIND - 1))

##

if [ "$hflag" ]; then
    help
fi

if [ "$annotate" ]; then
    printf "\n\tannotating\n\n"
fi

if [ "$vcf" ]; then
    printf "\n\tmake vcf files\n\n"
fi

#### START ####

for i in *.fa; do
     if [[ $i == *.fa ]]; then    
        echo "***.fa files found and changed to .fasta***"    
         mv "$i" "${i%.fa}.fasta"
    fi
done

for i in *.fna; do
    if [[ $i == *.fna ]]; then
        echo "***.fna files found and changed to .fasta***" 
         mv "$i" "${i%.fna}.fasta"
     fi
done

for i in *.fas; do
    if [[ $i == *.fas ]]; then
        echo "***.fas files found and changed to .fasta***" 
         mv "$i" "${i%.fas}.fasta"
     fi
done

all_files=$(ls * | grep -v "slurm.*.out" | wc -l)
fasta_files=$(ls *fasta | wc -l)
echo "all_files: $all_files"
echo "fasta_files: $fasta_files"

if [[ "$all_files" != "$fasta_files" ]]; then
    printf "\n\n ### ERROR: CHECK EXTENSIONS\n\n"
    printf "### Only .fa, .fna, .fas and .fasta allowed\n"
    exit 1
fi

for i in *.fasta; do
    name=${i%.fasta}
    new_name=$(echo $name |  sed 's/[ :,%^&*().]/_/g' | sed 's/__/_/g' | sed 's/__/_/g')
    printf "${i}\t${new_name}.fasta\n"
    mv "${i}" "${new_name}.fasta"
done

if [ -z "$kmer" ]; then
    printf "\nNo kmer supplied will run Kchooser\n"
else
    printf "\nkmer supplied: $kmer\n"
fi

if [ -z "$majority" ]; then
    majority="0.8"
    printf "Default majority will be used: $majority\n\n"
else
    printf "majority supplied: $majority\n\n"
fi

mkdir starting_files
mv *.fasta ./starting_files

# Make the kSNP needed list of file paths, pointing kSNP to files to use and there names
# A is for automatic, use S for semi-automatic and to change file names
MakeKSNP3infile starting_files/ in_list A

# Check if kmer was given
# If no kmer, then run Kchooser
if [ -z "$kmer" ]; then
    MakeFasta in_list my_fasta_master.fasta
    Kchooser my_fasta_master.fasta
    kmer=`grep "optimum value" Kchooser.report | sed 's/.*\([0-9][0-9]\).*/\1/'`
fi

echo "cpus: $NR_CPUS"
echo "kmer: $kmer"
echo "majority: $majority"

### Run kSNP3
if [ "$annotate" -a "$vcf" ]; then
    # Make a list of files that have GI numbers and can be used for annonation
    egrep -l "chromosome|complete" ./starting_files/* | sed 's/.\/starting_files\///' | sed 's/\..*//' > annotated_genomes
    kSNP3 -in in_list -outdir run -CPU ${NR_CPUS} -k ${kmer} -annotate annotated_genomes -vcf -ML -core -min_frac ${majority} | tee log.txt
elif [ "$annotate" ]; then
    # Make a list of files that have GI numbers and can be used for annonation
    #grep -l ">*gi" ./starting_files/* | sed 's/.\/starting_files\///' | sed 's/\..*//' > annotated_genomes
    # For 3.1 update (2017-09-08) change "gi" grep to "complete genome"
    egrep -l "chromosome|complete" ./starting_files/* | sed 's/.\/starting_files\///' | sed 's/\..*//' > annotated_genomes
    echo "############################ ANNOTATE #################"
    #read -p "$LINENO Enter"
    kSNP3 -in in_list -outdir run -CPU ${NR_CPUS} -k ${kmer} -annotate annotated_genomes -ML -core -min_frac ${majority} | tee log.txt
elif [ "$vcf" ]; then
    kSNP3 -in in_list -outdir run -CPU ${NR_CPUS} -k ${kmer} -vcf -ML -core -min_frac ${majority} | tee log.txt
else
    kSNP3 -in in_list -outdir run -CPU ${NR_CPUS} -k ${kmer} -ML -core -min_frac ${majority} | tee log.txt
    # kSNP3 -in in_list -outdir run -CPU ${NR_CPUS} -k ${kmer} -min_frac ${majority} | tee log.txt
fi

cp ./run/tree_tipAlleleCounts.ML.tre ./${dir_name}_tipAlleleCounts-ML.tre
cp ./run/tree_tipAlleleCounts.majority*.tre ./${dir_name}_tipAlleleCounts-majority-${majority}.tre

cat ${dir_name}_tipAlleleCounts-ML.tre | nw_display -s -w 1000 -v 20 -b 'opacity:0' -i 'font-size:8' -l 'font-family:serif;font-style:italic' -d 'stroke-width:2;stroke:blue' - > ${dir_name}-ML.svg && svg_to_pdf ${dir_name}-ML.svg ${dir_name}-ML.pdf

cat ${dir_name}_tipAlleleCounts-majority-${majority}.tre | nw_display -s -w 1000 -v 20 -b 'opacity:0' -i 'font-size:8' -l 'font-family:serif;font-style:italic' -d 'stroke-width:2;stroke:blue' - > ${dir_name}-majority-${majority}.svg && svg_to_pdf ${dir_name}-majority-${majority}.svg ${dir_name}-majority-${majority}.pdf

echo "" > mytempfile
echo "" >> mytempfile
echo "kmer: $kmer" >> mytempfile
echo "" >> mytempfile
echo "COUNT_SNPs" >> mytempfile
cat ./run/COUNT_SNPs >> mytempfile
echo "" >> mytempfile
echo "COUNT_coreSNPs" >> mytempfile
cat ./run/COUNT_coreSNPs >> mytempfile
echo "" >> mytempfile
echo "COUNT_Homoplastic_SNPs.core" >> mytempfile
cat ./run/COUNT_Homoplastic_SNPs.core >> mytempfile
echo "" >> mytempfile
echo "COUNT_Homoplastic_SNPs.majority${majority}" >> mytempfile
cat ./run/COUNT_Homoplastic_SNPs.majority${majority} >> mytempfile
echo "" >> mytempfile
if [ "$kmer" ]; then
    echo "From kSNP3 user manual:" >> mytempfile
    echo "FraCtion of Kmers (FCK) that are present in all genomes. FCK is a measure of sequence diversity, the lower is FCK the more diverse are the sequences. As diversity increases the efficiency with which kSNP3 detects SNPs decreases" >> mytempfile
    echo "Experience has shown when FCK is ≥ 0.05 SNP detection efficiency is adequate, and the accuracy of parsimony trees estimated by kSNP3 is > 97%; i.e. the trees can be considered to be reliable" >> mytempfile
    grep "FCK" Kchooser.report >> mytempfile
    echo "" >> mytempfile
fi

# remove the largest files, likely not needed
rm ./run/core_SNPs
rm ./run/nonCore_SNPs
rm ./run/SNPs_all

# created 2016-07-21 by stuber
