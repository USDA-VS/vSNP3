# Additional Bioinformatics Tools for Genomic Analysis

## Table of Contents
1. [Introduction](#introduction)
2. [Example Dataset](#example-dataset)
3. [Mashtree](#mashtree)
4. [kSNP](#ksnp)
5. [Kraken/Krona](#krakenkreona)
6. [SRA Tools](#sra-tools)

## Introduction

In genomic analysis, particularly when working with vSNP (variant calling and phylogenetic analysis tool), several complementary programs can significantly enhance your workflow. This guide focuses on three powerful tools: Mashtree, kSNP, and Kraken, along with instructions for using SRA Tools to obtain sequence data.

### Why use these tools?

- **Reference Selection**: vSNP performs best when samples are within 1,000 SNPs of a reference. Mashtree and kSNP can aid in selecting appropriate references.
- **Phylogenetic Analysis**: Both Mashtree and kSNP build reference-independent phylogenetic trees, offering different trade-offs between speed and accuracy.
- **Read Identification**: Kraken excels at rapid read identification, crucial for detecting contamination or unexpected sample composition.

## Example Dataset

Before we dive into the tools, let's set up an example dataset to work with.

### Preparing FASTA files for reference-free tree building

```bash
# Create and navigate to a working directory
cd ~
mkdir tree_test && cd tree_test

# Create a list of accession numbers
cat << EOF > accession_list.txt
NC_000962
NC_018143
NZ_CP017594
NZ_OW052188
NC_015758
NC_002945
NZ_CP039850
NZ_LR882497
EOF

# Download FASTA files using vSNP3
while read i; do
    vsnp3_download_fasta_gbk_gff_by_acc.py -a $i -f
done < accession_list.txt
```

Note: Ensure you have vSNP3 installed. If not, you can install it following these [instructions](https://github.com/USDA-VS/vSNP3).

## Mashtree

Mashtree is a rapid method for creating phylogenetic trees based on MinHash distances.

### Installation and Usage

```bash
# Create and activate a conda environment for Mashtree
conda create -n mashtree -c conda-forge -c bioconda mashtree
conda activate mashtree

# Navigate to the directory with test files
cd ~/tree_test

# Build a tree from FASTA files
mashtree --sketch-size 1000000 --numcpus 4 *.fasta > mashtree.tre
```

## kSNP

kSNP is a SNP-based approach to phylogenetic tree construction that doesn't require genome alignment or a reference genome.

### Installation

As of late 2023, kSNP4.1 needs to be downloaded from [SourceForge](https://sourceforge.net/projects/ksnp/files/).

1. Download the prebuilt binary for your environment.
2. Unzip the file and place it in your desired location (e.g., `${HOME}`).
3. Add kSNP to your PATH:
   ```bash
   echo 'export PATH="${HOME}/kSNP4/kSNP4.1pkg:$PATH"' >> ~/.zshrc
   source ~/.zshrc
   ```

### Usage

```bash
# Navigate to the directory with FASTA files
cd ~/tree_test

# Prepare input file
MakeKSNP4infile -indir ./ -outfile myInfile S

# Choose optimal k-mer size
Kchooser4 -in myInfile

# Run kSNP
kSNP4 -in myInfile -outdir ksnp_run -CPU 8 -k 21 -core -ML -min_frac 0.8
```

## Kraken/Krona

Kraken is a system for ultrafast metagenomic sequence classification using exact k-mer matches. Krona provides interactive visualization of the results.

### Installation

```bash
# Create and activate a conda environment for Kraken
conda create -n kraken -c conda-forge -c bioconda kraken2 krona krakentools wget pandas pigz
conda activate kraken

# Download Kraken database (example using standard-8 database)
cd ~
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz
mkdir k2_standard_08gb
tar -xzf k2_standard_08gb_*.tar.gz -C k2_standard_08gb

# Link database and update taxonomy (adjust paths as needed)
rm -rf ${HOME}/anaconda3/envs/kraken/opt/krona/taxonomy
ln -s ${HOME}/k2_standard_08gb ${HOME}/anaconda3/envs/kraken/opt/krona/taxonomy
ktUpdateTaxonomy.sh
```

Additional prebuilt Kraken Databases available [here](https://benlangmead.github.io/aws-indexes/k2)

### Usage

Here's an example using a wrapper script (adjust the path to your specific location):

```bash
./vsnp3/bin/vsnp3_kraken2_wrapper.py -r1 SRR6046640_R1.fastq.gz -r2 SRR6046640_R2.fastq.gz --database ~/k2_standard_08gb
```

## SRA Tools

SRA Tools allow you to access data from the NCBI Sequence Read Archive.

### Installation

```bash
conda create -n sra-tools -c conda-forge -c bioconda sra-tools
conda activate sra-tools
```

### Usage

#### Basic Usage

```bash
# Download and split FASTQ files
fasterq-dump --split-files -O . SRR26282520

# Alternative method
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6046640/SRR6046640
fastq-dump --split-files SRR6046640
```

#### Platform-Specific Instructions

##### macOS

If you've downloaded the SRA Toolkit directly:

```bash
~/sratoolkit.3.0.7-mac64/bin/fasterq-dump -S SRR6046640
```

##### Docker

Ensure Docker is installed and running, then:

```bash
docker pull ncbi/sra-tools
docker run -t --rm -v $PWD:/output:rw -w /output ncbi/sra-tools fasterq-dump -e 2 -p SRR6046640
```

##### Singularity

```bash
singularity pull docker://ncbi/sra-tools
singularity run sra-tools_latest.sif fasterq-dump -e 2 -p SRR6046640
```

## Conclusion

These tools form a powerful suite for genomic analysis, complementing vSNP3 and each other. By mastering Mashtree, kSNP, Kraken/Krona, and SRA Tools, you'll be well-equipped to handle a wide range of genomic analysis tasks efficiently.

Remember to always check for the latest versions and updates of these tools, as bioinformatics software evolves rapidly.

For more detailed information on each tool, please refer to their respective documentation:

- [Mashtree GitHub](https://github.com/lskatz/mashtree)
- [kSNP Documentation](https://sourceforge.net/projects/ksnp/files/)
- [Kraken2 Manual](https://github.com/DerrickWood/kraken2/wiki/Manual)
- [SRA Tools Documentation](https://github.com/ncbi/sra-tools/wiki)

Happy analyzing!
