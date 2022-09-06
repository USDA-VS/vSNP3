# Additional Programs

Many programs can be used to help identify reads.  Three programs that can be useful to use alongside vSNP are Mashtree, kSNP and Kraken.  

[Mashtree](https://github.com/lskatz/mashtree) and [kSNP](https://pubmed.ncbi.nlm.nih.gov/25913206) are reference independent phylogenetic tree building programs.  Mashtree is very fast, kSNP is slower but results may be more accurate and additional information is provide to help qualify results.

[Kraken](https://ccb.jhu.edu/software/kraken2/) uses kmers to identify reads.  If a sample is not behaving as expected or contamination is suspected Kraken is a powerful tool for determining read identification.  When used with Krona an easy to read, and share, HTML file is provided.

# Example Dataset

## FASTAs for reference free tree building

```
cd ~; mkdir tree_test; cd tree_test
```

Make `list` with the following

```
NC_000962
NC_018143
NZ_CP017594
NZ_OW052188
NC_015758
NC_002945
NZ_CP039850
NZ_LR882497
```

Download list

```
for i in `cat list`; do vsnp3_download_fasta_gbk_gff_by_acc.py -a $i -f; done
```

# Mashtree

New conda environment

```
conda create --name mashtree
```

```
conda activate mashtree
```

```
conda install mashtree -c conda-forge -c bioconda 
```

Change to directory with test files

```
cd ~/tree_test
```
Build tree from FASTAs

```
mashtree --sketch-size 1000000 --numcpus 4 *.fasta > mashtree.tre
```


# kSNP

New conda environment

```
conda create --name ksnp
```

```
conda activate ksnp
```

```
conda install -c hcc -c conda-forge -c bioconda ksnp fasttree # Only works from Linux environment (WSL works)
```
Build tree from FASTAs
```
cd ~/tree_test
```
Just an Example.  Supply your specific path to wrapper.
```
~/git/gitlab/vsnp3/bin/vsnp3_ksnp_wrapper.sh
```

# Kraken/Krona

Make new conda environment

```
conda create --name kraken
```

```
conda activate kraken
```

```
conda install kraken2 krona krakentools sra-tools=2.11.0 pandas -c conda-forge -c bioconda 
```

This asking to do additional setup for these programs.

[Download](https://benlangmead.github.io/aws-indexes/k2) Kraken database.

There are many Databases to choose from.  If unsure and download speeds allow try the standard database.  If a smaller database is necessary Standard-8 may be a good option.
```
cd ~; wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20220607.tar.gz
```
```
mkdir k2_standard_08gb; tar -xzf k2_standard_08gb_20220607.tar.gz -C k2_standard_08gb
```

After installation of Kraken and Krona the following message is shown.

```
Krona installed.  You still need to manually update the taxonomy
databases before Krona can generate taxonomic reports.  The update
script is ktUpdateTaxonomy.sh.  The default location for storing
taxonomic databases is /Users/tstuber/opt/anaconda3/envs/more/opt/krona/taxonomy

If you would like the taxonomic data stored elsewhere, simply replace
this directory with a symlink.  For example:

rm -rf /Users/tstuber/opt/anaconda3/envs/more/opt/krona/taxonomy
mkdir /path/on/big/disk/taxonomy
ln -s /path/on/big/disk/taxonomy /Users/tstuber/opt/anaconda3/envs/more/opt/krona/taxonomy
ktUpdateTaxonomy.sh
```
The default location should work unless disk space is an issue

```
ktUpdateTaxonomy.sh
```

## FASTQ for Kraken Testing

```
cd; mkdir kraken_test
```
```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6046640/SRR6046640 #fortuitum
```
```
fastq-dump --split-files SRR6046640
```
```
rm SRR6046640; mv SRR6046640_1.fastq SRR6046640_R1.fastq; mv SRR6046640_2.fastq SRR6046640_R2.fastq; pigz *fastq
```
Just an Example.  Supply your specific path to wrapper.
```
~/git/gitlab/vsnp3/bin/vsnp3_kraken2_wrapper.py -r1 SRR6046640_R1.fastq.gz -r2 SRR6046640_R2.fastq.gz --database ~/k2_standard_08gb
```