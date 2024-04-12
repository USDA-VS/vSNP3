# Additional Programs

Many programs can be used to help identify reads.  Three programs useful to use alongside vSNP are Mashtree, kSNP and Kraken.  

Best results from vSNP are provided when a sample is less than 1,000 SNPs from a reference.  If a sample is too distant from a reference the alignment error can cause time consuming corrections.  Good reference selection is important for best results.  Mashtree and kSNP can help in reference selection.  [Mashtree](https://github.com/lskatz/mashtree) and [kSNP](https://pubmed.ncbi.nlm.nih.gov/25913206) are reference independent phylogenetic tree building programs.  Mashtree is very fast, kSNP is slower but results may be more accurate and additional information is provide to help qualify results.  

[Kraken](https://ccb.jhu.edu/software/kraken2/) uses kmers to identify reads.  If a sample is not behaving as expected or contamination is suspected Kraken is a powerful tool for determining read identification quickly.  When used with Krona an easy to read HTML file is provided.

Below are brief installation and usage insturctions for these tools.  See their individual links for more detail.  The scripts provided for kSNP and Kraken are only for example.  Users should make updates as needed.

# Example Dataset

## FASTAs for reference-free tree building

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

vsnp3 available from github
```
cd ~; git clone https://github.com/USDA-VS/vsnp3.git
```

Building Mashtree, kSNP and Kraken in their own conda environments ensures installation dependencies do not conflict.  Scripts provided in the cloned vsnp3 repo above are needed since conda environments are independent.

# Mashtree

Create conda environment

```
conda create -n mashtree -c conda-forge -c bioconda mashtree
```
```
conda activate mashtree
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

As of late 2023 kSNP needs to be download from [sourceforge](https://sourceforge.net/projects/ksnp/files/).

There is a new version of kSNP as of 2023, kSNP4.1.

Choose the prebuild binary for your environment, download and unzip.

Place unzipped file in desired location (${HOME} will work)

Add to PATH, `PATH="${HOME}/kSNP4/kSNP4.1pkg":$PATH`

Change directory to location of FASTA files

```
MakeKSNP4infile -indir ./ -outfile myInfile S
```
```
Kchooser4 -in myInfile
```
```
kSNP4 -in myInfile -outdir run -CPU 8 -k 21 -core -ML -min_frac 0.8
```

# Kraken/Krona

Create conda environment

```
conda create -n kraken -c conda-forge -c bioconda kraken2 krona krakentools wget pandas pigz
```

After the conda install it will provide additional setup instructions for these programs.

[Download](https://benlangmead.github.io/aws-indexes/k2) Kraken database.

There are many Databases to choose from.  If unsure and download speeds allow try the standard database.  If a smaller database is necessary Standard-8 may be a good option.
```
cd ~; wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20220607.tar.gz
```
```
mkdir k2_standard_08gb; tar -xzf k2_standard_08gb_20220607.tar.gz -C k2_standard_08gb
```

If needed link database to conda environment and download taxonomy.

```
rm -rf ${HOME}/anaconda3/envs/kraken/opt/krona/taxonomy
ln -s ${HOME}/k2_standard_08gb ${HOME}/anaconda3/envs/kraken/opt/krona/taxonomy
ktUpdateTaxonomy.sh
```

Just an Example.  Supply your specific path to wrapper.
```
~/anaconda3/envs/vsnp3/bin/vsnp3_kraken2_wrapper.py -r1 SRR6046640_R1.fastq.gz -r2 SRR6046640_R2.fastq.gz --database ~/k2_standard_08gb
```

## SRA Tools
```
conda create -n sra-tools -c conda-forge -c bioconda -n sra-tools
```
```
fasterq-dump --split-files -O . SRR26282520
```
```
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6046640/SRR6046640 
```
```
fastq-dump --split-files SRR6046640
```
### macOS
```
~/sratoolkit.3.0.7-mac64/bin/fasterq-dump -S SRR6046640
```
### Docker
Download Docker.  It must be running.
```
docker pull ncbi/sra-tools
docker run -t --rm -v $PWD:/output:rw -w /output ncbi/sra-tools fasterq-dump -e 2 -p SRR6046640
```
### Singularity
```
singularity pull docker://ncbi/sra-tools
singularity run sra-tools_latest.sif fasterq-dump -e 2 -p SRR6046640
```
