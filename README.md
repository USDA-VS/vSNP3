# vSNP3: Variant Calling and SNP Analysis Tool

vSNP3 is a robust tool for high-resolution SNP analysis, designed for disease tracing and outbreak investigations. It generates BAM, VCF, and annotated SNP tables along with corresponding phylogenetic trees.

## Features

- Generates annotated SNP tables and phylogenetic trees
- Processes large-scale datasets
- Accommodates multiple references
- Accreditation-friendly with easy error correction and SNP validation
- Scalable reporting
- Compatible with Python 3.8 - 3.11

## Installation

### Prerequisites

- Anaconda or Miniconda

### Conda Installation

```bash
conda create -c conda-forge -c bioconda -n vsnp3 vsnp3=3.25
```

For detailed Anaconda setup instructions, see [conda instructions](./docs/instructions/conda_instructions.md).

### Verification

To verify the installation:

```bash
which vsnp3_step1.py
vsnp3_step1.py -h
vsnp3_step2.py -h
```

## Quick Start

1. Clone the test dataset:
   ```bash
   cd ~
   git clone https://github.com/USDA-VS/vsnp3_test_dataset.git
   ```

2. Add reference:
   ```bash
   cd ~/vsnp3_test_dataset/vsnp_dependencies
   vsnp3_path_adder.py -d `pwd`
   ```

3. Run test with AF2122 (Mycobacterium bovis):
   - Step 1:
     ```bash
     cd ~/vsnp3_test_dataset/AF2122_test_files/step1
     vsnp3_step1.py -r1 *_R1*.fastq.gz -r2 *_R2*.fastq.gz -t Mycobacterium_AF2122
     ```
   - Step 2:
     ```bash
     cd ~/vsnp3_test_dataset/AF2122_test_files/step2
     vsnp3_step2.py -a -t Mycobacterium_AF2122
     ```

## Usage

vSNP3 is divided into two main steps:

### Step 1

Main entry: `vsnp3_step1.py`

Additional scripts:
- vsnp3_alignment_vcf.py
- vsnp3_assembly.py
- vsnp3_best_reference_sourmash.py
- vsnp3_fastq_stats_seqkit.py
- vsnp3_group_reporter.py
- vsnp3_vcf_annotation.py
- vsnp3_zero_coverage.py

### Step 2

Main entry: `vsnp3_step2.py`

Additional scripts:
- vsnp3_fasta_to_snps_table.py
- vsnp3_group_on_defining_snps.py
- vsnp3_html_step2_summary.py
- vsnp3_remove_from_analysis.py

### Utility Scripts

- vsnp3_path_adder.py
- vsnp3_bruc_mlst.py
- vsnp3_download_fasta_gbk_gff_by_acc.py
- vsnp3_excel_merge_files.py
- vsnp3_filter_finder.py
- vsnp3_spoligotype.py

For detailed usage of each script, use the `-h` option.


## Version Enhancements

# vSNP3 Version Enhancements

1. Extended Python support for versions 3.9 to 3.11.

2. Improved modularity for enhanced function usage.

3. Implemented Sourmash for optimal reference selection.

4. Flexibility in dependency file provision:
   - File options can be explicit or based on reference directory.

5. Added customizable thresholds as user options.

6. Introduced utility script for downloading complete or partial GenBank files.

8. Updated and enhanced annotation descriptions:
   - Improved handling of cases with no provided annotation.

9. Refined indel calling:
   - Indels at group SNP positions are now called as 'N'.

10. Improved error handling:
    - Script halts when VCF files in a set have incorrect chromosome information.

11. Enhanced reproducibility:
    - Defining SNPs and filters are now zipped with starting VCF files.

12. Added capability to test defining SNPs individually.

13. Enabled generation of cascading SNP tables from FASTA alignments alone.

14. Improved version tracking and option capturing.

15. Enhanced file input options for greater flexibility.

16. Added beta support for Oxford Nanopore reads.

17. Refined naming conventions:
    - 'all_vcf' now used as the defining SNP column label.

18. Simplified installation process:
    - Removed Picard, eliminating Java requirements and facilitating easier Conda installation.

19. Implemented strict sample/VCF file naming requirements for metadata matching:
    - Two-column metadata Excel worksheet's first column must exactly match VCF filename, with fallback options.

20. Improved organization of cascading SNP tables.

21. Modified default behavior for unmapped reads:
    - Assembly of unmapped reads now requires a specific flag.

22. Enhanced logging:
    - Detailed steps are now recorded in 'run_log.txt'.

23. Changed spoligotype handling for TB complex isolates:
    - Now requires '-s' option with step 1, or use of 'vsnp3_spoligotype.py'.

24. Relocated Brucella MLST functionality:
    - Removed from step 1; now available through 'vsnp3_bruc_mlst.py'.

- And [more...](https://github.com/USDA-VS/vSNP3/releases)


## Additional Tools

For information on additional tools, see [Additional Tools](./docs/instructions/additional_tools.md).

## Contributing

[Citation](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-024-10437-5)

For more information or support, please open an issue on the GitHub repository.
