# vSNP3 Parabricks GPU Acceleration Implementation

This document describes the GPU acceleration implementation for vSNP3 using Parabricks, created in the `parabricks-gpu-acceleration` branch.

## Overview

The GPU acceleration implementation replaces CPU-intensive BWA alignment and FreeBayes variant calling with Parabricks GPU-accelerated equivalents:

- **BWA mem** → **Parabricks fq2bam** (includes alignment + sorting + duplicate marking)
- **FreeBayes** → **Parabricks HaplotypeCaller** (GATK-based variant calling)

## Key Benefits

- **Speed**: Significant acceleration for alignment and variant calling on large genomes
- **Quality**: Proven GATK HaplotypeCaller for SNP calling
- **Compatibility**: Same output files and formats as CPU version
- **Flexibility**: Optional GPU mode - can fall back to CPU processing

## Environment Setup

### Step 1: Create Conda Environment

```bash
# Create conda environment with vsnp3 (includes all dependencies automatically)
mamba create -n vsnp3_parabricks_3.33 -c conda-forge -c bioconda vsnp3

# This automatically installs: python, bwa, minimap2, samtools, bcftools,
# freebayes, pandas, numpy, biopython, openpyxl, seqkit, sourmash, spades,
# and all other vsnp3 dependencies

# Copy dependencies from existing vsnp3 environment
cp /project/shared/miniconda3/envs/vsnp3_3.25/dependencies \
   /project/shared/miniconda3/envs/vsnp3_parabricks_3.33/dependencies
```

### Step 2: Setup Parabricks Container (One-Time)

**Choose one option:**

**Option A: Shared Container (Recommended for multi-user HPC)**
```bash
# Setup shared location (admin or first user)
sudo mkdir -p /project/shared/containers
cd /project/shared/containers
module load apptainer/1.1.9
apptainer pull docker://nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1
```

**Option B: Personal Container**
```bash
# Setup in your home directory
mkdir -p ~/containers
cd ~/containers
module load apptainer/1.1.9
apptainer pull docker://nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1
```

### Step 3: Test Container

```bash
# Test container works (adjust path as needed)
apptainer exec --nv /project/shared/containers/clara-parabricks_4.3.1-1.sif pbrun version

# Test on GPU node
srun --partition=prod-compute-gpu --gres=gpu:1 --pty bash
module load apptainer/1.1.9
apptainer exec --nv /project/shared/containers/clara-parabricks_4.3.1-1.sif nvidia-smi
```

## Usage

After setup, using GPU acceleration is simple:

### Slurm Integration (Recommended)

```bash
# Navigate to your data directory
cd /project/mycobacteria_brucella/my_samples/

# Submit GPU job - script automatically finds container and adds --gpu flag
sbatch /home/tstuber/git/gitlab/vsnp3/internal/git_gpu_vsnp3_step1.sh -r reference.fasta
```

### Direct Command Line (Optional)

```bash
# Activate environment (use conda activate, not mamba activate)
conda activate vsnp3_parabricks_3.33

# Run with --gpu flag (container auto-detected)
vsnp3_step1.py -r1 sample_R1.fastq.gz -r2 sample_R2.fastq.gz --gpu -r reference.fasta
```

### How It Works

**The GPU batch script automatically:**
- ✅ **Finds containers** in common locations:
  - `/project/shared/containers` (recommended)
  - `$HOME/containers`
  - Current working directory
  - Other standard locations
- ✅ **Tests container functionality** before use
- ✅ **Adds `--gpu` flag** to vsnp3_step1.py automatically
- ✅ **Provides clear instructions** if container missing

**No manual container copying or path management needed!**

## Technical Implementation

### Files Modified

1. **vsnp3_step1.py**
   - Added `--gpu` argument
   - Pass GPU parameter to Alignment class

2. **vsnp3_alignment_vcf.py**
   - Added GPU processing logic in Alignment class
   - Parabricks fq2bam for alignment pipeline
   - Parabricks HaplotypeCaller for variant calling
   - Compatibility with existing output formats

3. **New Slurm Script: git_gpu_vsnp3_step1.sh**
   - Optimized for GPU nodes (partition: prod-compute-gpu)
   - Reduced CPU allocation (12 CPUs vs 48)
   - GPU resource allocation (#SBATCH --gres=gpu:1)
   - Uses vsnp3_parabricks_3.33 environment

### Processing Flow Comparison

#### CPU Processing (Traditional)
1. BWA index creation
2. BWA mem alignment → SAM
3. samtools fixmate → BAM
4. samtools sort → sorted BAM
5. samtools markdup → final BAM
6. samtools index
7. FreeBayes variant calling → VCF
8. Quality filtering (QUAL>20)

#### GPU Processing (New)
1. BWA index creation (for compatibility)
2. **Parabricks fq2bam** → final BAM (includes alignment, sorting, duplicate marking)
3. samtools index
4. **Parabricks HaplotypeCaller** → VCF
5. Quality filtering (QUAL>20)

### Compatibility Considerations

- **Nanopore reads**: GPU acceleration disabled for nanopore (uses minimap2 + bcftools as before)
- **Output files**: Identical structure and naming to CPU version
- **Statistics**: Fake markduplicate_stats.txt created for compatibility
- **VCF format**: GATK HaplotypeCaller produces GATK-compatible VCF (no MQM→MQ conversion needed)
- **Coverage analysis**: Uses same samtools depth + zero coverage workflow

### Resource Requirements

#### GPU Node Requirements
- **Partition**: prod-compute-gpu
- **GPUs**: 1 GPU per job
- **CPUs**: 12 CPUs (reduced from 48)
- **Memory**: Standard allocation
- **Node**: Must support CUDA-compatible GPU

#### Software Requirements
- **Parabricks**: clara-parabricks-pipelines package
- **CUDA**: Compatible GPU drivers
- **Traditional tools**: Maintained for fallback and compatibility

## Performance Considerations

### When to Use GPU

**Recommended for**:
- Large eukaryotic genomes (originally designed for bacterial/viral)
- High-throughput processing
- Illumina paired-end or single-end reads
- Workloads requiring speed optimization

**Not recommended for**:
- Nanopore reads (uses CPU pathway automatically)
- Very small genomes where GPU overhead exceeds benefits
- Limited GPU availability

### Expected Performance

- **Alignment**: Significant speedup for large reference genomes
- **Variant Calling**: GATK HaplotypeCaller generally faster than FreeBayes
- **Overall**: Best performance gains on larger genomes and higher coverage

## Validation

### Output Verification

The GPU implementation produces identical output structure:
```
alignment_<reference_name>/
├── <sample>_nodup.bam
├── <sample>_nodup.bam.bai
├── <sample>_zc.vcf
├── <sample>_annotated.vcf (if GBK provided)
├── <sample>_run_log.txt
├── <sample>_stats.xlsx
├── <sample>.pdf
└── unmapped_reads/
```

### Quality Assurance

- GATK HaplotypeCaller has proven accuracy for bacterial variant calling
- Same quality filtering thresholds (QUAL>20)
- Identical downstream processing (coverage analysis, annotation)

## Troubleshooting

### Common Issues

1. **Container not found**:
   ```bash
   # Check container exists and is accessible
   ls -la *parabricks*.sif *clara*.sif
   apptainer exec --nv container.sif pbrun version
   ```

2. **GPU not available**:
   ```bash
   # Check GPU access on compute node
   srun --partition=prod-compute-gpu --gres=gpu:1 --pty bash
   nvidia-smi
   ```

3. **Container fails to run**:
   ```bash
   # Check Apptainer module is loaded
   module list | grep apptainer
   # Test simple container command
   apptainer exec container.sif echo "Hello World"
   ```

4. **Network restrictions during container pull**:
   ```bash
   # Try from login node (not compute node)
   # Or download manually and transfer to HPC
   ```

5. **Memory/storage issues**:
   ```bash
   # Container is ~15GB, check available space
   df -h .
   # Reduce concurrent jobs per GPU node
   ```

### Container Troubleshooting

**Container Pull Issues:**
```bash
# If docker pull fails, try alternatives:
apptainer pull docker://nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1  # different version
apptainer pull docker://nvcr.io/nvidia/clara/clara-parabricks:latest
```

**GPU Access Issues:**
```bash
# Verify GPU is available to container
apptainer exec --nv container.sif nvidia-smi
# If GPU not accessible, check:
# 1. --nv flag is used
# 2. GPU drivers are compatible
# 3. Running on GPU-enabled node
```

**File Permission Issues:**
```bash
# Ensure container file has correct permissions
chmod +r container.sif
# Ensure working directory is accessible
```

### Debugging

- Use `--debug` flag to retain intermediate files
- Check `<sample>_run_log.txt` for detailed processing log including container commands
- Verify GPU availability: `nvidia-smi`
- Test container separately: `apptainer exec --nv container.sif pbrun version`
- Check Slurm job output for container-specific errors

## Future Enhancements

### Potential Improvements

1. **Automatic GPU detection**: Auto-enable GPU if available
2. **Hybrid processing**: CPU for small genomes, GPU for large
3. **Multi-GPU support**: Scale across multiple GPUs
4. **Benchmarking**: Systematic performance comparison

### Compatibility

- Designed to eventually merge back to main branch
- Maintains backward compatibility with all existing workflows
- Optional feature - does not affect standard CPU processing

---

## Summary

This implementation provides significant performance improvements for vSNP3 step1 processing while maintaining complete compatibility with existing workflows. The GPU acceleration is particularly beneficial for processing large genomes and high-throughput scenarios, making vSNP3 suitable for eukaryotic genome analysis while preserving its accuracy for bacterial and viral genomics.