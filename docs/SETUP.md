# Environment Setup Guide

## System Requirements

- **Windows 10+**: WSL2 (Windows Subsystem for Linux)
- **macOS/Linux**: Native terminal
- **Python**: 3.10+ (managed via conda)
- **Disk Space**: ~5GB for full data, ~500MB for slice

## Step 1: Install Miniconda (if needed)

### Windows Users
1. Install WSL2:
   ```powershell
   wsl --install
   ```
2. In WSL terminal, download and install Miniconda:
   ```bash
   curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

### macOS Users
```bash
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

### Linux Users
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

## Step 2: Create Conda Environment

Navigate to the project root and run:

```bash
cd cse284-hmm
conda env create -f environment.yml
conda activate hmm_env
```

Verify installation:
```bash
python --version  # Should be 3.10.x
python -c "import pysam; import pandas; print('All dependencies installed!')"
```

## Step 3: Download Data

### Option A: Quick Start (Recommended for first-time users)

Small slice of chromosome 22 (~50MB):
```bash
cd data
bash download_data_slice.sh
cd ..
```

### Option B: Full Dataset

Complete chromosome 22 (~1GB):
```bash
cd data
bash download_data_full22.sh
cd ..
```

**Expected files after download:**
- `data/raw/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- `data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel`
- `data/raw/maps/genetic_map_GRCh37_chr22.txt`

## Step 4: Verify Setup

Run the sanity check to confirm everything works:
```bash
cd src
python ../scripts/01_sanity_check.py
```

Expected output:
- Reference sample (NA19625) ancestry painted
- PNG file saved to `outputs/sanity_check/`
- Console output with statistics

If this succeeds, your environment is ready!

## Troubleshooting

### Issue: `pysam` module not found
```bash
conda activate hmm_env
conda install -c bioconda pysam
```

### Issue: `tabix` not found
```bash
conda install -c bioconda tabix
```

### Issue: VCF file corrupted/incomplete
```bash
# Delete and re-download
cd data
rm -f raw/vcf/*.vcf.gz*
bash download_data_slice.sh
```

### Issue: Memory error with full dataset
- Use slice version instead: `bash download_data_slice.sh`
- Or reduce sample count in analysis scripts

## Environment Variables (Optional)

For faster analysis, set cache directory:
```bash
export HMM_CACHE_DIR="$(pwd)/data/cache"
```

## Next Steps

1. Read [USAGE.md](USAGE.md) to run analyses
2. Check [DATA_GUIDE.md](DATA_GUIDE.md) for file descriptions
3. Review [IMPLEMENTATION.md](IMPLEMENTATION.md) for technical details
