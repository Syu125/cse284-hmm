# Data Folder Structure

This folder contains all data files needed for the HMM local ancestry inference project. It uses a consolidated script for downloading and preparing data.

## Folder Structure

```
data/
├── prepare_data.sh          # Consolidated download & prep script
├── raw/                     # Original downloaded files (do not modify)
│   ├── vcf/                 # Variant Call Format files
│   ├── panels/              # Sample population assignments
│   └── maps/                # Genetic position mappings
├── processed/               # Preprocessed/sliced data ready for analysis
└── cache/                   # Cached frequency files (auto-generated)
```

## Quick Start

### Option 1: Download Slice (Quick Testing)
```bash
bash prepare_data.sh slice    # ~50MB, 5-10 minutes
```

### Option 2: Download Full Dataset (Production)
```bash
bash prepare_data.sh full     # ~1GB, 20-30 minutes
```

### Option 3: Download Both
```bash
bash prepare_data.sh all      # Downloads both slice and full
```

## Files

### `raw/vcf/` - VCF Files
Variant Call Format files containing genotypes.

**Full Dataset**:
- `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` - All samples, full chromosome 22
- `reference_yri_ceu_chr22.vcf.gz` - YRI + CEU reference samples
- `query_asw_chr22.vcf.gz` - ASW query (admixed) samples

**Indexed** (`.vcf.gz.tbi` files) - Index for fast random access

### `raw/panels/` - Sample Information
Population assignments and sample metadata.

- `integrated_call_samples_v3.20130502.ALL.panel` - 1000 Genomes official panel
- `asw_samples.txt`, `yri_ceu_samples.txt` - Sample lists
- `rfmix_sample_map.txt` - RFMix-compatible sample mapping

### `raw/maps/` - Genetic Maps
Physical to genetic position mappings.

- `genetic_map_GRCh37_chr22.txt` - Standard genetic map (primary)
- `genetic_map_chr22_rfmix.txt` - RFMix-formatted genetic map

### `processed/` - Ready-to-Use Data
Sliced and preprocessed versions of raw data.

- `chr22_slice.vcf*` - Fast testing subset (~50MB)
- `reference_yri_ceu_slice.vcf.gz` - Sliced reference samples
- `query_asw_slice.vcf.gz` - Sliced query samples
- `rfmix_slice_output.*` - RFMix comparison outputs

### `scripts/` - Data Utilities
Bash and Python scripts for data downloading and processing.

- `download_data_full22.sh` - Download full chromosome 22
- `download_data_slice.sh` - Download small slice
- `prepare_rfmix_vcfs.sh` - Prepare files for RFMix comparison
- `run_rfmix_slice.sh` - Run RFMix on slice
- `clean_genetic_map.py` - Clean genetic map format
- `generate_rfmix_map.py` - Generate RFMix-compatible maps

### `cache/` - Auto-Generated Files
Cached computation results (safe to delete and regenerate).

- `yri_freqs_chr22.pkl` - YRI allele frequencies (auto-generated)
- `ceu_freqs_chr22.pkl` - CEU allele frequencies (auto-generated)

**Note**: These are generated on first run; subsequent runs use cached versions for speed.

## Important Notes

### Do Not Modify
- All files in `raw/` are downloaded from 1000 Genomes and should not be modified
- Modifying these will invalidate results

### Downloaded Files (~1-2 GB total)
- VCF files are large (up to 1 GB) and not tracked by Git
- If you don't have them, run the download scripts

### Folder Sizes
- `raw/vcf/` full: ~1 GB
- `raw/vcf/` slice: ~50 MB
- `raw/panels/` + `raw/maps/`: ~10 MB
- `cache/`: ~50 MB (auto-generated)
- `processed/`: ~500 MB (sliced versions)

## For More Information

See [../docs/DATA_GUIDE.md](../docs/DATA_GUIDE.md) for detailed file descriptions and usage.
