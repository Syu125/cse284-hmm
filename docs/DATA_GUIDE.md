# Data Guide

## Overview

This project uses genetic data from the **1000 Genomes Project Phase 3** to train and validate an HMM for local ancestry inference. The data is organized hierarchically: raw downloads, processed versions, and cached computations.

---

## Data Organization

```
data/
├── raw/                     # Downloaded from 1000 Genomes (do not modify)
│   ├── vcf/                 # Variant Call Format files
│   ├── panels/              # Sample population assignments
│   ├── maps/                # Genetic position mappings
│   └── scripts/             # Download and prep scripts
├── processed/               # Cleaned/sliced versions
├── scripts/                 # Data preparation utilities
└── cache/                   # Pickled frequency files
```

---

## Source Data Files

### 1. VCF Files (Variant Call Format)

**Location**: `data/raw/vcf/`

**File Options**:

| File | Size | Variants | Description |
|------|------|----------|-------------|
| `chr22_slice.vcf` | ~50 MB | ~20k | Subset of chromosome 22 (fast testing) |
| `ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz` | ~1 GB | ~490k | Full chromosome 22 for production |
| `reference_yri_ceu_chr22.vcf.gz` | ~200 MB | Subset | Reference populations (YRI + CEU) |
| `query_asw_chr22.vcf.gz` | ~100 MB | Subset | Query population (ASW samples) |

**Download Command**:
```bash
cd data
bash scripts/download_data_slice.sh     # Quick version (~50MB)
bash scripts/download_data_full22.sh    # Full version (~1GB)
```

**Reference Source**:
- Project: [1000 Genomes Project Phase 3](http://www.1000genomes.org/)
- FTP: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

### 2. Sample Panel File

**File**: `data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel`

**Purpose**: Maps sample IDs to populations

**Format**:
```
sample_id   pop   super_pop   gender
NA19625     ASW   AFR         male
NA19626     ASW   AFR         female
NA18507     YRI   AFR         male
NA12878     CEU   EUR         female
...
```

**Key Populations**:
- `YRI` - Yoruba in Ibadan, Nigeria (West African)
- `CEU` - Utah residents with Northern/Western European ancestry
- `ASW` - African ancestry in Southwest USA (admixed target)

**Source**: 1000 Genomes sample information file

### 3. Genetic Map Files

**Files**: 
- `data/raw/maps/genetic_map_GRCh37_chr22.txt` (primary)
- `data/raw/maps/genetic_map_chr22_rfmix.txt` (reformatted for RFMix)

**Purpose**: Maps physical positions (base pairs) to genetic positions (centiMorgans)

**Format**:
```
Physical_Position       COMBINED_rate(cM/Mb)    Genetic_Map_Position(cM)
50000                   0.0                     0.0
50002                   0.0                     0.0
50100                   0.000283333             0.0000283333
...
```

**Why it matters**:
- Recombination rate depends on genetic distance, not physical distance
- Genetic distance varies across the genome
- Critical for HMM transition model

---

## Processed Data Files

### Created During Execution

These files are **generated** by analysis scripts (not downloaded):

| File | Created By | Purpose |
|------|-----------|---------|
| `data/cache/yri_freqs_chr22.pkl` | `scripts/01_sanity_check.py` | Cached YRI allele frequencies |
| `data/cache/ceu_freqs_chr22.pkl` | `scripts/01_sanity_check.py` | Cached CEU allele frequencies |
| `outputs/sanity_check/NA19625_karyogram.png` | `01_sanity_check.py` | Ancestry painting PNG |
| `outputs/real_samples/asw_ancestry_results.csv` | `03_real_sample_analysis.py` | Results table |

---

## Sample Selection

### Slice Dataset (Quick Testing)
Used by default if full dataset isn't available.

**YRI Reference** (~50 samples):
- NA18507, NA18508, ... (Yoruba from Nigeria)

**CEU Reference** (~50 samples):
- NA12878, NA12889, ... (Utah European ancestry)

**ASW Query** (~20 samples):
- NA20509, NA20510, ... (Southwest USA admixed)

### Full Dataset (Production)
All samples from chromosome 22 (Phase 3):

**YRI**: 108 samples  
**CEU**: 99 samples  
**ASW**: 61 samples  
**TOTAL**: 1000+ samples (all populations)

---

## Data Preparation Pipeline

### Download Step
```bash
# Run this once to get data
bash data/scripts/download_data_full22.sh
```

### Processing Step (Automatic)
Scripts automatically:
1. **Extract frequencies**: Read VCF files for YRI/CEU allele counts
2. **Interpolate positions**: Map SNP positions to genetic distances
3. **Cache results**: Save frequencies to pickle files for reuse

### Manual Preprocessing (Optional)
```bash
# Slice the full VCF (if you want a smaller working copy)
bcftools view -r 22:1-5000000 input.vcf.gz -o sliced.vcf

# Recalculate genetic map (advanced)
python data/scripts/clean_genetic_map.py

# Prepare RFMix-specific files (for benchmarking)
bash data/scripts/prepare_rfmix_vcfs.sh
```

---

## Caching Strategy

### First Run (Frequency Calculation)
```
Time: 5-10 minutes (full dataset)
Output: yri_freqs_chr22.pkl, ceu_freqs_chr22.pkl
```

### Subsequent Runs
```
Time: <1 second (loads from cache)
```

### Cache Location
```bash
data/cache/
├── yri_freqs_chr22.pkl          # YRI allele frequencies
└── ceu_freqs_chr22.pkl          # CEU allele frequencies
```

### Invalidate Cache
```bash
rm data/cache/*.pkl
# Next run will recalculate
```

---

## Data Statistics

### Chromosome 22 Coverage

| Metric | Value |
|--------|-------|
| **Physical Length** | 51.3 Mb |
| **Genetic Length** | ~44 cM |
| **Variants** | ~490k (full), ~20k (slice) |
| **Average SNP spacing** | ~100 bp |
| **Recombination rate** | ~0.86 cM/Mb (average) |

### Population Allele Frequencies

Example from 1000 Genomes Phase 3:

| SNP | Position | YRI Freq | CEU Freq | ASW Freq |
|-----|----------|----------|----------|----------|
| rs1 | 50000000 | 0.45 | 0.32 | 0.38 |
| rs2 | 50001000 | 0.89 | 0.05 | 0.52 |
| rs3 | 50002000 | 0.12 | 0.78 | 0.45 |

**Interpretation**:
- SNP rs1: Fairly similar between populations
- SNP rs2: Strong YRI allele (ancestry informative)
- SNP rs3: Strong CEU allele (ancestry informative)

---

## Benchmarking Files

### RFMix Outputs
If comparing against RFMix:

```
evaluation/results/
├── rfmix_slice_output.msp.tsv       # Ancestry state calls (RFMix)
├── rfmix_slice_output.sis.tsv       # Posterior probabilities
└── rfmix_comparison_summary.csv     # Comparison metrics
```

### Expected Comparison Results

Typical concordance with RFMix:
- **Full dataset**: 92-96% agreement at SNP level
- **Slice dataset**: 90-95% agreement

Variance due to:
- Different emission models
- Different transition models
- Different parameter tuning

---

## Data Privacy & Attribution

### Citation
When publishing results using this data, cite:

> The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.

### Data License
- Public data from 1000 Genomes Project
- No restrictions for research use
- Appropriate to publish results

### Sample Identifiers
- Sample IDs are anonymized (e.g., NA19625)
- Can be linked to demographic data via 1000 Genomes database
- Respect privacy when publishing

---

## Troubleshooting

### VCF File Currupted
```bash
# Check file integrity
bcftools info data/raw/vcf/file.vcf.gz

# Re-download if corrupted
rm data/raw/vcf/file.vcf.gz*
bash data/scripts/download_data_full22.sh
```

### Missing Genetic Map
```bash
# Check if file exists
ls data/raw/maps/genetic_map*.txt

# Download if missing
cd data/scripts
bash download_data_full22.sh
```

### Cache Issues
```bash
# Clear cache and recalculate
rm -f data/cache/*.pkl

# Rebuild on next run
python ../scripts/01_sanity_check.py
```

### Frequency Calculation Taking Too Long
```bash
# Use slice dataset instead
bash data/scripts/download_data_slice.sh
```

---

## Advanced: Custom Data

To use different populations or chromosomes:

1. **Download** different VCF files from 1000 Genomes
2. **Update** sample panel file with population assignments
3. **Download** genetic map for your chromosome
4. **Modify** analysis scripts to reference new paths
