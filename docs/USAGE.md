# Running Analyses

This guide explains how to run each analysis in sequence. Start from the root directory (`cse284-hmm/`).

## Quick Start

All scripts should be run from the project root directory (`cse284-hmm/`).

### Step 0: Download Data (if not already done)

```bash
bash data/prepare_data.sh slice    # Download ~50MB slice (quick)
# OR
bash data/prepare_data.sh full     # Download ~1GB full dataset
```

### Step 1-4: Run Analyses

---

## Analysis Workflows

### 1. Sanity Check (Quick Test)
**Purpose**: Validate the HMM pipeline on a single known sample  
**Sample**: NA19625 (YRI/CEU admixed)  
**Runtime**: 5-10 minutes  
**Data**: Works with both slice and full datasets

```bash
cd src
python ../scripts/01_sanity_check.py
```

**Expected Output**:
- `outputs/sanity_check/NA19625_karyogram.png` - Ancestry painting for the sample
- Console output with ancestry statistics

**What it does**:
- Loads reference populations (YRI, CEU)
- Computes allele frequencies
- Runs HMM inference on a single admixed individual
- Generates visualization

---

### 2. Simulated Analysis (Validation)
**Purpose**: Test HMM on synthetically created admixed individual  
**Runtime**: 5-15 minutes  
**Data**: Works with both slice and full datasets

```bash
cd src
python ../scripts/02_simulated_analysis.py
```

**Expected Output**:
- `outputs/simulated/simulated_admixed_karyogram.png` - Ancestry painting
- Console output showing detected vs true ancestry switches

**What it does**:
- Creates synthetic admixed individual (e.g., 50% YRI, 50% CEU)
- Runs HMM WITHOUT knowing the true ancestry
- Compares inferred ancestry to known solution
- Validates algorithm accuracy

---

### 3. Real Sample Analysis (Production)
**Purpose**: Analyze real admixed samples from 1000 Genomes  
**Samples**: ASW (African ancestry in Southwest USA)  
**Runtime**: 30-60 minutes (full dataset), 2-5 minutes (slice)  
**Data**: Recommended to use full dataset for robust results

```bash
cd src
python ../scripts/03_real_sample_analysis.py
```

**Expected Output**:
```
outputs/real_samples/
├── asw_ancestry_results.csv      # Per-sample statistics
├── asw_karyograms/
│   ├── sample_1.png
│   ├── sample_2.png
│   └── ...
└── asw_population_histogram.png  # Population-level summary
```

**What it does**:
- Loads all ASW samples from the dataset
- Infers local ancestry for each individual
- Aggregates statistics (switch rate, tract lengths, ancestry proportions)
- Generates individual plots and population summary

**Output Files**:
- `asw_ancestry_results.csv` columns:
  - `sample_id`: Sample identifier
  - `yri_proportion`: Estimated YRI ancestry %
  - `ceu_proportion`: Estimated CEU ancestry %
  - `switch_rate_per_mb`: Ancestry switches per megabase
  - `median_tract_length_kb`: Median ancestry tract size

---

### 4. Population Analysis (Advanced)
**Purpose**: Detailed population-level statistics  
**Runtime**: 1-2 hours (full dataset)  
**Data**: Requires full dataset

```bash
cd src
python ../scripts/04_population_analysis.py
```

**Expected Output**:
```
outputs/analysis/
├── population_stats.csv              # Aggregate statistics
├── ancestry_distribution.png         # Histogram of proportions
├── tract_length_distribution.png     # Genome-wide tract sizes
└── switch_rate_distribution.png      # Ancestry switch rates
```

---

## Benchmark Against RFMix

If you have RFMix results for comparison:

### Step 1: Export HMM Predictions

```bash
python benchmark/export_model_predictions.py \
  --vcf data/raw/vcf/query_asw_chr22.vcf.gz \
  --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
  --query-pop ASW \
  --out benchmark/predictions/model_predictions.csv
```

### Step 2: Compare with RFMix

```bash
python benchmark/compare_with_rfmix.py \
  --model benchmark/predictions/model_predictions.csv \
  --rfmix rfmix_predictions_chr22.csv \
  --sample-col sample_id \
  --position-col position \
  --model-label-col label \
  --rfmix-label-col label \
  --out benchmark/results/rfmix_comparison.csv
```

**Expected Output**: `rfmix_comparison.csv` with metrics:
- Concordance (% agreement)
- Cohen's kappa
- Switches per Mb
- Tract length statistics

---

## Data Options

All scripts automatically detect available data:

### Using Slice (Default)
- Data file: `data/processed/chr22_slice.vcf`
- Size: ~50MB
- Useful for: Testing, fast iteration
- Runtime: 5x faster than full

### Using Full Dataset
- Data file: `data/raw/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz`
- Size: ~1GB
- Useful for: Production analysis, robust results
- Runtime: Full (2-10 minutes per sample)

Scripts automatically use full dataset if available, fall back to slice.

---

## Running Tests

Validate data loading and components:

```bash
cd src
python -m pytest ../tests/
```

Or run individual test:
```bash
python ../tests/test_data_loading.py
```

---

## Output Files Location

| Analysis | Output Location |
|----------|-----------------|
| Sanity Check | `outputs/sanity_check/` |
| Simulated | `outputs/simulated/` |
| Real Samples | `outputs/real_samples/` |
| Population Stats | `outputs/analysis/` |
| RFMix Predictions | `evaluation/predictions/` |
| RFMix Comparison | `evaluation/results/` |

---

## Performance Tuning

### Speed Up Analysis
```bash
# Use process caching for frequencies (first run only)
python ../scripts/01_sanity_check.py --cache
```

### Reduce Memory Usage
```bash
# Use slice dataset
cd data
bash download_data_slice.sh
```

### Parallel Processing (if implemented)
```bash
# Some scripts support multiprocessing
python ../scripts/03_real_sample_analysis.py --jobs 4
```

---

## Troubleshooting

### "VCF file not found"
```bash
# Download data
cd data
bash download_data_slice.sh  # or download_data_full22.sh
```

### "Memory error"
- Reduce sample count in analysis script
- Use slice dataset instead
- Run analysis on machine with more RAM

### "Segmentation fault" (uncommon)
- Update pysam: `conda update -c bioconda pysam`
- Try slice dataset first
- Check VCF file integrity: `bcftools view -H <file.vcf.gz>`

### Script runs but produces no output
- Check `outputs/` folder structure exists
- Verify write permissions: `touch outputs/test.txt`
- Check console for error messages
