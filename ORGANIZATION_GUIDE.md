# Repository Organization Guide

## Current State Analysis

Your CSE284 HMM project implements a **Hidden Markov Model for Local Ancestry Inference** using genetic data from the 1000 Genomes Project. The project is functional but has organizational opportunities for improved clarity and reproducibility.

---

## Recommended Folder Structure

```
cse284-hmm/
├── README.md                           # Main project overview
├── ORGANIZATION_GUIDE.md               # This file
├── environment.yml                     # Conda environment (KEEP HERE)
├── .gitignore                          # Version control exclusions
│
├── docs/                               # NEW: Documentation
│   ├── SETUP.md                        # Environment setup instructions
│   ├── USAGE.md                        # How to run analyses
│   ├── DATA_GUIDE.md                   # Data sources and file descriptions
│   └── IMPLEMENTATION.md               # Technical implementation details
│
├── data/                               # Raw data and data scripts (REFACTOR)
│   ├── README.md                       # Data manifest
│   ├── .gitkeep                        # Ensure directory tracks
│   ├── raw/                            # NEW: Raw downloaded files
│   │   ├── vcf/                        # VCF files
│   │   ├── maps/                       # Genetic maps
│   │   └── panels/                     # Population panels
│   ├── processed/                      # NEW: Preprocessed data
│   │   ├── chr22_slice.vcf*
│   │   └── sliced_data/
│   ├── scripts/                        # NEW: Data prep scripts
│   │   ├── download_data_full22.sh
│   │   ├── download_data_slice.sh
│   │   ├── prepare_rfmix_vcfs.sh
│   │   ├── run_rfmix_slice.sh
│   │   ├── clean_genetic_map.py
│   │   ├── generate_rfmix_map.py
│   │   └── precalculate_frequencies.py
│   └── cache/                          # NEW: Cached frequency files (.pkl)
│
├── src/                                # Source code
│   ├── __init__.py
│   ├── hmm/                            # Core HMM implementation
│   │   ├── __init__.py
│   │   ├── emission.py
│   │   ├── transition.py
│   │   └── viterbi.py
│   ├── data/                           # Data I/O utilities
│   │   ├── __init__.py
│   │   └── data_parser.py
│   ├── visualization/                  # Plotting utilities
│   │   ├── __init__.py
│   │   └── karyogram.py
│   ├── utils.py                        # Shared utilities
│   └── tests/                          # Unit tests
│       ├── __init__.py
│       └── test_data_loading.py
│
├── scripts/                            # NEW: Analysis workflows
│   ├── 01_simulated_analysis.py        # Test on synthetic admixed
│   ├── 02_real_sample_analysis.py      # Analyze ASW samples
│   ├── 03_population_analysis.py       # Population-level stats
│   └── README.md                       # Running order and descriptions
│
├── benchmark/                          # NEW: Benchmarking against RFMix
│   ├── README.md
│   ├── export_predictions.py           # Export HMM predictions
│   ├── compare_with_rfmix.py
│   ├── convert_rfmix_to_snp_csv.py
│   ├── metrics.py
│   ├── predictions/
│   │   ├── model_predictions_chr22.csv
│   │   └── model_slice_predictions.csv
│   └── results/
│       └── rfmix_comparison_summary.csv
│
├── outputs/                            # NEW: Results and figures
│   ├── simulated/
│   ├── real_samples/
│   │   ├── asw_ancestry_results.csv
│   │   ├── karyogram_*.png
│   │   └── asw_population_histogram.png
│   └── analysis/
│       └── population_level_stats.csv
│
└── requirements.txt                    # NEW: pip requirements (optional backup)
```

---

## File Reorganization Details

### 1. **Documentation (NEW: `docs/` folder)**
Move documentation from scattered locations into dedicated files:

- **SETUP.md**: Conda environment creation, WSL setup for Windows users
- **USAGE.md**: How to run each analysis script in sequence
- **DATA_GUIDE.md**: Description of all data files, their sources, and expected formats
- **IMPLEMENTATION.md**: Technical details on HMM model, equations, and algorithm

**Benefits**: Clear separation of concerns; easy for new users to get started

---

### 2. **Data Reorganization (within `data/` folder)**

**Current issue**: Raw data, scripts, and outputs are mixed together

**New structure**:
- **`raw/`**: Files downloaded from 1000 Genomes (gitignored)
- **`processed/`**: Cleaned/sliced versions ready for analysis
- **`scripts/`**: All data preparation bash/Python scripts
- **`cache/`**: Cached pickled frequency files for faster loading

**Update `.gitignore`**:
```
# Data files (too large for Git)
data/raw/**/*.vcf
data/raw/**/*.vcf.gz*
data/processed/**/*.vcf
data/processed/**/*.vcf.gz*
data/cache/*.pkl

# Python cache
src/__pycache__/
src/**/__pycache__/

# Outputs
outputs/
```

---

### 3. **Script Organization (NEW: `scripts/` folder)**

Instead of loose Python files in `src/`:
- **`01_simulated_analysis.py`** ← rename from `simulated_admixed.py`
- **`02_real_sample_analysis.py`** ← rename from `analyze_real.py`
- **`03_population_analysis.py`** ← rename from `population_analysis.py`

Add **`scripts/README.md`** with execution instructions:
```
# Running Analyses

1. **Simulated Analysis** (2-5 min):
   python 01_simulated_analysis.py
   Output: outputs/simulated/

2. **Real Sample Analysis** (30-60 min):
   python 02_real_sample_analysis.py
   Output: outputs/real_samples/

3. **Population Analysis** (varies):
   python 03_population_analysis.py
   Output: outputs/analysis/
```

---

### 4. **Evaluation Reorganization (rename `eval/` → `benchmark/`)**

- Consolidate benchmarking code
- Create `predictions/` and `results/` subdirectories for clarity
- Add README explaining comparison workflow

---

### 5. **Outputs Organization (NEW: `outputs/` folder)**

**Current issue**: Output files (`.png`, `.csv`) scattered in `src/`

**New structure**:
```
outputs/
├── simulated/
│   └── (simulated_admixed outputs)
├── real_samples/
│   ├── asw_ancestry_results.csv
│   ├── karyogram_*.png
│   └── population_histograms.png
└── analysis/
    └── population_level_stats.csv
```

**Update code** to write here instead of `src/`

---

## File Cleanup Recommendations

### Remove or Archive:
- **`checker.py`**: Unclear purpose—document or delete
- **`population_analysis.py`** in `src/`: Move to `scripts/03_population_analysis.py`
- **`Miniconda3-latest-Linux-x86_64.sh`**: Delete or move to `docs/` if needed for reference

### Rename:
- `evaluation/` → `benchmark/` (clearer naming)
- `simulated_admixed.py` → `01_simulated_analysis.py`
- `analyze_real.py` → `02_real_sample_analysis.py`

---

## Quick Start for Reproducibility

After reorganization, a new user should:

1. Clone repo
2. Follow `docs/SETUP.md` to create conda environment
3. Run `bash data/prepare_data.sh slice` (or `full` for complete dataset)
4. Follow `scripts/README.md` to run analyses
5. Check `outputs/` for results

---

## Implementation Priority

### Phase 1 (Essential - 2-3 hours):
1. Create `docs/` folder with SETUP.md, USAGE.md
2. Create `scripts/` folder and move analysis files
3. Create `outputs/` folder structure
4. Update `.gitignore` and add `.gitkeep` files

### Phase 2 (Important - 2-3 hours):
1. Reorganize `data/` folder structure
2. Update all Python imports to reflect new paths
3. Create `DATA_GUIDE.md`
4. Create `IMPLEMENTATION.md`

### Phase 3 (Polish - 1-2 hours):
1. Rename `eval/` → `benchmark/`
2. Add docstrings to scripts
3. Create `scripts/README.md` with examples
4. Final validation: run scripts from new locations

---

## Summary of Benefits

| Aspect | Before | After |
|--------|--------|-------|
| **Clarity** | Mixed files | Clear purpose for each folder |
| **Reproducibility** | Scattered docs | Centralized docs with setup/usage |
| **Scalability** | Hard to add new analyses | Easy: add to `scripts/` folder |
| **Collaboration** | Unclear how to run | Step-by-step documented workflow |
| **Output tracking** | Files everywhere | Organized in `outputs/` |
| **Data management** | Raw and processed mixed | Separated with cache support |

---

## Version Control Tips

```bash
# After reorganization
git add .
git commit -m "refactor: reorganize repository structure for clarity and reproducibility"

# If moving large files
git mv src/old_script.py scripts/01_new_script.py
git mv data/output.csv outputs/results/output.csv
```

---

## Notes

- **Don't delete files immediately** — comment them out first to ensure nothing breaks
- **Update relative paths** in all Python scripts after moving files
- **Test each script** from its new location before considering reorganization complete
- **Add `.gitkeep`** files to empty directories so Git tracks folder structure
- Consider adding a **`Makefile`** later for common tasks (setup, run tests, clear cache)
