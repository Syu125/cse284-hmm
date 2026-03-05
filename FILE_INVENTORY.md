# File Inventory & Cleanup Guide

This document catalogs every file in the repository and explains its purpose. Use this to identify files for removal or retention.

---

## 🟢 KEEP - Core Project Files

### Root Level Configuration

| File | Purpose | Status |
|------|---------|--------|
| `README.md` | Main project documentation and overview | **KEEP** |
| `environment.yml` | Conda environment specification with dependencies | **KEEP** |
| `.gitignore` | Git ignore patterns for generated files | **KEEP** |
| `ORGANIZATION_GUIDE.md` | Repository organization guidelines (this document) | **KEEP** |

---

## 📚 Documentation (`docs/`)

| File | Purpose | Status |
|------|---------|--------|
| `docs/SETUP.md` | Environment setup and installation instructions | **KEEP** |
| `docs/USAGE.md` | Step-by-step guide for running analyses | **KEEP** |
| `docs/DATA_GUIDE.md` | Detailed description of data files and sources | **KEEP** |
| `docs/IMPLEMENTATION.md` | Technical HMM implementation details with equations | **KEEP** |

**Recommendation**: All essential for user onboarding.

---

## 🎯 Analysis Scripts (`scripts/`)

| File | Purpose | Status |
|------|---------|--------|
| `scripts/README.md` | Guide for running analysis scripts in sequence | **KEEP** |
| `scripts/01_real_sample_analysis.py` | Analyze real ASW (admixed) samples | **KEEP** |
| `scripts/02_population_analysis.py` | Population-level statistics and aggregation | **KEEP** |

**Recommendation**: Core workflow files - all essential.

---

## 🧬 Data Management (`data/`)

### Documentation

| File | Purpose | Status |
|------|---------|--------|
| `data/README.md` | Data folder structure and file descriptions | **KEEP** |

### Download & Preprocessing Scripts (`data/scripts/`)

| File | Purpose | Status |
|------|---------|--------|
| `download_data_slice.sh` | Download small chromosome 22 slice (~50MB, **primary**) | **KEEP** |
| `download_data_full22.sh` | Download full chromosome 22 (~1GB) | **KEEP** |
| `prepare_rfmix_vcfs.sh` | Prepare VCF files for RFMix comparison | **KEEP** (if benchmarking) |
| `run_rfmix_slice.sh` | Run RFMix on slice data for comparison | **KEEP** (if benchmarking) |
| `clean_genetic_map.py` | Clean/reformat genetic map file | **KEEP** (utility) |
| `generate_rfmix_map.py` | Generate RFMix-compatible genetic map | **KEEP** (if benchmarking) |

**Recommendation**: Keep all unless you don't plan to use RFMix for benchmarking. If no RFMix:
- Can remove: `prepare_rfmix_vcfs.sh`, `run_rfmix_slice.sh`, `generate_rfmix_map.py`

---

## 💻 Source Code (`src/`)

### Core HMM Implementation (`src/hmm/`)

| File | Purpose | Status |
|------|---------|--------|
| `hmm/emission.py` | Emission probability model (P(genotype \| ancestry)) | **KEEP** |
| `hmm/transition.py` | Transition probability model (recombination-based) | **KEEP** |
| `hmm/viterbi.py` | Viterbi algorithm for ancestry inference | **KEEP** |

**Recommendation**: Core algorithm - all essential.

### Data Handling (`src/data/`)

| File | Purpose | Status |
|------|---------|--------|
| `data/__init__.py` | Package marker (empty or minimal) | **KEEP** |
| `data/data_parser.py` | VCF parsing, frequency calculation, genetic map handling | **KEEP** |

**Recommendation**: Essential data I/O.

### Visualization (`src/visualization/`)

| File | Purpose | Status |
|------|---------|--------|
| `visualization/__init__.py` | Package initialization, exports `plot_ancestry` | **KEEP** |
| `visualization/karyogram.py` | Ancestry painting plot generation | **KEEP** |

**Recommendation**: Essential for results visualization.

### Testing (`src/tests/`)

| File | Purpose | Status |
|------|---------|--------|
| `tests/__init__.py` | Package marker for test module | **KEEP** |
| `tests/test_data_loading.py` | Unit tests for data loading and emission model | **KEEP** |

**Recommendation**: Good practice to keep tests.

---

## ⚠️ REVIEW - Potentially Redundant Files

### `src/utils.py`
- **Purpose**: Caching utility for allele frequencies
- **Content**: `get_cached_frequencies()` function
- **Status**: **KEEP** (useful utility, used by scripts)
- **Note**: Contains proper caching logic with directory creation

### `src/visualization.py`
- **Purpose**: DEPRECATED - Old visualization file
- **Content**: Only a deprecation comment, no actual code
- **Status**: **🗑️ REMOVE** (functionality moved to `visualization/` package)
- **Reason**: No longer needed since visualization package exists

### `src/precalculate_frequencies.py`
- **Purpose**: Standalone script to pre-cache frequencies
- **Content**: Example usage of caching (not a module)
- **Status**: **⚠️ OPTIONAL** - Can remove or keep as example
- **Reason**: Scripts already handle caching via `utils.py`
- **Recommendation**: **Remove** since caching is automatic in analysis scripts

### `src/checker.py`
- **Purpose**: Debug/development script to check data loading
- **Content**: Tests panel parsing, frequency calculation, emission model
- **Status**: **🗑️ REMOVE** (covered by `tests/test_data_loading.py`)
- **Reason**: Redundant with proper test file
- **Differences from test**: Less comprehensive, older code

---

## 📊 Evaluation & Benchmarking (`benchmark/`)

| File | Purpose | Status |
|------|---------|--------|
| `benchmark/README.md` | Guide for running RFMix comparison | **KEEP** |
| `benchmark/__init__.py` | Package marker | **KEEP** |
| `benchmark/export_model_predictions.py` | Export HMM predictions to CSV for comparison | **KEEP** |
| `benchmark/compare_with_rfmix.py` | Compare HMM vs RFMix at SNP level | **KEEP** |
| `benchmark/convert_rfmix_to_snp_csv.py` | Convert RFMix output format to CSV | **KEEP** (if using RFMix) |
| `benchmark/metrics.py` | Evaluation metrics (concordance, kappa, switches) | **KEEP** |

**Recommendation**: Keep all if you plan to benchmark against RFMix. Otherwise:
- Keep: `export_model_predictions.py`, `metrics.py` (useful for any evaluation)
- Optional: `compare_with_rfmix.py`, `convert_rfmix_to_snp_csv.py` (RFMix-specific)

---

## 🗂️ Empty Folders (with `.gitkeep`)

These folders exist to receive generated files:

```
data/raw/vcf/
data/raw/panels/
data/raw/maps/
data/processed/
data/cache/
outputs/real_samples/
outputs/analysis/
benchmark/predictions/
benchmark/results/
```

**Status**: **KEEP** - Essential directory structure

---

## 📋 Summary & Recommendations

### 🗑️ Files to Remove (3 files)

1. **`src/visualization.py`** - Deprecated, empty stub
2. **`src/checker.py`** - Debug script, redundant with `tests/test_data_loading.py`
3. **`src/precalculate_frequencies.py`** - Example script, redundant (caching is automatic)

### ⚠️ Optional to Remove (if not using RFMix benchmarking)

4. `data/scripts/prepare_rfmix_vcfs.sh`
5. `data/scripts/run_rfmix_slice.sh`
6. `data/scripts/generate_rfmix_map.py`
7. `benchmark/convert_rfmix_to_snp_csv.py`
8. `benchmark/compare_with_rfmix.py`

### ✅ Essential Files to Keep (26 core files)

- **Root**: 4 files (README, environment.yml, .gitignore, ORGANIZATION_GUIDE.md)
- **Documentation**: 5 files (4 guides + data README)
- **Scripts**: 5 files (4 analysis scripts + README)
- **Source code**: 10 files (HMM, data, visualization, tests, utils)
- **Evaluation**: 2 files (export predictions, metrics)

---

## 🚀 Cleanup Commands

### Remove Confirmed Redundant Files
```bash
# From repository root
rm src/visualization.py
rm src/checker.py
rm src/precalculate_frequencies.py
```

### Optional: Remove RFMix-specific Files (if not benchmarking)
```bash
rm data/scripts/prepare_rfmix_vcfs.sh
rm data/scripts/run_rfmix_slice.sh
rm data/scripts/generate_rfmix_map.py
rm benchmark/convert_rfmix_to_snp_csv.py
rm benchmark/compare_with_rfmix.py
```

---

## 🎯 Final Clean Repository Stats

After removing the 3 redundant files:

- **Total essential files**: ~26-31 (depending on RFMix decision)
- **Lines of code**: ~2500 (core implementation)
- **Documentation**: 5 comprehensive guides
- **Ready for**: Publication, sharing, collaboration

---

## 📝 Notes

- All files have been reviewed against current usage in analysis scripts
- Redundant files identified based on:
  - Deprecated code
  - Duplicate functionality
  - Covered by better implementations
- **No data files** are present (clean state for first-time users)
- All folder structures preserved for generated outputs
