# Analysis Scripts

This folder contains the main analysis workflows for the HMM local ancestry inference project. Run them in sequence from the project root.

## Quick Start

```bash
cd cse284-hmm
conda activate hmm_env
python scripts/01_sanity_check.py
```

## Scripts

### 1. Sanity Check - `01_sanity_check.py`
**Purpose**: Quick validation test on a known sample  
**Runtime**: 5-10 minutes  
**Output**: `outputs/sanity_check/`

Tests the HMM pipeline on NA19625, a known admixed individual.

```bash
python scripts/01_sanity_check.py
```

### 2. Simulated Analysis - `02_simulated_analysis.py`
**Purpose**: Validate HMM performance on synthetic data  
**Runtime**: 5-15 minutes  
**Output**: `outputs/simulated/`

Creates a synthetic admixed individual and compares inferred vs true ancestry.

```bash
python scripts/02_simulated_analysis.py
```

### 3. Real Sample Analysis - `03_real_sample_analysis.py`
**Purpose**: Analyze real ASW (admixed) samples  
**Runtime**: 30-60 minutes (full data), 2-5 minutes (slice)  
**Output**: `outputs/real_samples/`

Processes all ASW samples from the dataset.

```bash
python scripts/03_real_sample_analysis.py
```

### 4. Population Analysis - `04_population_analysis.py`
**Purpose**: Population-level statistics and visualizations  
**Runtime**: 1-2 hours  
**Output**: `outputs/analysis/`

Aggregates statistics across all analyzed samples.

```bash
python scripts/04_population_analysis.py
```

## Running All Analyses

```bash
# Sequential execution (recommended for first-time)
python scripts/01_sanity_check.py
python scripts/02_simulated_analysis.py
python scripts/03_real_sample_analysis.py
python scripts/04_population_analysis.py
```

## Troubleshooting

### "Module not found" error
Ensure you're in the project root directory when running:
```bash
cd /path/to/cse284-hmm
python scripts/01_sanity_check.py  # ✓ Correct
```

Not from the scripts directory:
```bash
cd scripts
python 01_sanity_check.py  # ✗ Will fail with import errors
```

### "VCF file not found"
Download data first:
```bash
cd data
bash scripts/download_data_slice.sh  # or download_data_full22.sh
cd ..
python scripts/01_sanity_check.py
```

### Script crashes mid-execution
Check `outputs/` folder exists and is writable:
```bash
ls outputs/
touch outputs/test.txt  # Should work
```

## Output Structure

All scripts save results to `outputs/`:

```
outputs/
├── sanity_check/         # 01_sanity_check.py
│   └── *.png
├── simulated/            # 02_simulated_analysis.py
│   └── *.png
├── real_samples/         # 03_real_sample_analysis.py
│   ├── *.csv
│   └── *.png
└── analysis/             # 04_population_analysis.py
    ├── *.csv
    └── *.png
```

## For More Information

- **Setup**: See [docs/SETUP.md](../docs/SETUP.md)
- **Data**: See [docs/DATA_GUIDE.md](../docs/DATA_GUIDE.md)
- **Implementation**: See [docs/IMPLEMENTATION.md](../docs/IMPLEMENTATION.md)
- **Usage**: See [docs/USAGE.md](../docs/USAGE.md)
