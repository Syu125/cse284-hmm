# Analysis Scripts

This folder contains the main analysis workflows for the HMM local ancestry inference project. Run them in sequence from the project root.

## Quick Start

```bash
cd cse284-hmm
conda activate hmm_env
python scripts/01_real_sample_analysis.py
```

## Scripts

### 1. Real Sample Analysis - `01_real_sample_analysis.py`
**Purpose**: Analyze real ASW (admixed) samples  
**Runtime**: 30-60 minutes (full data), 2-5 minutes (slice)  
**Output**: `outputs/real_samples/`

Processes all ASW samples from the dataset and generates ancestry paintings.

```bash
python scripts/01_real_sample_analysis.py
```

### 2. Population Analysis - `02_population_analysis.py`
**Purpose**: Population-level statistics and visualizations  
**Runtime**: 1-2 hours  
**Output**: `outputs/analysis/`

Aggregates statistics across all analyzed samples for population-wide insights.

```bash
python scripts/02_population_analysis.py
```

## Running All Analyses

```bash
# Sequential execution (recommended for first-time)
python scripts/01_real_sample_analysis.py
python scripts/02_population_analysis.py
```

## Troubleshooting

### "Module not found" error
Ensure you're in the project root directory when running:
```bash
cd /path/to/cse284-hmm
python scripts/01_real_sample_analysis.py  # ✓ Correct
```

Not from the scripts directory:
```bash
cd scripts
python 01_real_sample_analysis.py  # ✗ Will fail with import errors
```

### "VCF file not found"
Download data first:
```bash
bash data/prepare_data.sh slice  # or prepare_data.sh full
python scripts/01_real_sample_analysis.py
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
├── real_samples/         # 01_real_sample_analysis.py
│   ├── *.csv
│   └── *.png
└── analysis/             # 02_population_analysis.py
    ├── *.csv
    └── *.png
```

## For More Information

- **Setup**: See [docs/SETUP.md](../docs/SETUP.md)
- **Data**: See [docs/DATA_GUIDE.md](../docs/DATA_GUIDE.md)
- **Implementation**: See [docs/IMPLEMENTATION.md](../docs/IMPLEMENTATION.md)
- **Usage**: See [docs/USAGE.md](../docs/USAGE.md)
