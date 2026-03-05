# CSE284 HMM: Local Ancestry Inference

Implementation of a Hidden Markov Model for inferring local ancestry using the Viterbi algorithm on chromosome 22 from the 1000 Genomes Project.

## Quick Start

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate hmm_env

# 2. Download data (~50MB slice recommended for testing)
bash data/prepare_data.sh slice

# 3. Run analysis
python scripts/01_real_sample_analysis.py
```

## Project Structure

```
cse284-hmm/
├── data/                        # Data management
│   ├── prepare_data.sh         # Consolidated download & prep script
│   ├── raw/                    # Original 1000 Genomes files
│   ├── processed/              # Sliced/preprocessed data
│   └── cache/                  # Auto-generated frequency caches
│
├── src/                         # HMM Implementation
│   ├── hmm/                    # Core algorithm (emission, transition, viterbi)
│   ├── data/                   # VCF parsing & frequency calculation
│   ├── visualization/          # Ancestry painting plots
│   ├── tests/                  # Unit tests
│   └── utils.py                # Caching utilities
│
├── scripts/                     # Analysis Workflows
│   ├── 01_real_sample_analysis.py
│   └── 02_population_analysis.py
│
├── benchmark/                   # RFMix Benchmarking
│   ├── export_model_predictions.py
│   ├── compare_with_rfmix.py
│   ├── metrics.py
│   ├── predictions/
│   └── results/
│
└── docs/                        # Documentation
    ├── SETUP.md
    ├── USAGE.md
    ├── DATA_GUIDE.md
    └── IMPLEMENTATION.md
```

## Three Main Components

### 1. **Data (`data/`)**
- Single consolidated script: `prepare_data.sh`
- Handles downloading from 1000 Genomes Project
- Organizes files into `raw/`, `processed/`, and `cache/` folders

### 2. **Source (`src/`)**
- Core HMM implementation
- Data I/O utilities
- Visualization and testing

### 3. **Benchmark (`benchmark/`)**
- RFMix comparison tools
- Model prediction export
- Evaluation metrics

## Getting Started

See [docs/SETUP.md](docs/SETUP.md) for detailed installation and setup instructions.

## Running Analyses

1. **Real Samples** (30 min) - Analyze ASW (admixed) individuals
   ```bash
   python scripts/01_real_sample_analysis.py
   ```

2. **Population** (1+ hr) - Aggregate population statistics
   ```bash
   python scripts/02_population_analysis.py
   ```

See [docs/USAGE.md](docs/USAGE.md) for complete step-by-step guide.

## Benchmarking Against RFMix

Export predictions and compare:
```bash
python benchmark/export_model_predictions.py \
  --vcf data/raw/vcf/query_asw_chr22.vcf.gz \
  --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
  --query-pop ASW \
  --out benchmark/predictions/model_predictions.csv

python benchmark/compare_with_rfmix.py \
  --model benchmark/predictions/model_predictions.csv \
  --rfmix rfmix_predictions.csv \
  --out benchmark/results/comparison.csv
```

## Documentation

- [Setup & Installation](docs/SETUP.md)
- [Running Analyses](docs/USAGE.md)
- [Data Guide](docs/DATA_GUIDE.md)
- [Technical Implementation](docs/IMPLEMENTATION.md)
- [Organization Guide](ORGANIZATION_GUIDE.md)
- [File Inventory & Cleanup](FILE_INVENTORY.md)

## Key Features

- **Hardy-Weinberg Emission Model**: P(genotype | ancestry state)
- **Recombination-based Transitions**: Uses genetic distance for realistic transitions
- **Viterbi Algorithm**: Efficient DP in log-space for ancestry inference
- **Caching**: Auto-caches frequency calculations for speed
- **RFMix Comparison**: Benchmark against industry-standard tool

## Citation

If using this for research, cite the 1000 Genomes Project:

> The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.

## References

Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*.