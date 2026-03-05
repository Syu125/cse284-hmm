# CSE284 HMM: Local Ancestry Inference

Implementation of a Hidden Markov Model for inferring local ancestry using the Viterbi algorithm on chromosome 22 from the 1000 Genomes Project.

## Quick Start

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate hmm_env

# 2. Download data
bash data/prepare_data.sh slice

# 3. Run HMM analysis
python scripts/01_real_sample_analysis.py
python scripts/02_population_analysis.py
```
## Benchmark
More information can be found [here](benchmark/README.md)

## Project Structure

```
cse284-hmm/
├── benchmark/                   # RFMix Benchmarking
│   ├── export_model_predictions.py
│   ├── compare_with_rfmix.py
│   ├── metrics.py
│   ├── predictions/
│   └── results/
|
├── data/                        # Data management
│   ├── prepare_data.sh         # Consolidated download & prep script
│   ├── raw/                    # Original 1000 Genomes files
│   ├── processed/              # Sliced/preprocessed data
│   └── cache/                  # Auto-generated frequency caches
│
├── scripts/                     # Analysis Workflows
│   ├── 01_real_sample_analysis.py
│   └── 02_population_analysis.py
│
├── src/                         # HMM Implementation
│   ├── data/                   # VCF parsing & frequency calculation
│   ├── hmm/                    # Core algorithm (emission, transition, viterbi)
│   ├── visualization/          # Ancestry painting plots
│   ├── tests/                  # Unit tests
│   └── utils.py                # Caching utilities
│
└── docs/                        # Documentation
    ├── SETUP.md                # How to setup environment and download data
    └── IMPLEMENTATION.md
```

## Citation

> The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.

## References

Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*.