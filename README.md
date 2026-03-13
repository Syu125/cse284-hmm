# CSE284 HMM: Local Ancestry Inference

Implementation of a phased Hidden Markov Model for inferring local ancestry using the Viterbi algorithm on chromosome 22 from the 1000 Genomes Project.

## Quick Start

Before creating the environment, please make sure you have all the necessary system requirements needed to run the following commands. More information can be found [here](docs/SETUP.md).

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate hmm_env

# 2. Download a slice of the data (you can download the full version too, but for testing purposes stick to the slice)
bash data/prepare_data.sh slice

# 3. HMM analysis scripts
python general-analysis/01_real_sample_analysis.py
python general-analysis/02_population_analysis.py
```
## Benchmark
Benchmark now defaults to haplotype-level comparison (`YRI`, `CEU`) so all methods are evaluated in the same state space.
For complete benchmark setup and analysis, see [BENCHMARK.md](benchmark/BENCHMARK.md).

## Applying to a Dataset

This project is applied to a **public dataset** from the **1000 Genomes Project** (chromosome 22), including admixed **ASW** individuals.

### What is implemented

- **Phased haplotype HMM inference (`K`-state)**
	- Core model runs as a single-haplotype `K=2` ancestry decoder (`CEU`, `YRI`) with Viterbi decoding per haplotype.
	- Diploid combined labels are retained only as compatibility outputs for downstream tooling.
- **Per-sample real-data inference** via `general-analysis/01_real_sample_analysis.py`
	- Runs local ancestry inference on real ASW samples
	- Produces per-individual ancestry summaries and haplotype ancestry painting plots in `general-analysis/output/`
- **Population-level analysis** via `general-analysis/02_population_analysis.py`
	- Runs inference across ASW samples
	- Exports population summary statistics to `general-analysis/output/asw_ancestry_results.csv`
	- Produces cohort-level distribution plots in `general-analysis/output/`
- **Benchmark comparison against RFMix and FLARE** (same public data slice) via `benchmark/`
	- Produces concordance/kappa and tract-level metrics in `benchmark/results/`

### Interpretation summary
- On the chr22 slice benchmark, haplotype-level agreement (`YRI`, `CEU`) is high, with strong concordance across model pairs.
- Chance-adjusted agreement (Cohen's kappa) is moderate-to-strong depending on the pair, with best kappa in `flare_vs_rfmix` and strongest raw concordance in `hmm_vs_rfmix`.
- In repeated sample-size sweeps, metric variance decreases as sample size increases (especially by `N=50`), indicating more stable estimates with larger cohorts.
- Current evaluation prioritizes haplotype-consistent comparisons and stability trends over legacy diploid/heterozygous compatibility labels.

## Project Structure

```
cse284-hmm/
├── benchmark/                   # HMM vs RFMix/FLARE benchmarking
│   ├── predictions/
│   ├── results/
│   ├── BENCHMARK.md            # Benchmark runner usage and latest results
│   ├── run_benchmark.sh         # Single benchmark entrypoint
│   └── scripts/                 # Benchmark python scripts

├── data/                           # Data management
│   ├── cache/                     # Auto-generated frequency caches
│   ├── processed/                 # Sliced/preprocessed data
│   ├── raw/                       # Original 1000 Genomes files
│   └── prepare_data.sh            # Consolidated download & prep script
│
├── docs/                            # Documentation
│   └── SETUP.md                   # How to setup environment and download data
│
├── general-analysis/                # Analysis workflows on real ASW samples
│   ├── output/                    # Output folder for analysis scripts
│   ├── 01_real_sample_analysis.py # Per-sample haplotype ancestry inference
│   └── 02_population_analysis.py  # Population-level ancestry summaries
│
├── src/                             # HMM Implementation
│   ├── data/                      # VCF parsing and frequency calculation
│   ├── hmm/                       # Core algorithm (emission, transition, viterbi)
│   ├── visualization/             # Ancestry painting plots
│   ├── tests/                     # Unit tests
│   └── utils.py                   # Caching utilities for larger datasets

├── .gitignore                     # Includes the large data files
├── environment.yml                # Necessary libraries to install in environment
└── README.md                      # This
```

## LLM Usage
I used LLM assistance to help with the following:
- Reviewing my hmm implementation and helping to debug issues/make improvements
- Refactoring the model and analyses to a phased haplotype `K`-state workflow
- Revising documentation structure and tone
- Organizing my files and adding comments

## Citation

> The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.

## References


Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*.
