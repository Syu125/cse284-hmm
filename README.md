# CSE284 HMM: Local Ancestry Inference

Implementation of a Hidden Markov Model for inferring local ancestry using the Viterbi algorithm on chromosome 22 from the 1000 Genomes Project.

## Quick Start

Before creating the environment, please make sure you have all the necessary system requirements. They can be found [here](docs/SETUP.md).

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate hmm_env

# 2. Download a slice of the data (you can download the full version too, but for testing purposes stick to the slice)
bash data/prepare_data.sh slice

# 3. HMM analysis scripts
python scripts/01_real_sample_analysis.py
python scripts/02_population_analysis.py
```
## Benchmark
Benchmark uses a 3-class comparison (`YRI`, `CEU`, `HET`) by default.
For complete benchmark setup, commands, options, and troubleshooting, see [benchmark/README.md](benchmark/README.md).
Use `--valid-labels YRI,CEU` only as a homozygous-only diagnostic mode.

## Applying to a Dataset

This project is applied to a **public dataset** from the **1000 Genomes Project** (chromosome 22), including admixed **ASW** individuals.

### What is implemented
- **Per-sample real-data inference** via `scripts/01_real_sample_analysis.py`
	- Runs local ancestry inference on real ASW samples
	- Produces per-individual ancestry summaries and ancestry painting plots in `scripts/output/`
- **Population-level analysis** via `scripts/02_population_analysis.py`
	- Runs inference across ASW samples
	- Exports population summary statistics to `scripts/output/asw_ancestry_results.csv`
	- Produces a distribution plot `scripts/output/asw_population_histogram.png`
- **Benchmark comparison against RFMix** (same public data slice) via `benchmark/`
	- Produces concordance/kappa and tract-level metrics in `benchmark/results/`

### Interpretation summary
- On the chr22 slice benchmark, 3-class agreement (`YRI`, `CEU`, `HET`) is strong overall, with mean concordance around **0.78–0.84** depending on sample set.
- Chance-adjusted agreement (Cohen's kappa) is in the **moderate-to-substantial** range (~**0.48–0.68**), indicating non-trivial agreement beyond class prevalence effects.
- In repeated sample-size sweeps, metric variance decreases as sample size increases (notably at `N=50`), suggesting more stable evaluation at larger cohort sizes.
- For this reason, kappa and stability trends are treated as primary evidence, while binary-only labels are used only as diagnostics.

See `benchmark/BENCHMARK.md` for full commands, tables, and detailed interpretation.

## Project Structure

```
cse284-hmm/
├── benchmark/                   # RFMix Benchmarking
│   |── data/
│   ├── predictions/
│   |── results/
│   ├── compare_with_rfmix.py
│   ├── convert_rfmix_to_snp_csv.py
│   ├── export_model_predictions.py
│   ├── fix_genetic_map.py
│   ├── metrics.py
│   └── prepare_benchmark_data.py
|
├── data/                           # Data management
│   |── cache/                     # Auto-generated frequency caches
│   |── processed/                 # Sliced/preprocessed data
│   ├── raw/                       # Original 1000 Genomes files
│   └── prepare_data.sh            # Consolidated download & prep script
│
├── docs/                            # Documentation
|   ├── SETUP.md                   # How to setup environment and download data
|   └── IMPLEMENTATION.md          # Breakdown of the HMM
│
├── scripts/                         # Analysis Workflows
│   ├── output/                    # Output folder for the two scripts
│   ├── 01_real_sample_analysis.py # Analyze some samples
│   └── 02_population_analysis.py  # Analyze full population
│
├── src/                             # HMM Implementation
|   ├── data/                      # VCF parsing & frequency calculation
|   ├── hmm/                       # Core algorithm (emission, transition, viterbi)
|   ├── visualization/             # Ancestry painting plots
|   ├── tests/                     # Unit tests
|   └── utils.py                   # Caching utilities when running full datasets
|
├── .gitignore                     # Includes the large data files
├── environment.yml                # Necessary libraries to install in environment
└── README.md                      # This
```

## Citation

> The 1000 Genomes Project Consortium. (2015). A global reference for human genetic variation. *Nature*, 526(7571), 68-74.

## References

Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*.