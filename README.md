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
*Note: I'm still working on this. The benchmark script is complete, but I still need to analyze and evaluate the results.*

More information can be found [here](benchmark/README.md).

## Applying to a Dataset
*Note: This is also on my TODO - not implemented yet.*

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