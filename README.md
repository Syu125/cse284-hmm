# cse284-hmm

Implementation of HMM for Local Ancestry Inference using the Viterbi algorithm.

## Project Structure

```
src/
├── hmm/                          # Core HMM implementation
│   ├── __init__.py
│   ├── emission.py              # Emission probability model
│   ├── transition.py            # Transition probability model  
│   └── viterbi.py               # Viterbi algorithm for inference
│
├── data/                         # Data loading and parsing
│   ├── __init__.py
│   └── data_parser.py           # Functions for parsing VCF, panels, genetic maps
│
├── visualization/               # Visualization utilities
│   ├── __init__.py
│   └── karyogram.py             # Ancestry painting plots
│
├── tests/                        # Test and validation scripts
│   ├── __init__.py
│   └── test_data_loading.py     # Validate data loading and components
│
├── utils.py                      # Utility functions (caching, etc.)
│
├── sanity_check_na19625.py       # Quick test on a single sample  
├── simulated_admixed.py          # Test on simulated admixed individual
├── analyze_real.py               # Analyze real admixed samples
└── population_analysis.py        # Population-level analysis
```

## Key Components

### HMM Model (`src/hmm/`)
- **Emission Model**: Calculates `P(genotype | ancestry state)` using Hardy-Weinberg equilibrium
- **Transition Model**: Calculates `P(state[i+1] | state[i])` based on genetic distance  
- **Viterbi Algorithm**: Finds the most likely sequence of ancestry states

### Data Handling (`src/data/`)
- Parses 1000 Genomes VCF files
- Extracts allele frequencies by population
- Interpolates genetic map positions

## Setup

### Windows Users
If using Windows, install WSL first:
```bash
wsl --install
```

### Dependencies
Create conda environment:
```bash
cd cse284-hmm
conda env create -f environment.yml
conda activate hmm_env
```

## Getting Data

The pipeline requires three data files from 1000 Genomes Phase 3:

1. **Genetic Map** - Physical to genetic position mapping
2. **Sample Panel** - Population assignments for samples (YRI, CEU, ASW, etc.)  
3. **VCF File** - Genotypes for chromosome 22

To download data:
```bash
cd data/
bash download_data_slice.sh    # Downloads limited slice (faster for testing)
# OR
bash download_data_full22.sh   # Downloads full chromosome 22
```

## Running Analyses

### 1. Sanity Check (Quick Start)
```bash
cd src
python sanity_check_na19625.py
```
Tests the pipeline on a single known individual (NA19625).

### 2. Simulated Analysis
```bash
python simulated_admixed.py
```
Creates a synthetic admixed individual and validates HMM detection of ancestry switches.

### 3. Real Sample Analysis
```bash
python analyze_real.py
```
Analyzes a few real admixed (ASW) samples from the dataset.

### 4. Population Analysis  
```bash
python population_analysis.py
```
Runs inference on all ASW samples and generates population-level statistics.

### 5. Data Validation
```bash
python -m pytest tests/test_data_loading.py
# OR
python tests/test_data_loading.py
```
Validates that data loading, frequency calculation, and emission model work correctly.

### 6. Benchmark Against RFMix
Export your HMM per-SNP predictions:
```bash
python eval/export_model_predictions.py --query-pop ASW --limit-samples 5 --out model_predictions_chr22.csv
```

Compare to an RFMix SNP-level CSV with matching columns:
```bash
python eval/compare_with_rfmix.py \
	--model model_predictions_chr22.csv \
	--rfmix rfmix_predictions_chr22.csv \
	--sample-col sample_id \
	--position-col position \
	--model-label-col label \
	--rfmix-label-col label \
	--out rfmix_comparison_summary.csv
```

Expected columns in both files:
- `sample_id`
- `position`
- `label` (supported labels: `YRI/CEU`, `AFR/EUR`, `0/1`)

Summary outputs include concordance, Cohen's kappa, switches per Mb, median tract length, and global YRI percentage differences.

## Performance Notes

- **Frequency Caching**: Use `utils.get_cached_frequencies()` to cache expensive VCF frequency calculations
- **Memory**: Processing full chromosome 22 with all samples requires significant memory - start with the slice
- **Runtime**: Full dataset analysis takes ~5-10 minutes depending on sample count

## Implementation Details

### Hardy-Weinberg Emission Model
For allele frequency $p$ and genotype $(a_1, a_2)$:
$$P(\text{genotype} | \text{state}) = P(a_1) \cdot P(a_2) = p^{a_1}(1-p)^{1-a_1} \cdot p^{a_2}(1-p)^{1-a_2}$$

### Transition Model  
States transition based on recombination rate:
$$P(\text{switch}) = 1 - e^{-G \cdot d}$$
where $G$ = generations since admixture and $d$ = genetic distance in Morgans

### Viterbi
Dynamic programming in log-space to find:
$$\arg\max_{\text{path}} \prod_i P(g_i | s_i) \cdot P(s_i | s_{i-1})$$

## References

- Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm
- The 1000 Genomes Project Consortium (2015). A global reference for human genetic variation