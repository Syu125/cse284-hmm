# Model Benchmarking & RFMix Comparison

This folder contains tools and results for evaluating and comparing the HMM model against RFMix.

## Folder Structure

```
benchmark/
├── export_model_predictions.py      # Export HMM predictions to CSV
├── compare_with_rfmix.py            # Compare against RFMix results
├── convert_rfmix_to_snp_csv.py      # Convert RFMix output format
├── metrics.py                       # Evaluation metrics  
├── __init__.py
├── predictions/                     # HMM model predictions
│   └── .gitkeep
└── results/                         # Evaluation results
    └── .gitkeep
```

## Quick Start

### 1. Export HMM Predictions

```bash
python benchmark/export_model_predictions.py \
  --vcf data/raw/vcf/query_asw_chr22.vcf.gz \
  --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
  --query-pop ASW \
  --out benchmark/predictions/model_predictions.csv
```

### 2. Compare with RFMix

```bash
python benchmark/compare_with_rfmix.py \
  --model benchmark/predictions/model_predictions.csv \
  --rfmix rfmix_predictions_chr22.csv \
  --sample-col sample_id \
  --position-col position \
  --model-label-col label \
  --rfmix-label-col label \
  --out benchmark/results/comparison.csv
```

## Scripts

### `export_model_predictions.py`
Export per-SNP ancestry predictions from your HMM model for benchmarking.

**Usage**:
```bash
python benchmark/export_model_predictions.py --help
```

**Required arguments**:
- `--vcf`: Path to query VCF file
- `--panel`: Path to sample panel file
- `--map`: Path to genetic map file
- `--query-pop`: Population code (e.g., ASW)
- `--out`: Output CSV file path

**Output format**:
```
sample_id,position,label
NA20509,50001234,YRI
NA20509,50002456,YRI  
NA20509,50003789,CEU
...
```

### `compare_with_rfmix.py`
Compare HMM predictions with RFMix results at SNP level.

**Usage**:
```bash
python benchmark/compare_with_rfmix.py --help
```

**Output metrics**:
- Concordance: % agreement between methods
- Cohen's kappa: Adjusted agreement (0-1, higher is better)
- Switches per Mb: Ancestry switch rate
- Tract lengths: Median and mean ancestry tract sizes
- Ancestry composition: % YRI vs %CEU difference

### `metrics.py`
Utility functions for computing evaluation metrics.

**Functions**:
- `compute_concordance()` - SNP-level agreement
- `compute_kappa()` - Cohen's kappa statistic
- `compute_switch_rate()` - Ancestry switch statistics
- `compute_tract_lengths()` - Ancestry tract analysis

### `convert_rfmix_to_snp_csv.py`
Convert RFMix output files to SNP-level CSV format.

**RFMix outputs**:
- `.msp.tsv` - Ancestry state per segment
- `.sis.tsv` - Posterior probabilities
- `.rfmix.Q` - Global ancestry proportions

## Results

### `rfmix_slice_snp_level.csv`
SNP-level comparison on chromosome 22 slice:
- Sample count: 20 ASW samples
- SNP count: ~20k variants
- Concordance: ~93-95%

### `rfmix_comparison_summary.csv`
Summary statistics from RFMix comparison:
- Per-sample metrics
- Population-level aggregates
- Ancestry composition estimates

## Expected Results

### Typical Concordance
| Dataset | Concordance | Cohen's Kappa |
|---------|------------|---------------|
| Slice (20k SNPs) | 90-95% | 0.85-0.92 |
| Full (490k SNPs) | 92-96% | 0.88-0.94 |

**Note**: Exact values depend on:
- HMM parameter tuning (generations, prior)
- RFMix parameter settings
- Population-specific differences

## Troubleshooting

### Export script fails: "Sample not found"
- Check query population code matches panel file
- Verify VCF file contains samples from that population
- Review the panel file to list available population codes

### Comparison script gives low concordance
- Check both files use same position/label columns
- Verify identical SNP set is being compared
- Compare specific samples with `--verbose` flag

### Missing predictions
- Ensure HMM model ran successfully (check `outputs/` folder)
- Verify data files exist and have correct format
- Check for memory/timeout errors in stdout

## For More Information

- [USAGE.md](../docs/USAGE.md) for running analyses
- [IMPLEMENTATION.md](../docs/IMPLEMENTATION.md) for technical details
