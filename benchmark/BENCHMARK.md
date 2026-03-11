# Model Benchmarking: RFMix & FLARE Comparison

This folder contains tools and results for evaluating and comparing the HMM model against external local-ancestry methods (RFMIX).

## Prepping the Benchmark

### First, let's prepare the benchmark data

```bash
python benchmark/prepare_benchmark_data.py
```

For this benchmark, I first compared my HMM implementation to FLARE because both are HMM-based local ancestry methods. However, they produced very different outputs. When taking a deeper dive, I noticed it was because of how the implementation:
- My implementation infers diploid states directly from population allels frequencies with a simple transition structure.
- FLARE's implementation is haplotype-oriented and uses EM parameters estimations to fine-tune the transitions.
As a result of this difference, when comparing my implementation to FLARE, the similarity scores were very low.

So, I decided to test with RFMix, which is a more advanced ("up-to-date") tool for LAI. I compared both implementations (my HMM and FLARE) to RFMix, and the results show that my HMM was slightly more aligned with RFMix's results. 

Below, I walk through setting up both FLARE and RFmix, the steps to run the benchmark, and the results.

### FLARE

**Install FLARE:**
```bash
# Java 11+ required
conda install -c conda-forge "openjdk>=11"

# Download executable jar
wget https://faculty.washington.edu/browning/flare.jar -O benchmark/flare.jar
```

**Create FLARE ref-panel file:**
```bash
python benchmark/make_flare_ref_panel.py \
  --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --reference-vcf benchmark/data/reference_yri_ceu_chr22.vcf.gz \
  --pops YRI,CEU \
  --out benchmark/data/flare_ref_panel_yri_ceu.txt
```

**Create FLARE genetic map:**
```bash
python benchmark/make_flare_map.py \
  --input data/raw/maps/genetic_map_GRCh37_chr22.txt \
  --output benchmark/data/genetic_map_GRCh37_chr22.flare.map
```

**Run FLARE:**
```bash
java -Xmx8g -jar benchmark/flare.jar \
  ref=benchmark/data/reference_yri_ceu_chr22.vcf.gz \
  ref-panel=benchmark/data/flare_ref_panel_yri_ceu.txt \
  gt=benchmark/data/query_asw_chr22.vcf.gz \
  map=benchmark/data/genetic_map_GRCh37_chr22.flare.map \
  out=benchmark/results/flare_output \
  nthreads=4
```

Output: `benchmark/results/flare_output.anc.vcf.gz`

**Convert FLARE Output to SNP-Level CSV:**
```bash
tabix -p vcf benchmark/results/flare_output.anc.vcf.gz

python benchmark/convert_flare_to_snp_csv.py \
  --flare benchmark/results/flare_output.anc.vcf.gz \
  --tie-policy unknown \
  --out benchmark/predictions/flare_predictions.csv \
  --chromosome 22
```

### RFMix

**Install RFMix:**
```bash
conda install -c bioconda rfmix
```

**Create RFMix genetic map:**
```bash
python benchmark/prepare_benchmark_data.py
python benchmark/fix_genetic_map.py
```

**Run RFMix:**
```bash
rfmix \
  -f benchmark/data/query_asw_chr22.vcf.gz \
  -r benchmark/data/reference_yri_ceu_chr22.vcf.gz \
  -g benchmark/data/genetic_map_GRCh37_chr22.strict.txt \
  -m data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --chromosome=22 \
  -o benchmark/results/rfmix_output
```

**Output files:**
- `benchmark/results/rfmix_output.msp.tsv` - Ancestry state per segment
- `benchmark/results/rfmix_output.sis.tsv` - Segment info (sample, start, end)
- `benchmark/results/rfmix_output.rfmix.Q` - Posterior probabilities
- `benchmark/results/rfmix_output.fb.tsv` - Forward-backward details

**Convert RFMix Output to SNP-Level CSV:**
```bash
python benchmark/convert_rfmix_to_snp_csv.py \
  --msp benchmark/results/rfmix_output.msp.tsv \
  --vcf benchmark/data/query_asw_chr22.vcf.gz \
  --tie-policy unknown \
  --out benchmark/predictions/rfmix_predictions.csv \
  --chromosome 22
```
## Running the Benchmark
### Export My HMM Predictions

```bash
python benchmark/export_model_predictions.py \
  --vcf data/processed/chr22_slice.vcf.gz \
  --panel data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel \
  --map data/raw/maps/genetic_map_GRCh37_chr22.txt \
  --query-pop ASW \
  --out benchmark/predictions/model_predictions.csv
```

### 4. Sample-Size Stability Sweep

The sweep now reports all three pairwise comparisons per run:
- `hmm_vs_rfmix`
- `hmm_vs_flare`
- `flare_vs_rfmix`

To estimate variability across sample composition, run repeated seeded subsamples:

```bash
python benchmark/run_sample_size_sweep.py \
  --sample-sizes 5,20,50 \
  --seeds 0,1,2,3,4 \
  --sample-strategy random \
  --rfmix benchmark/predictions/rfmix_predictions.csv \
  --flare benchmark/predictions/flare_predictions.csv \
  --out-prefix sample_sweep
```

Plots and tables have already been generated (can be found in the [plots](/benchmark/plots/) folder). For reference, below is the command to generate them:

```bash
python benchmark/format_sweep_table.py \
  --input benchmark/results/sample_sweep_summary.csv \
  --out-csv benchmark/results/sample_sweep_summary_wide.csv \
```

## Benchmark Results

The following results use full 3-class evaluation (`YRI`, `CEU`, `HET`) on the chr22 slice.

### Sample-Size Sweep (5 Random Repeats per N)

| Sample Size | Repeats | Mean Aligned Rows | Mean Samples Compared | Mean Concordance | Std Concordance | Mean Kappa | Std Kappa |
|------------:|--------:|------------------:|----------------------:|-----------------:|----------------:|-----------:|----------:|
| 5           | 5       | 79,730.0          | 5.0                   | 0.795023         | 0.080396        | 0.542547   | 0.215686  |
| 20          | 5       | 318,920.0         | 20.0                  | 0.810998         | 0.058768        | 0.502608   | 0.101310  |
| 50          | 5       | 797,300.0         | 50.0                  | 0.820949         | 0.010196        | 0.483218   | 0.013899  |

Interpretation:
- Variance drops substantially as sample size increases (especially at `N=50`), indicating more stable benchmark estimates.
- Kappa remains in a similar range across sample sizes, suggesting consistent chance-adjusted performance.


## Folder Structure

```
benchmark/
├── export_model_predictions.py      # Export HMM predictions to CSV
├── compare_with_reference.py        # Compare against any SNP-level reference CSV
├── convert_rfmix_to_snp_csv.py      # Convert RFMix output format
├── make_flare_ref_panel.py          # Create FLARE 2-column ref-panel file
├── make_flare_map.py                # Create FLARE PLINK-format map file
├── convert_flare_to_snp_csv.py      # Convert FLARE .anc.vcf.gz to SNP-level CSV
├── metrics.py                       # Evaluation metrics  
├── __init__.py
├── predictions/                     # HMM model predictions
│   └── .gitkeep
└── results/                         # Evaluation results
    └── .gitkeep
```
