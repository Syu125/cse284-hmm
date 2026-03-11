# Model Benchmarking: RFMix & FLARE Comparison

This folder contains tools and results for evaluating and comparing the HMM model against external local-ancestry methods (RFMIX).

## Prepping the Benchmark

### First, let's prepare the benchmark data

For this benchmark, I first compared my HMM implementation to FLARE because both are HMM-based local ancestry methods. However, they produced very different outputs. When taking a deeper dive, I noticed it was because of how the implementation:
- My implementation infers diploid states directly from population allels frequencies with a simple transition structure.
- FLARE's implementation is haplotype-oriented and uses EM parameters estimations to fine-tune the transitions.
As a result of this difference, when comparing my implementation to FLARE, the similarity scores were very low.

So, I decided to test with RFMix, which is a more advanced ("up-to-date") tool for LAI. I compared both implementations (my HMM and FLARE) to RFMix, and the results show that my HMM was slightly more aligned with RFMix's results.

Below, I walk through setting up both FLARE and RFmix, the steps to run the benchmark, and the results.

## Running the Benchmark

This command runs everything:
```bash
bash benchmark/run_benchmark.sh all
```

If you only need to run a specific stage, use the following commands:
```bash
bash benchmark/run_benchmark.sh prepare
bash benchmark/run_benchmark.sh flare
bash benchmark/run_benchmark.sh rfmix
bash benchmark/run_benchmark.sh convert
bash benchmark/run_benchmark.sh model
bash benchmark/run_benchmark.sh sweep
bash benchmark/run_benchmark.sh format
```

### Stage Mapping

- `prepare`: split query/reference files and build tool-specific maps/panel files
- `flare`: run FLARE and write raw outputs to `benchmark/results/`
- `rfmix`: run RFMix and write raw outputs to `benchmark/results/`
- `convert`: convert FLARE/RFMix raw outputs to SNP-level CSVs under `benchmark/predictions/`
- `model`: export HMM SNP-level predictions under `benchmark/predictions/model/`
- `sweep`: run sample-size sweeps and pairwise comparisons
- `format`: generate summary tables and plots

## Benchmark Results

The following results use full 3-class evaluation (`YRI`, `CEU`, `HET`) on the chr22 slice.

### Sample-Size Sweep (5 Random Repeats per N)

| Method Pair | Sample Size | Repeats | Mean Aligned Rows | Mean Samples Compared | Mean Concordance | Std Concordance | Mean Kappa | Std Kappa |
|------------|------------:|--------:|------------------:|----------------------:|-----------------:|----------------:|-----------:|----------:|
| `hmm_vs_rfmix` | 5  | 5 | 79,730.0 | 5.0  | 0.795023 | 0.080396 | 0.542547 | 0.215686 |
| `hmm_vs_rfmix` | 20 | 5 | 318,920.0 | 20.0 | 0.810998 | 0.058768 | 0.502608 | 0.101310 |
| `hmm_vs_rfmix` | 50 | 5 | 797,300.0 | 50.0 | 0.820949 | 0.010196 | 0.483218 | 0.013899 |
| `flare_vs_rfmix` | 5  | 5 | 3,880.0  | 5.0  | 0.731082 | 0.113002 | 0.334463 | 0.244155 |
| `flare_vs_rfmix` | 20 | 5 | 15,520.0 | 20.0 | 0.781108 | 0.033241 | 0.454204 | 0.068581 |
| `flare_vs_rfmix` | 50 | 5 | 38,800.0 | 50.0 | 0.794907 | 0.006684 | 0.490487 | 0.017056 |
| `hmm_vs_flare` | 5  | 5 | 3,585.0  | 5.0  | 0.622817 | 0.076610 | 0.122853 | 0.170209 |
| `hmm_vs_flare` | 20 | 5 | 14,340.0 | 20.0 | 0.703068 | 0.054431 | 0.210507 | 0.102977 |
| `hmm_vs_flare` | 50 | 5 | 35,850.0 | 50.0 | 0.720123 | 0.008110 | 0.236270 | 0.011355 |

As you can see, the variance drops with larger sample sizes across all method pairs, with the strongest stability at `N=50`.  
`hmm_vs_rfmix` remains the top pair on both concordance and kappa, while `hmm_vs_flare` remains the lowest.

### Concordance And Kappa Trends

<p align="center">
	<img src="plots/sample_sweep_concordance.png" alt="Sample Sweep Concordance" width="49%" />
	<img src="plots/sample_sweep_kappa.png" alt="Sample Sweep Kappa" width="49%" />
</p>

Concordance (how similar the sets agree):
- `hmm_vs_rfmix` has the highest concordance across all sample sizes, indicating the closest overall agreement with RFMix among the evaluated pairs.
- `flare_vs_rfmix` is intermediate and improves as sample size increases.
- `hmm_vs_flare` remains the lowest, showing that the two HMM-family methods produce systematically different calls despite both being HMM-based.

Kappa (how similar the sets agree, accounting for chance):
- `hmm_vs_rfmix` also has the strongest chance-adjusted agreement (highest kappa), reinforcing the concordance trend.
- `flare_vs_rfmix` improves with larger sample sizes and approaches moderate agreement.
- `hmm_vs_flare` stays substantially lower, indicating disagreement beyond class prevalence effects.

## Folder Structure

```
benchmark/
├── BENCHMARK.md                   # Benchmark writeup and usage
├── run_benchmark.sh               # Single benchmark entrypoint
├── scripts/                       # Benchmark Python scripts
├── predictions/                   # Method predictions (my HM, FLARE, and RFMix)
├── results/                       # Raw outputs and comparison files
├── plots/                         # Formatted benchmark plots/tables
├── data/
└── flare.jar
```
