from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_int_list(value: str) -> list[int]:
    values = [part.strip() for part in value.split(",") if part.strip()]
    if not values:
        raise ValueError("Expected at least one integer value")
    return [int(part) for part in values]


def run_command(command: list[str]) -> None:
    completed = subprocess.run(command, text=True)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed with code {completed.returncode}: {' '.join(command)}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run repeated benchmark evaluations across sample sizes and random seeds."
    )
    parser.add_argument("--sample-sizes", default="5,20,50", help="Comma-separated sample sizes")
    parser.add_argument("--seeds", default="0,1,2,3,4", help="Comma-separated random seeds")
    parser.add_argument("--sample-strategy", choices=["first", "random"], default="random")

    parser.add_argument("--vcf", default="data/processed/chr22_slice.vcf.gz")
    parser.add_argument("--panel", default="data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel")
    parser.add_argument("--map", default="data/raw/maps/genetic_map_GRCh37_chr22.txt")
    parser.add_argument("--query-pop", default="ASW")
    parser.add_argument("--generations", type=float, default=100.0)

    parser.add_argument("--rfmix", default="benchmark/predictions/rfmix_predictions.csv")
    parser.add_argument("--valid-labels", default="YRI,CEU,HET")

    parser.add_argument("--predictions-dir", default="benchmark/predictions")
    parser.add_argument("--results-dir", default="benchmark/results")
    parser.add_argument("--out-prefix", default="sample_sweep")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    sample_sizes = parse_int_list(args.sample_sizes)
    seeds = parse_int_list(args.seeds)

    predictions_dir = Path(args.predictions_dir)
    results_dir = Path(args.results_dir)
    predictions_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    python_exe = sys.executable
    run_rows: list[dict[str, float | int | str]] = []

    for sample_size in sample_sizes:
        for seed in seeds:
            model_out = predictions_dir / f"model_predictions_n{sample_size}_seed{seed}.csv"
            compare_out = results_dir / f"comparison_n{sample_size}_seed{seed}.csv"

            export_cmd = [
                python_exe,
                "benchmark/export_model_predictions.py",
                "--vcf",
                args.vcf,
                "--panel",
                args.panel,
                "--map",
                args.map,
                "--query-pop",
                args.query_pop,
                "--limit-samples",
                str(sample_size),
                "--sample-strategy",
                args.sample_strategy,
                "--seed",
                str(seed),
                "--generations",
                str(args.generations),
                "--out",
                str(model_out),
            ]
            print(f"[RUN] export n={sample_size}, seed={seed}")
            run_command(export_cmd)

            compare_cmd = [
                python_exe,
                "benchmark/compare_with_rfmix.py",
                "--model",
                str(model_out),
                "--rfmix",
                args.rfmix,
                "--valid-labels",
                args.valid_labels,
                "--out",
                str(compare_out),
            ]
            print(f"[RUN] compare n={sample_size}, seed={seed}")
            run_command(compare_cmd)

            compare_df = pd.read_csv(compare_out)
            run_rows.append(
                {
                    "sample_size": sample_size,
                    "seed": seed,
                    "aligned_rows": int(compare_df["n_sites"].sum()),
                    "samples_compared": int(compare_df["sample_id"].nunique()),
                    "mean_concordance": float(compare_df["concordance"].mean()),
                    "mean_kappa": float(compare_df["cohen_kappa"].mean()),
                }
            )

    runs_df = pd.DataFrame(run_rows).sort_values(["sample_size", "seed"]).reset_index(drop=True)

    summary_df = (
        runs_df.groupby("sample_size", as_index=False)
        .agg(
            repeats=("seed", "count"),
            mean_aligned_rows=("aligned_rows", "mean"),
            mean_samples_compared=("samples_compared", "mean"),
            mean_concordance=("mean_concordance", "mean"),
            std_concordance=("mean_concordance", "std"),
            mean_kappa=("mean_kappa", "mean"),
            std_kappa=("mean_kappa", "std"),
        )
        .sort_values("sample_size")
        .reset_index(drop=True)
    )

    runs_out = results_dir / f"{args.out_prefix}_runs.csv"
    summary_out = results_dir / f"{args.out_prefix}_summary.csv"

    runs_df.to_csv(runs_out, index=False)
    summary_df.to_csv(summary_out, index=False)

    print(f"Saved per-run results: {runs_out.resolve()}")
    print(f"Saved summary results: {summary_out.resolve()}")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
