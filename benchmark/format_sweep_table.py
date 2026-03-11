#!/usr/bin/env python3
"""Format sample sweep summary into publication-ready tables and plots.

Input expected: benchmark/results/*_summary.csv from run_sample_size_sweep.py
(with method_pair, sample_size, repeats, mean/std columns).

Outputs:
- Wide CSV with one row per sample size and pairwise metrics as columns
- Pivoted CSV tables for concordance and kappa (mean and std)
- PNG plots for concordance, kappa, and aligned rows vs sample size
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


METHOD_ORDER = ["hmm_vs_rfmix", "hmm_vs_flare", "flare_vs_rfmix"]
METHOD_PAIR_COLORS = {
    "hmm_vs_rfmix": "#829eeb",  # HMM
    "hmm_vs_flare": "#93c4e0",  # FLARE
    "flare_vs_rfmix": "#9ce3e9",  # RFMix
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Format sample sweep summary into tables + plots")
    parser.add_argument(
        "--input",
        default="benchmark/results/sample_sweep_summary.csv",
        help="Input summary CSV from run_sample_size_sweep.py",
    )
    parser.add_argument(
        "--out-csv",
        default="benchmark/plots/sample_sweep_summary_wide.csv",
        help="Output wide-format CSV",
    )
    parser.add_argument(
        "--out-dir",
        default="benchmark/plots",
        help="Output directory for generated tables and plots",
    )
    return parser


def build_wide_table(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []

    for sample_size, group in df.groupby("sample_size", sort=True):
        row: dict[str, float | int] = {
            "sample_size": int(sample_size),
            "repeats": int(group["repeats"].max()),
        }

        for method in METHOD_ORDER:
            method_rows = group[group["method_pair"] == method]
            if method_rows.empty:
                continue
            r = method_rows.iloc[0]
            row[f"{method}_aligned_rows"] = float(r["mean_aligned_rows"])
            row[f"{method}_concordance_mean"] = float(r["mean_concordance"])
            row[f"{method}_concordance_std"] = float(r["std_concordance"])
            row[f"{method}_kappa_mean"] = float(r["mean_kappa"])
            row[f"{method}_kappa_std"] = float(r["std_kappa"])

        rows.append(row)

    return pd.DataFrame(rows).sort_values("sample_size").reset_index(drop=True)


def save_metric_tables(df: pd.DataFrame, out_dir: Path) -> list[Path]:
    saved: list[Path] = []

    concordance_mean = df.pivot(index="sample_size", columns="method_pair", values="mean_concordance")
    concordance_std = df.pivot(index="sample_size", columns="method_pair", values="std_concordance")
    kappa_mean = df.pivot(index="sample_size", columns="method_pair", values="mean_kappa")
    kappa_std = df.pivot(index="sample_size", columns="method_pair", values="std_kappa")

    for metric_df, name in [
        (concordance_mean, "sample_sweep_concordance_mean.csv"),
        (concordance_std, "sample_sweep_concordance_std.csv"),
        (kappa_mean, "sample_sweep_kappa_mean.csv"),
        (kappa_std, "sample_sweep_kappa_std.csv"),
    ]:
        ordered_cols = [c for c in METHOD_ORDER if c in metric_df.columns]
        metric_df = metric_df[ordered_cols].reset_index()
        out_path = out_dir / name
        metric_df.to_csv(out_path, index=False)
        saved.append(out_path)

    return saved


def plot_metric_with_error(df: pd.DataFrame, mean_col: str, std_col: str, ylabel: str, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))

    for method in METHOD_ORDER:
        method_df = df[df["method_pair"] == method].sort_values("sample_size")
        if method_df.empty:
            continue
        plt.errorbar(
            method_df["sample_size"],
            method_df[mean_col],
            yerr=method_df[std_col],
            marker="o",
            capsize=3,
            linewidth=1.8,
            color=METHOD_PAIR_COLORS.get(method),
            label=method,
        )

    plt.xlabel("Sample Size")
    plt.ylabel(ylabel)
    plt.title(f"{ylabel} Across Sample Sizes")
    plt.grid(alpha=0.3)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def plot_aligned_rows(df: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(8, 5))

    for method in METHOD_ORDER:
        method_df = df[df["method_pair"] == method].sort_values("sample_size")
        if method_df.empty:
            continue
        plt.plot(
            method_df["sample_size"],
            method_df["mean_aligned_rows"],
            marker="o",
            linewidth=1.8,
            color=METHOD_PAIR_COLORS.get(method),
            label=method,
        )

    plt.xlabel("Sample Size")
    plt.ylabel("Mean Aligned Rows")
    plt.title("Aligned Rows Across Sample Sizes")
    plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()


def main() -> None:
    args = build_parser().parse_args()

    input_path = Path(args.input)
    df = pd.read_csv(input_path)

    required = {
        "method_pair",
        "sample_size",
        "repeats",
        "mean_aligned_rows",
        "mean_concordance",
        "std_concordance",
        "mean_kappa",
        "std_kappa",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Input file missing required columns: {sorted(missing)}")

    wide = build_wide_table(df)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    wide.to_csv(out_csv, index=False)
    print(f"Saved wide CSV: {out_csv.resolve()}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    table_paths = save_metric_tables(df, out_dir)
    for path in table_paths:
        print(f"Saved table: {path.resolve()}")

    concordance_plot = out_dir / "sample_sweep_concordance.png"
    kappa_plot = out_dir / "sample_sweep_kappa.png"
    rows_plot = out_dir / "sample_sweep_aligned_rows.png"

    plot_metric_with_error(df, "mean_concordance", "std_concordance", "Concordance", concordance_plot)
    plot_metric_with_error(df, "mean_kappa", "std_kappa", "Cohen Kappa", kappa_plot)
    plot_aligned_rows(df, rows_plot)

    print(f"Saved plot: {concordance_plot.resolve()}")
    print(f"Saved plot: {kappa_plot.resolve()}")
    print(f"Saved plot: {rows_plot.resolve()}")


if __name__ == "__main__":
    main()
