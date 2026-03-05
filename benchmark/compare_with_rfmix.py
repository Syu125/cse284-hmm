from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from metrics import align_predictions, compare_sample


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare your HMM ancestry calls against RFMix calls on overlapping SNP positions."
        )
    )
    parser.add_argument("--model", required=True, help="CSV with your model predictions")
    parser.add_argument("--rfmix", required=True, help="CSV with RFMix predictions")
    parser.add_argument("--sample-col", default="sample_id", help="Sample column present in both files")
    parser.add_argument("--position-col", default="position", help="Position column present in both files")
    parser.add_argument("--model-label-col", default="label", help="Label column in your model CSV")
    parser.add_argument("--rfmix-label-col", default="label", help="Label column in RFMix CSV")
    parser.add_argument(
        "--out",
        default="rfmix_comparison_summary.csv",
        help="Output CSV path for per-sample metrics",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    model_df = pd.read_csv(args.model)
    rfmix_df = pd.read_csv(args.rfmix)

    aligned = align_predictions(
        model_df=model_df,
        ref_df=rfmix_df,
        sample_col=args.sample_col,
        pos_col=args.position_col,
        model_label_col=args.model_label_col,
        ref_label_col=args.rfmix_label_col,
    )

    if aligned.empty:
        raise ValueError(
            "No overlapping rows after alignment. Confirm sample/position columns and preprocessing match."
        )

    per_sample_rows = []
    for sample_id, group in aligned.groupby("sample_id", sort=True):
        summary = compare_sample(group, model_col="model_label", ref_col="ref_label")
        per_sample_rows.append(
            {
                "sample_id": sample_id,
                "n_sites": summary.n_sites,
                "concordance": summary.concordance,
                "cohen_kappa": summary.cohen_kappa,
                "switches_per_mb_model": summary.switches_per_mb_model,
                "switches_per_mb_rfmix": summary.switches_per_mb_ref,
                "median_tract_bp_model": summary.median_tract_bp_model,
                "median_tract_bp_rfmix": summary.median_tract_bp_ref,
                "yri_pct_model": summary.yri_pct_model,
                "yri_pct_rfmix": summary.yri_pct_ref,
            }
        )

    report_df = pd.DataFrame(per_sample_rows).sort_values("sample_id").reset_index(drop=True)
    report_df.to_csv(args.out, index=False)

    print(f"Aligned rows: {len(aligned):,}")
    print(f"Samples compared: {report_df['sample_id'].nunique()}")
    print(f"Mean concordance: {report_df['concordance'].mean():.4f}")
    print(f"Mean kappa: {report_df['cohen_kappa'].mean():.4f}")
    print(f"Saved summary: {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
