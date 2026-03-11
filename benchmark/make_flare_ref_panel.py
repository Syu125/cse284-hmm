#!/usr/bin/env python3
"""Create a FLARE-compatible 2-column reference panel file.

FLARE expects `ref-panel` to have exactly two whitespace-delimited fields per line:
1) sample ID (present in reference VCF)
2) ancestry label

This script filters the project's 1000G panel file to selected populations and
optionally intersects with samples present in a reference VCF.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


def read_vcf_samples(vcf_path: str) -> list[str]:
    """Return sample IDs from the #CHROM header in a gzipped VCF."""
    with gzip.open(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                return fields[9:]
    raise ValueError(f"Could not find #CHROM header in VCF: {vcf_path}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create FLARE ref-panel file with two columns: sample and ancestry label"
    )
    parser.add_argument(
        "--panel",
        default="data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel",
        help="Input 1000G panel TSV with columns including sample and pop",
    )
    parser.add_argument(
        "--reference-vcf",
        default="benchmark/data/reference_yri_ceu_chr22.vcf.gz",
        help="Reference VCF used by FLARE; if provided, output is restricted to its samples",
    )
    parser.add_argument(
        "--pops",
        default="YRI,CEU",
        help="Comma-separated populations to include as ancestry labels (default: YRI,CEU)",
    )
    parser.add_argument(
        "--out",
        default="benchmark/data/flare_ref_panel_yri_ceu.txt",
        help="Output path for FLARE ref-panel file",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    requested_pops = [p.strip() for p in args.pops.split(",") if p.strip()]
    if not requested_pops:
        raise ValueError("--pops must include at least one population")

    panel_df = pd.read_csv(args.panel, sep="\t")
    required_cols = {"sample", "pop"}
    missing = required_cols.difference(panel_df.columns)
    if missing:
        raise ValueError(f"Panel file missing required columns: {sorted(missing)}")

    filtered = panel_df[panel_df["pop"].isin(requested_pops)].copy()
    if filtered.empty:
        raise ValueError(
            f"No rows found for requested populations: {','.join(requested_pops)}"
        )

    # Keep only populations present in the reference VCF to prevent FLARE input mismatch.
    if args.reference_vcf:
        vcf_samples = set(read_vcf_samples(args.reference_vcf))
        filtered = filtered[filtered["sample"].isin(vcf_samples)].copy()

    if filtered.empty:
        raise ValueError("No overlapping samples after VCF intersection/filtering")

    filtered = filtered[["sample", "pop"]].drop_duplicates().sort_values("sample")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Delimiter is a single tab; FLARE accepts any whitespace delimiter.
    filtered.to_csv(out_path, sep="\t", header=False, index=False)

    print(f"Requested populations: {','.join(requested_pops)}")
    print(f"Wrote {len(filtered):,} rows to {out_path.resolve()}")


if __name__ == "__main__":
    main()
