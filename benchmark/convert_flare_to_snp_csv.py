#!/usr/bin/env python3
"""Convert FLARE .anc.vcf.gz output to benchmark SNP-level CSV.

Reads AN1/AN2 per sample at each variant and maps ancestry indices to labels using
##ANCESTRY metadata from the VCF header. Output schema matches benchmark tools:
- sample_id
- position
- label

Label resolution from haplotypes:
- same ancestry on both haplotypes -> ancestry label (for example CEU, YRI)
- mixed ancestries -> according to --tie-policy (default: HET)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import pysam


def resolve_consensus_label(an1: int, an2: int, idx_to_name: dict[int, str], tie_policy: str) -> str | None:
    """Resolve diploid hap ancestry calls to one benchmark label."""
    label1 = idx_to_name.get(an1)
    label2 = idx_to_name.get(an2)

    if label1 is None or label2 is None:
        return None

    if label1 == label2:
        return label1

    if tie_policy == "drop":
        return None
    if tie_policy == "ceu":
        return "CEU"
    if tie_policy == "yri":
        return "YRI"
    return "HET"


def parse_ancestry_header(header_record: str) -> dict[int, str]:
    """Parse FLARE ancestry metadata into index -> ancestry name.

    Supported header forms:
    - ##ANCESTRY=<ID=0,Name=YRI>
    - ##ANCESTRY=<YRI=0,CEU=1>
    """
    payload = header_record.split("<", 1)[1].rsplit(">", 1)[0]
    fields = {}
    for item in payload.split(","):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        fields[key.strip()] = value.strip()

    # Form 1: ID/Name pair.
    if "ID" in fields and "Name" in fields:
        idx = int(fields["ID"])
        name = fields["Name"].strip()
        if not name:
            raise ValueError(f"Malformed ancestry header entry: {header_record}")
        return {idx: name}

    # Form 2: ancestry-name -> index pairs (for example YRI=0,CEU=1).
    idx_to_name: dict[int, str] = {}
    for name, idx_str in fields.items():
        try:
            idx = int(idx_str)
        except ValueError:
            continue
        idx_to_name[idx] = name

    if not idx_to_name:
        raise ValueError(f"Malformed ancestry header entry: {header_record}")
    return idx_to_name


def load_ancestry_map(vcf: pysam.VariantFile) -> dict[int, str]:
    """Return ancestry index -> ancestry name mapping from VCF metadata."""
    idx_to_name: dict[int, str] = {}

    # Access raw header text because ancestry metadata key is tool-specific.
    for line in str(vcf.header).splitlines():
        if line.startswith("##ANCESTRY=<"):
            idx_to_name.update(parse_ancestry_header(line))

    if not idx_to_name:
        raise ValueError("No ##ANCESTRY metadata found in FLARE VCF header")

    return idx_to_name


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert FLARE .anc.vcf.gz to SNP-level CSV")
    parser.add_argument("--flare", required=True, help="FLARE output .anc.vcf.gz path")
    parser.add_argument("--out", default="benchmark/predictions/flare_predictions.csv", help="Output CSV path")
    parser.add_argument("--chromosome", default="22", help="Chromosome to process (default: 22)")
    parser.add_argument(
        "--tie-policy",
        choices=["unknown", "drop", "ceu", "yri"],
        default="unknown",
        help="How to resolve mixed haplotype ancestries (default: unknown -> HET)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    with pysam.VariantFile(args.flare) as vcf:
        idx_to_name = load_ancestry_map(vcf)
        samples = list(vcf.header.samples)

        rows = []
        tie_sites = 0
        dropped_ties = 0

        # Use indexed random access when possible; otherwise fall back to
        # sequential iteration so the script works without a .tbi index.
        try:
            record_iter = vcf.fetch(args.chromosome)
            chrom_filter = None
        except (ValueError, OSError):
            record_iter = iter(vcf)
            chrom_filter = args.chromosome

        for rec in record_iter:
            if chrom_filter is not None:
                if rec.chrom != chrom_filter and rec.chrom != f"chr{chrom_filter}":
                    continue
            pos = int(rec.pos)

            for sample_id in samples:
                sample_data = rec.samples[sample_id]
                an1 = sample_data.get("AN1")
                an2 = sample_data.get("AN2")
                if an1 is None or an2 is None:
                    continue

                # AN fields are scalar integers.
                try:
                    a1 = int(an1)
                    a2 = int(an2)
                except (TypeError, ValueError):
                    continue

                if a1 != a2:
                    tie_sites += 1

                label = resolve_consensus_label(a1, a2, idx_to_name, args.tie_policy)
                if label is None:
                    if a1 != a2:
                        dropped_ties += 1
                    continue

                rows.append({"sample_id": sample_id, "position": pos, "label": label})

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        raise ValueError("No rows written. Check chromosome or FLARE AN1/AN2 fields.")

    out_df = out_df.sort_values(["sample_id", "position"]).reset_index(drop=True)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)

    print("Ancestry index map:")
    for idx in sorted(idx_to_name):
        print(f"  {idx}: {idx_to_name[idx]}")
    print(f"Tie policy: {args.tie_policy}")
    print(f"Mixed haplotype sites observed: {tie_sites:,}")
    if args.tie_policy == "drop":
        print(f"Mixed haplotype sites dropped: {dropped_ties:,}")
    print(f"Wrote {len(out_df):,} rows across {out_df['sample_id'].nunique()} samples")
    print(f"Saved: {out_path.resolve()}")


if __name__ == "__main__":
    main()
