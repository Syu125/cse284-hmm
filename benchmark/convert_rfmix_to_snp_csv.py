#!/usr/bin/env python3
"""Convert RFMix .msp.tsv segment output to per-SNP ancestry CSV."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import pysam


def convert_msp_to_snp_level(
    msp_file: str, vcf_file: str, output_file: str, chromosome: str = "22"
) -> None:
    """
    Expand RFMix segment-level haplotype predictions to per-SNP sample predictions.

    RFMix .msp.tsv format:
    - Line 1: subpopulation codes (CEU=0, YRI=1)
    - Line 2: header with columns: #chm spos epos sgpos egpos n snps SAMPLE.0 SAMPLE.1 ...
    - Data rows: one per segment with haplotype calls (0 or 1)

    This function:
    1. Parses haplotype-level calls
    2. For each segment, finds SNPs in that range
    3. Assigns consensus ancestry (majority of two haplotypes) to sample
    4. Writes SNP-level CSV: sample_id, position, label
    """
    print(f"[-] Reading RFMix segments from {msp_file}...")

    # Manually parse header and data
    with open(msp_file) as f:
        # Skip metadata line
        f.readline()
        # Read header line
        header_line = f.readline().strip()

    # Parse header - columns are tab-separated
    header = header_line.split("\t")
    
    # Remove leading # from first column if present
    header[0] = header[0].lstrip("#")
    
    print(f"    Header columns: {header[:10]}...")  # Show first 10
    
    # Find sample mapping from haplotype columns (skip first 6 cols: chm, spos, epos, sgpos, egpos, n snps)
    haplotype_cols = header[6:]
    samples = {}
    for col in haplotype_cols:
        sample_id = col.rsplit(".", 1)[0]  # Remove .0 or .1 suffix
        if sample_id not in samples:
            samples[sample_id] = []
        samples[sample_id].append(col)
    
    print(f"    Found {len(samples)} samples with haplotype data")

    # Read segment data using pandas
    # Skip first 2 lines, don't use pandas header parsing
    msp_df = pd.read_csv(msp_file, sep="\t", skiprows=2, header=None)
    
    # Manually assign column names from the header we parsed
    msp_df.columns = header
    
    print(f"    DataFrame columns: {list(msp_df.columns[:10])}")
    print(f"    Processed {len(msp_df)} segments")

    # Load VCF and extract all SNP positions
    print(f"[-] Reading SNP positions from {vcf_file}...")
    snp_positions = []
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf.fetch(chromosome):
            snp_positions.append(record.pos)

    snp_positions = sorted(set(snp_positions))
    print(f"    Found {len(snp_positions):,} unique SNP positions")

    # Ancestry mappings: 0=CEU, 1=YRI (from RFMix header)
    ancestry_map = {0: "CEU", 1: "YRI"}

    # Expand segments to SNPs
    snp_rows = []
    for _, segment in msp_df.iterrows():
        start_pos = int(segment["spos"])
        end_pos = int(segment["epos"])

        # Find SNPs in this segment
        snps_in_segment = [pos for pos in snp_positions if start_pos <= pos < end_pos]

        # For each sample, get consensus ancestry call
        for sample_id, hap_cols in samples.items():
            # Get haplotype calls (0 or 1) for this sample in this segment
            hap_calls = []
            for hap_col in hap_cols:
                if hap_col in msp_df.columns:
                    try:
                        call = int(segment[hap_col])
                        hap_calls.append(call)
                    except (ValueError, TypeError):
                        continue

            # Consensus: if both haplotypes match, use that; else use majority
            if len(hap_calls) == 2:
                consensus = int(round(sum(hap_calls) / len(hap_calls)))
            elif len(hap_calls) == 1:
                consensus = hap_calls[0]
            else:
                continue  # Skip if no valid calls

            ancestry = ancestry_map.get(consensus, "unknown")

            # Assign this ancestry to all SNPs in segment
            for snp_pos in snps_in_segment:
                snp_rows.append({"sample_id": sample_id, "position": snp_pos, "label": ancestry})

    output_df = pd.DataFrame(snp_rows).sort_values(["sample_id", "position"]).reset_index(drop=True)
    output_df.to_csv(output_file, index=False)

    print(f"[+] Converted {len(snp_rows):,} SNPs across {output_df['sample_id'].nunique()} samples")
    print(f"[+] Saved to {Path(output_file).resolve()}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert RFMix .msp.tsv to per-SNP CSV")
    parser.add_argument("--msp", required=True, help="RFMix .msp.tsv file")
    parser.add_argument("--vcf", required=True, help="VCF file (to get SNP positions)")
    parser.add_argument("--out", default="rfmix_snp_level.csv", help="Output CSV path")
    parser.add_argument("--chromosome", default="22", help="Chromosome identifier")
    return parser


if __name__ == "__main__":
    args = build_parser().parse_args()
    convert_msp_to_snp_level(args.msp, args.vcf, args.out, args.chromosome)
