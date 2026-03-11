#!/usr/bin/env python3
"""Create a FLARE-compatible PLINK map file.

Input expected: genetic map with columns similar to
Chromosome, Position(bp), Rate(cM/Mb), Map(cM)

Output (no header), 4 whitespace-delimited fields per row:
1) chromosome (without 'chr' prefix)
2) marker id (CHROM:POS)
3) genetic distance in cM
4) base-pair position (integer)
"""

from __future__ import annotations

import argparse
from pathlib import Path


def split_fields(line: str) -> list[str]:
    if "\t" in line:
        return line.rstrip("\n").split("\t")
    return line.rstrip("\n").split()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert project genetic map to FLARE PLINK map format")
    parser.add_argument(
        "--input",
        default="data/raw/maps/genetic_map_GRCh37_chr22.txt",
        help="Input map file with columns Chromosome/Position(bp)/Map(cM)",
    )
    parser.add_argument(
        "--output",
        default="benchmark/data/genetic_map_GRCh37_chr22.flare.map",
        help="Output PLINK map file for FLARE",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    skipped = 0

    with in_path.open("r", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8", newline="\n") as fout:
        header = fin.readline()
        if not header:
            raise ValueError("Input map file is empty")

        for line in fin:
            raw = line.strip()
            if not raw:
                continue

            fields = split_fields(line)
            if len(fields) < 4:
                skipped += 1
                continue

            try:
                chrom = fields[0].replace("chr", "")
                pos = int(fields[1])
                cm = float(fields[3])
            except ValueError:
                skipped += 1
                continue

            marker = f"{chrom}:{pos}"
            fout.write(f"{chrom}\t{marker}\t{cm:.6f}\t{pos}\n")
            kept += 1

    print(f"Wrote FLARE map: {out_path.resolve()}")
    print(f"Kept rows: {kept}")
    print(f"Skipped rows: {skipped}")


if __name__ == "__main__":
    main()
