#!/usr/bin/env python3
"""Create an RFMix-compatible map: chromosome, pos, gpos (Morgans)."""

from __future__ import annotations

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert genetic map to RFMix format with strict monotonicity."
    )
    parser.add_argument(
        "--input",
        default="data/raw/maps/genetic_map_GRCh37_chr22.txt",
        help="Input genetic map file path",
    )
    parser.add_argument(
        "--output",
        default="benchmark/data/genetic_map_GRCh37_chr22.strict.txt",
        help="Output path for strict RFMix genetic map",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Write output without the 'chromosome\tpos\tgpos' header",
    )
    return parser.parse_args()


def split_fields(line: str) -> list[str]:
    if "\t" in line:
        return line.rstrip("\n").split("\t")
    return line.rstrip("\n").split()


def main() -> int:
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    kept = 0
    dropped = 0
    previous_gen = -1.0
    previous_phys = -1

    with in_path.open("r", encoding="utf-8") as infile, out_path.open(
        "w", encoding="utf-8", newline="\n"
    ) as outfile:
        header = infile.readline()
        if not header:
            raise ValueError("Input map file is empty")

        header_fields = split_fields(header)
        if len(header_fields) < 4:
            raise ValueError("Invalid header in genetic map file")

        if not args.no_header:
            outfile.write("chromosome\tpos\tgpos\n")

        for line in infile:
            stripped = line.strip()
            if not stripped:
                continue

            fields = split_fields(line)
            if len(fields) < 4:
                continue

            try:
                chromosome = fields[0].replace("chr", "")
                physical_position = int(fields[1])
                genetic_position = float(fields[3]) / 100.0
            except ValueError:
                dropped += 1
                continue

            if physical_position > previous_phys and genetic_position > previous_gen:
                outfile.write(
                    f"{chromosome}\t{physical_position}\t{genetic_position:.10f}\n"
                )
                previous_phys = physical_position
                previous_gen = genetic_position
                kept += 1
            else:
                dropped += 1

    print(f"Wrote strict map: {out_path}")
    print(f"Kept rows: {kept}")
    print(f"Dropped rows: {dropped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
