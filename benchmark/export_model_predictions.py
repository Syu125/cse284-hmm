from __future__ import annotations

import argparse
from pathlib import Path
import sys

import pandas as pd
import pysam

SRC_ROOT = Path(__file__).resolve().parent.parent / 'src'
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from data.data_parser import (
    get_allele_frequencies,
    get_genetic_map,
    get_population_dict,
    interpolate_genetic_position,
)
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export per-SNP ancestry predictions from your HMM model for downstream benchmarking."
    )
    parser.add_argument(
        "--vcf",
        default="../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
    )
    parser.add_argument("--panel", default="../data/integrated_call_samples_v3.20130502.ALL.panel")
    parser.add_argument("--map", default="../data/genetic_map_GRCh37_chr22.txt")
    parser.add_argument(
        "--query-pop",
        default="ASW",
        help="Population to infer for and export (for example ASW)",
    )
    parser.add_argument("--limit-samples", type=int, default=5)
    parser.add_argument("--generations", type=float, default=100.0)
    parser.add_argument("--out", default="model_predictions_chr22.csv")
    return parser


def main() -> None:
    args = build_parser().parse_args()

    pops = get_population_dict(args.panel)
    if args.query_pop not in pops or len(pops[args.query_pop]) == 0:
        raise ValueError(f"Population {args.query_pop} not found or empty in panel file")

    phys, gen = get_genetic_map(args.map)
    map_min, map_max = int(phys[0]), int(phys[-1])

    yri_freqs = get_allele_frequencies(args.vcf, pops["YRI"])
    ceu_freqs = get_allele_frequencies(args.vcf, pops["CEU"])
    valid_snps = sorted(
        pos
        for pos in set(yri_freqs.keys()) & set(ceu_freqs.keys())
        if map_min <= pos <= map_max
    )
    valid_set = set(valid_snps)

    emission = EmissionModel(yri_freqs, ceu_freqs)
    transition = TransitionModel(generations=args.generations)
    engine = InferenceEngine(emission, transition)

    # Get actual samples from the VCF file
    vcf_check = pysam.VariantFile(args.vcf)
    available_samples = list(vcf_check.header.samples)
    vcf_check.close()
    
    # Filter to only samples that exist in both panel and VCF
    candidate_samples = pops[args.query_pop][: args.limit_samples]
    selected_samples = [s for s in candidate_samples if s in available_samples]
    
    if not selected_samples:
        raise ValueError(f"No {args.query_pop} samples found in VCF. Available: {available_samples[:5]}")
    
    print(f"Processing {len(selected_samples)} {args.query_pop} samples from VCF")
    rows = []

    for sample_id in selected_samples:
        vcf = pysam.VariantFile(args.vcf)
        snp_positions = []
        genotypes = []

        # Try to fetch chromosome first; fall back to iterating all records if no index
        try:
            records = vcf.fetch("22")
        except ValueError:
            # No index available, iterate all records
            records = vcf

        for record in records:
            if record.pos in valid_set:
                snp_positions.append(record.pos)
                genotypes.append(record.samples[sample_id].allele_indices)

        get_cm = lambda position: interpolate_genetic_position(position, phys, gen)
        states = engine.run_viterbi(snp_positions, genotypes, get_cm)

        rows.extend(
            {
                "sample_id": sample_id,
                "position": pos,
                "label": state,
            }
            for pos, state in zip(snp_positions, states)
        )

        yri_pct = 100.0 * sum(state == "YRI" for state in states) / len(states)
        print(f"Processed {sample_id}: {len(states):,} SNPs, YRI={yri_pct:.2f}%")

    output_path = Path(args.out)
    pd.DataFrame(rows).to_csv(output_path, index=False)
    print(f"Saved model predictions to {output_path.resolve()}")


if __name__ == "__main__":
    main()
