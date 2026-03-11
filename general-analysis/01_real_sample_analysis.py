"""
Analyze real admixed individuals from the 1000 Genomes dataset.

This script:
1. Loads ASW (African ancestry in Southwest USA) samples
2. Infers local ancestry for each individual using the HMM
3. Prints summary statistics for each person
4. Generates ancestry paintings for each individual

ASW samples are known to have admixed ancestry, making them suitable
for validating ancestry inference on real data.
"""

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent.parent / 'src'
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

OUTPUT_DIR = Path("general-analysis/output")
OUTPUT_DIR.mkdir(exist_ok=True)

import pysam
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization import plot_ancestry

def main():
    # 1. Setup - Auto-detect available dataset
    slice_path = "data/processed/chr22_slice.vcf.gz"
    full_path = "data/raw/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    
    if Path(slice_path).exists():
        vcf_path = slice_path
        print(f"[+] Using slice dataset: {slice_path}")
    elif Path(full_path).exists():
        vcf_path = full_path
        print(f"[+] Using full dataset: {full_path}")
    else:
        print("[!] Error: No VCF data found.")
        print("    Run: bash data/prepare_data.sh slice")
        sys.exit(1)
    
    panel_path = "data/raw/panels/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "data/raw/maps/genetic_map_GRCh37_chr22.txt"
    
    pops = get_population_dict(panel_path)
    phys, gen = get_genetic_map(map_path)
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # 2. Pick some admixed sample (ASW)
    asw_sample_ids = pops['ASW'][:3] # Take the first person in the ASW group
    print(f"[-] Analyzing real samples: {asw_sample_ids}")
    
    # 3. Collect their data
    for sample_id in asw_sample_ids:
        print(f"\n[!] Processing individual: {sample_id}")
        
        # 3. Collect data for THIS specific sample
        vcf = pysam.VariantFile(vcf_path)
        common_pos = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
        common_set = set(common_pos)
        
        snp_positions = []
        genotypes = []
        
        for record in vcf:
            if record.pos in common_set:
                sample_data = record.samples[sample_id]
                gt = sample_data.allele_indices
                if (not sample_data.phased) or gt is None or len(gt) < 2 or any(a is None for a in gt):
                    continue
                snp_positions.append(record.pos)
                genotypes.append(gt)

        if not snp_positions:
            print("    Warning: no phased SNPs available for this sample; skipping")
            continue

        # 4. Run Inference
        emissions = EmissionModel(yri_freqs, ceu_freqs)
        transitions = TransitionModel(generations=100)
        engine = InferenceEngine(emissions, transitions)
        
        get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
        states = engine.run_viterbi(snp_positions, genotypes, get_cm)
        
        # 5. Global Stats
        total = len(states)
        yri_yri_pct = (states.count("YRI_YRI") / total) * 100
        yri_ceu_pct = (states.count("YRI_CEU") / total) * 100
        ceu_yri_pct = (states.count("CEU_YRI") / total) * 100
        ceu_ceu_pct = (states.count("CEU_CEU") / total) * 100
        print(
            "    Summary: "
            f"{yri_yri_pct:.1f}% YRI_YRI | "
            f"{yri_ceu_pct:.1f}% YRI_CEU | "
            f"{ceu_yri_pct:.1f}% CEU_YRI | "
            f"{ceu_ceu_pct:.1f}% CEU_CEU"
        )

        # 6. Visualize - ensure unique filename per sample!
        plot_ancestry(snp_positions, states, save_path=OUTPUT_DIR / f"ancestry_{sample_id}.png")

if __name__ == "__main__":
    main()