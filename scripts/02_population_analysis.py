"""
Population-level analysis of ancestry in admixed samples.

This script:
1. Loads all ASW (Southwest African ancestry) samples from 1000 Genomes
2. Infers ancestry for each individual in parallel
3. Computes statistics on ancestry composition across the population
4. Generates visualizations of the ancestry distribution

Results are saved to CSV for downstream analysis.
"""

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent.parent / 'src'
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

OUTPUT_DIR = Path("scripts/output")
OUTPUT_DIR.mkdir(exist_ok=True)

import pysam
import pandas as pd
import matplotlib.pyplot as plt
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
def main():
    # 1. Setup Data - Auto-detect available dataset
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
    
    asw_samples = pops.get('ASW', [])
    print(f"[-] Total ASW samples to process: {len(asw_samples)}")

    # 2. Initialize HMM components outside the loop (for speed!)
    emissions = EmissionModel(yri_freqs, ceu_freqs)
    transitions = TransitionModel(generations=10)
    engine = InferenceEngine(emissions, transitions)
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    
    # Identify common SNPs once
    common_pos = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
    common_set = set(common_pos)
    
    population_data = []

    # 3. Process the population
    for i, sample_id in enumerate(asw_samples):
        vcf = pysam.VariantFile(vcf_path)
        snp_positions, genotypes = [], []
        
        for record in vcf:
            if record.pos in common_set:
                gt = record.samples[sample_id].allele_indices
                snp_positions.append(record.pos)
                genotypes.append(gt)
        
        # Run inference
        results = engine.run_viterbi(snp_positions, genotypes, get_cm)
        labels = [
            "CEU" if state == "CEU_CEU" else "YRI" if state == "YRI_YRI" else "HET"
            for state in results
        ]
        
        # Calculate stats
        yri_count = labels.count("YRI")
        ceu_count = labels.count("CEU")
        het_count = labels.count("HET")
        total = len(labels)
        yri_pct = (yri_count / total) * 100
        ceu_pct = (ceu_count / total) * 100
        het_pct = (het_count / total) * 100
        
        population_data.append({"sample_id": sample_id, "yri_pct": yri_pct, "ceu_pct": ceu_pct, "het_pct": het_pct})
        
        if i % 10 == 0:
            print(f"    Progress: {i}/{len(asw_samples)} samples completed...")

    # 4. Save and Plot
    df = pd.DataFrame(population_data)
    df.to_csv(OUTPUT_DIR / "asw_ancestry_results.csv", index=False)
    
    plt.figure(figsize=(10, 6))
    plt.hist(df['yri_pct'], bins=15, color='skyblue', edgecolor='black')
    plt.title("Distribution of African (YRI) Ancestry in ASW Population (Chr 22 Slice)")
    plt.xlabel("Percentage of YRI Ancestry")
    plt.ylabel("Number of Individuals")
    plt.savefig(OUTPUT_DIR / "asw_population_histogram.png")
    print(f"\n[+] Done! Results saved to {OUTPUT_DIR / 'asw_ancestry_results.csv'} and {OUTPUT_DIR / 'asw_population_histogram.png'}")

if __name__ == "__main__":
    main()