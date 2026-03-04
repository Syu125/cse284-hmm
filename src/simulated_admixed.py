"""
Simulation test for HMM ancestry inference.

This script:
1. Creates a synthetic admixed individual by combining genotypes from real YRI and CEU samples
2. Uses the midpoint to split ancestry (roughly 50/50)
3. Runs Viterbi inference to recover the simulated ancestry pattern
4. Checks if the HMM detects the switch point around the midpoint

Use this to validate HMM performance on known ancestry boundaries.
"""

import pysam
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization import plot_ancestry


def main():
    # 1. Load Data
    vcf_path = "../data/chr22_slice.vcf.gz"
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "../data/genetic_map_GRCh37_chr22.txt"
    
    pops = get_population_dict(panel_path)
    phys, gen = get_genetic_map(map_path)
    
    # 2. Get Frequencies (Our Training Data)
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # 3. Create a "Hybrid" Person  
    # Grab one real YRI sample and one real CEU sample
    yri_sample_id = pops['YRI'][0]
    ceu_sample_id = pops['CEU'][0]
    
    # Load VCF records into memory for efficient lookup
    vcf = pysam.VariantFile(vcf_path)
    vcf_records = {record.pos: record for record in vcf if record.pos in yri_freqs and record.pos in ceu_freqs}
    
    snp_positions = []
    hybrid_genotypes = []
    
    # Combine common SNPs only
    common_pos = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
    midpoint = len(common_pos) // 2
    
    print(f"[-] Building hybrid from {yri_sample_id} (0-{midpoint}) and {ceu_sample_id} ({midpoint}+)")

    for i, pos in enumerate(common_pos):
        if pos not in vcf_records:
            continue
        record = vcf_records[pos]
        
        if i < midpoint:
            gt = record.samples[yri_sample_id].allele_indices
        else:
            gt = record.samples[ceu_sample_id].allele_indices
            
        snp_positions.append(pos)
        hybrid_genotypes.append(gt)

    # 4. Run the Engine
    emissions = EmissionModel(yri_freqs, ceu_freqs)
    transitions = TransitionModel(generations=10)
    engine = InferenceEngine(emissions, transitions)
    
    # Lambda for genetic distance lookup
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    
    print("[-] Running Viterbi Inference...")
    results = engine.run_viterbi(snp_positions, hybrid_genotypes, get_cm)
    
    # 5. Check the "Switch"
    print(f"\nResults (First 5): {results[:5]}")
    print(f"Results (Last 5):  {results[-5:]}")
    
    # Find where the HMM thinks the switch happened
    for i in range(1, len(results)):
        if results[i] != results[i-1]:
            print(f"\n[!] ANCESTRY SWITCH detected at index {i} (Position {snp_positions[i]})")
            print(f"Changed from {results[i-1]} to {results[i]}")
    
    plot_ancestry(snp_positions, results)

if __name__ == "__main__":
    main()