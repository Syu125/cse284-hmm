"""
Simulation test for HMM ancestry inference.

This script:
1. Creates a synthetic admixed individual by combining genotypes from real YRI and CEU samples
2. Uses the midpoint to split ancestry (roughly 50/50)
3. Runs Viterbi inference to recover the simulated ancestry pattern
4. Checks if the HMM detects the switch point around the midpoint

Use this to validate HMM performance on known ancestry boundaries.
"""

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent.parent / 'src'
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import pysam
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization import plot_ancestry


def main():
    # 1. Load Data
    vcf_path = "data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    panel_path = "data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "data/genetic_map_GRCh37_chr22.txt"
    
    pops = get_population_dict(panel_path)
    phys, gen = get_genetic_map(map_path)
    
    # Filter to valid genetic map range
    map_min_pos = phys[0]
    map_max_pos = phys[-1]
    print(f"[-] Genetic map coverage: {map_min_pos:.0f} to {map_max_pos:.0f} bp")
    
    # 2. Get Frequencies (Our Training Data)
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # Keep only SNPs in valid map range
    valid_snps = {pos for pos in set(yri_freqs.keys()) & set(ceu_freqs.keys()) 
                  if map_min_pos <= pos <= map_max_pos}
    print(f"[-] SNPs in both populations and genetic map: {len(valid_snps)}")
    
    # 3. Create a "Hybrid" Person  
    # Grab one real YRI sample and one real CEU sample
    yri_sample_id = pops['YRI'][0]
    ceu_sample_id = pops['CEU'][0]
    
    # Load VCF records into memory for efficient lookup
    vcf = pysam.VariantFile(vcf_path)
    vcf_records = {record.pos: record for record in vcf if record.pos in valid_snps}
    
    snp_positions = []
    hybrid_genotypes = []
    
    # Combine common SNPs only
    common_pos = sorted(valid_snps)
    midpoint = len(common_pos) // 2
    
    print(f"[-] Building hybrid from {yri_sample_id} (SNPs 0-{midpoint}) and {ceu_sample_id} (SNPs {midpoint}+)")

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
    transitions = TransitionModel(generations=100)  # Increase from 10 to 100 for larger blocks
    engine = InferenceEngine(emissions, transitions)
    
    # Lambda for genetic distance lookup
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    
    print("[-] Running Viterbi Inference...")
    results = engine.run_viterbi(snp_positions, hybrid_genotypes, get_cm)
    
    # 5. Check the "Switch"
    print(f"\n[-] Summary:")
    print(f"    Total SNPs: {len(results)}")
    print(f"    True ancestry boundary at SNP index: {midpoint}")
    print(f"    True boundary physical position: {snp_positions[midpoint] if midpoint < len(snp_positions) else 'N/A'}")
    
    # Find where the HMM thinks the switch happened
    print(f"\n[-] Inferred ancestry switches:")
    switch_count = 0
    inferred_switches = []
    for i in range(1, len(results)):
        if results[i] != results[i-1]:
            switch_count += 1
            inferred_switches.append(i)
            if switch_count <= 5:  # Show first 5 switches
                print(f"    Switch at index {i} (pos {snp_positions[i]}): {results[i-1]} → {results[i]}")
    
    if switch_count > 5:
        print(f"    ... and {switch_count - 5} more switches")
    
    # Check if the switch happened near the true boundary
    if inferred_switches:
        closest_switch = min(inferred_switches, key=lambda x: abs(x - midpoint))
        distance_to_boundary = abs(closest_switch - midpoint)
        print(f"\n[-] Switch detection accuracy:")
        print(f"    Closest inferred switch is {distance_to_boundary} SNPs away from true boundary")
        if distance_to_boundary < 100:
            print(f"    ✓ Good! Switch detected near expected location")
        elif distance_to_boundary < 500:
            print(f"    ~ Fair. Switch detected but somewhat offset")
        else:
            print(f"    ✗ Poor. Switch detected far from true location")
    
    plot_ancestry(snp_positions, results, true_switch_idx=midpoint)

if __name__ == "__main__":
    main()