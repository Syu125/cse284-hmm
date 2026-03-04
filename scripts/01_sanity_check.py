"""
Sanity check script for the HMM ancestry inference pipeline.

This script:
1. Loads a real individual from the 1000 Genomes dataset (NA19625)
2. Examines emission probabilities for the first 10 SNPs
3. Runs the Viterbi algorithm to infer ancestry
4. Validates that the inference is working as expected

Use this to test that your HMM implementation is functioning correctly
before running on larger datasets.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parent.parent / 'src'
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import pysam
import numpy as np
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization import plot_ancestry


def main():
    vcf_path = "data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    panel_path = "data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "data/genetic_map_GRCh37_chr22.txt"
    
    # Use an ASW (admixed) sample instead of pure YRI
    sample_id = None
    
    pops = get_population_dict(panel_path)
    
    # Try to find an ASW sample
    if 'ASW' in pops and pops['ASW']:
        sample_id = pops['ASW'][0]
        print(f"[-] Using ASW sample {sample_id} (admixed ancestry)")
    elif 'ACB' in pops and pops['ACB']:
        sample_id = pops['ACB'][0]
        print(f"[-] Using ACB sample {sample_id} (admixed ancestry)")
    else:
        # Fallback - but this won't have admixed ancestry
        sample_id = "NA19625"
        print(f"[!] ASW/ACB populations not found. Trying {sample_id} (may be pure YRI)")
    phys, gen = get_genetic_map(map_path)
    
    # Check genetic map coverage
    print(f"[-] Genetic map coverage:")
    print(f"    Physical positions: {phys[0]:.0f} bp to {phys[-1]:.0f} bp")
    print(f"    Genetic positions: {gen[0]:.4f} cM to {gen[-1]:.4f} cM")
    
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # Identify common SNPs that are WITHIN the genetic map range
    common_set = set(yri_freqs.keys()) & set(ceu_freqs.keys())
    map_min_pos = phys[0]
    map_max_pos = phys[-1]
    valid_snps = {pos for pos in common_set if map_min_pos <= pos <= map_max_pos}
    
    print(f"    Valid SNPs in both populations AND genetic map: {len(valid_snps)}")
    if not valid_snps:
        print("[!] ERROR: No SNPs found in valid range!")
        return
    
    # 2. Optimized Fetch - only valid SNPs
    vcf = pysam.VariantFile(vcf_path)
    snp_positions, genotypes = [], []
    
    print(f"[-] Fetching data for {sample_id}...")
    for record in vcf.fetch("22"):
        if record.pos in valid_snps:
            snp_positions.append(record.pos)
            genotypes.append(record.samples[sample_id].allele_indices)
    
    print(f"    Fetched {len(snp_positions)} SNPs for analysis")
    if len(snp_positions) < 100:
        print(f"[!] WARNING: Only {len(snp_positions)} SNPs - this may not be enough for reliable inference")
    
    # Check genotype composition
    gt_counts = {"(0, 0)": 0, "(0, 1)": 0, "(1, 1)": 0, "other": 0}
    for gt in genotypes:
        if gt == (0, 0):
            gt_counts["(0, 0)"] += 1
        elif gt == (0, 1) or gt == (1, 0):
            gt_counts["(0, 1)"] += 1
        elif gt == (1, 1):
            gt_counts["(1, 1)"] += 1
        else:
            gt_counts["other"] += 1
    
    total_gts = len(genotypes)
    print(f"\n[-] Genotype composition:")
    print(f"    (0,0): {gt_counts['(0, 0)']:6d} ({100*gt_counts['(0, 0)']/total_gts:5.1f}%)")
    print(f"    (0,1): {gt_counts['(0, 1)']:6d} ({100*gt_counts['(0, 1)']/total_gts:5.1f}%)")
    print(f"    (1,1): {gt_counts['(1, 1)']:6d} ({100*gt_counts['(1, 1)']/total_gts:5.1f}%)")
    
    if gt_counts['(0, 1)'] + gt_counts['(1, 1)'] < 100:
        print(f"\n[!] WARNING: Sample has very few alternate alleles - may be pure ancestry population")

    # 3. Check the "Emission" evidence for the first 10 SNPs
    emissions = EmissionModel(yri_freqs, ceu_freqs)
    
    print("\n[-] Checking evidence for the first 10 SNPs:")
    print(f"{'Position':<12} | {'GT':<6} | {'Log-P YRI':<10} | {'Log-P CEU':<10} | {'Winner'}")
    print("-" * 60)
    
    for i in range(10):
        pos = snp_positions[i]
        gt = genotypes[i]
        p_yri = emissions.get_log_emission("YRI", gt, pos)
        p_ceu = emissions.get_log_emission("CEU", gt, pos)
        winner = "CEU" if p_ceu > p_yri else "YRI"
        print(f"{pos:<12} | {str(gt):<6} | {p_yri:<10.2f} | {p_ceu:<10.2f} | {winner}")
    
    # Debug: Print emission prob for the first few SNPs
    # If YRI prob is way higher than CEU here, the HMM will never switch to CEU
    test_gt = genotypes[0]
    test_pos = snp_positions[0]
    p_yri = emissions.get_log_emission("YRI", test_gt, test_pos)
    p_ceu = emissions.get_log_emission("CEU", test_gt, test_pos)
    print(f"[-] SNP {test_pos} Log-Probs: YRI={p_yri:.2f}, CEU={p_ceu:.2f}")

    # 4. Debug: Check transition matrices at different distances
    print("\n[-] Transition model diagnostics:")
    transitions = TransitionModel(generations=10)
    
    # Test at small genetic distance (adjacent SNPs)
    trans_001 = transitions.get_transition_matrix(0, 0.01)  # 0.01 cM
    print(f"At 0.01 cM distance (very close SNPs):")
    print(f"  P(stay in same state) = {trans_001[0, 0]:.6f}")
    print(f"  P(switch states) = {trans_001[0, 1]:.6f}")
    
    # Test at larger distance
    trans_1 = transitions.get_transition_matrix(0, 1.0)  # 1.0 cM
    print(f"At 1.0 cM distance:")
    print(f"  P(stay in same state) = {trans_1[0, 0]:.6f}")
    print(f"  P(switch states) = {trans_1[0, 1]:.6f}")
    
    # Check actual genetic distances in our data
    print(f"\n[-] Genetic distance diagnostics:")
    genetic_dists = []
    for i in range(min(10, len(snp_positions)-1)):
        pos1, pos2 = snp_positions[i], snp_positions[i+1]
        cm1 = interpolate_genetic_position(pos1, phys, gen)
        cm2 = interpolate_genetic_position(pos2, phys, gen)
        dist_cm = abs(cm2 - cm1)
        genetic_dists.append(dist_cm)
        print(f"  SNP {i} to {i+1}: {dist_cm:.4f} cM (physical distance: {pos2-pos1} bp)")
    
    print(f"  Average cM distance between SNPs: {np.mean(genetic_dists):.4f}")
    
    # Check emission differences
    print(f"\n[-] Emission probability diagnostics:")
    print(f"{'Index':<6} | {'GT':<6} | {'Log YRI':<10} | {'Log CEU':<10} | {'Diff':<10}")
    print("-" * 50)
    for i in range(min(20, len(snp_positions))):
        pos = snp_positions[i]
        gt = genotypes[i]
        p_yri = emissions.get_log_emission("YRI", gt, pos)
        p_ceu = emissions.get_log_emission("CEU", gt, pos)
        diff = abs(p_yri - p_ceu)
        print(f"{i:<6} | {str(gt):<6} | {p_yri:<10.2f} | {p_ceu:<10.2f} | {diff:<10.2f}")
    
    # 4. Run Viterbi
    engine = InferenceEngine(emissions, transitions)
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    results = engine.run_viterbi(snp_positions, genotypes, get_cm)
    
    plot_ancestry(snp_positions, results, save_path=f"sanity_{sample_id}_full.png")
    print(f"[+] Plot saved to sanity_{sample_id}_full.png")

if __name__ == "__main__":
    main()