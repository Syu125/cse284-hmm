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
import numpy as np
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine


def main():
    vcf_path = "../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "../data/genetic_map_GRCh37_chr22.txt"
    sample_id = "NA19625"

    pops = get_population_dict(panel_path)
    phys, gen = get_genetic_map(map_path)
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # Identify common SNPs
    common_set = set(yri_freqs.keys()) & set(ceu_freqs.keys())
    
    # 2. Optimized Fetch
    vcf = pysam.VariantFile(vcf_path)
    snp_positions, genotypes = [], []
    
    print(f"[-] Fetching data for {sample_id}...")
    for record in vcf.fetch("22"):
        if record.pos in common_set:
            snp_positions.append(record.pos)
            genotypes.append(record.samples[sample_id].allele_indices)

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

    # 4. Run Viterbi
    transitions = TransitionModel(generations=100)
    engine = InferenceEngine(emissions, transitions)
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    results = engine.run_viterbi(snp_positions, genotypes, get_cm)
    
    plot_ancestry(snp_positions, results, save_path=f"sanity_{sample_id}_full.png")
    print(f"[+] Plot saved to sanity_{sample_id}_full.png")

if __name__ == "__main__":
    main()