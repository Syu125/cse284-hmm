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

import pysam
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization import plot_ancestry

def main():
    # 1. Setup
    vcf_path = "../data/chr22_slice.vcf"
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "../data/genetic_map_GRCh37_chr22.txt"
    
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
                # This line now gets a single string (sample_id)
                gt = record.samples[sample_id].allele_indices
                snp_positions.append(record.pos)
                genotypes.append(gt)

        # 4. Run Inference
        emissions = EmissionModel(yri_freqs, ceu_freqs)
        transitions = TransitionModel(generations=100)
        engine = InferenceEngine(emissions, transitions)
        
        get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
        results = engine.run_viterbi(snp_positions, genotypes, get_cm)
        
        # 5. Global Stats
        total = len(results)
        yri_pct = (results.count("YRI") / total) * 100
        print(f"    Summary: {yri_pct:.1f}% YRI | {100-yri_pct:.1f}% CEU")

        # 6. Visualize - ensure unique filename per sample!
        plot_ancestry(snp_positions, results, save_path=f"ancestry_{sample_id}.png")

if __name__ == "__main__":
    main()