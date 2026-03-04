import pysam
import pandas as pd
import matplotlib.pyplot as plt
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
def main():
    # 1. Setup Data
    vcf_path = "../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "../data/genetic_map_GRCh37_chr22.txt"
    
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
        
        # Calculate stats
        yri_count = results.count("YRI")
        total = len(results)
        yri_pct = (yri_count / total) * 100
        
        population_data.append({"sample_id": sample_id, "yri_pct": yri_pct, "ceu_pct": 100 - yri_pct})
        
        if i % 10 == 0:
            print(f"    Progress: {i}/{len(asw_samples)} samples completed...")

    # 4. Save and Plot
    df = pd.DataFrame(population_data)
    df.to_csv("asw_ancestry_results.csv", index=False)
    
    plt.figure(figsize=(10, 6))
    plt.hist(df['yri_pct'], bins=15, color='skyblue', edgecolor='black')
    plt.title("Distribution of African (YRI) Ancestry in ASW Population (Chr 22 Slice)")
    plt.xlabel("Percentage of YRI Ancestry")
    plt.ylabel("Number of Individuals")
    plt.savefig("asw_population_histogram.png")
    print("\n[+] Done! Results saved to asw_ancestry_results.csv and asw_population_histogram.png")

if __name__ == "__main__":
    main()