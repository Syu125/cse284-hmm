import pysam
from data.data_parser import get_population_dict, get_allele_frequencies, get_genetic_map, interpolate_genetic_position
from hmm.emission import EmissionModel
from hmm.transition import TransitionModel
from hmm.viterbi import InferenceEngine
from visualization.karyogram import plot_ancestry

def main():
    # 1. Setup
    vcf_path = "../data/chr22_slice.vcf"
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    map_path = "../data/genetic_map_GRCh37_chr22.txt"
    
    pops = get_population_dict(panel_path)
    phys, gen = get_genetic_map(map_path)
    yri_freqs = get_allele_frequencies(vcf_path, pops['YRI'])
    ceu_freqs = get_allele_frequencies(vcf_path, pops['CEU'])
    
    # 2. Pick a "Real" admixed sample (ASW)
    asw_sample_id = pops['ASW'][0] # Take the first person in the ASW group
    print(f"[-] Analyzing real sample: {asw_sample_id}")
    
    # 3. Collect their data
    vcf = pysam.VariantFile(vcf_path)
    common_pos = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
    
    snp_positions = []
    genotypes = []
    
    # Note: Using the linear scan logic we discussed to avoid index errors
    common_set = set(common_pos)
    for record in vcf:
        if record.pos in common_set:
            gt = record.samples[asw_sample_id].allele_indices
            snp_positions.append(record.pos)
            genotypes.append(gt)

    # 4. Run Inference
    emissions = EmissionModel(yri_freqs, ceu_freqs)
    transitions = TransitionModel(generations=10) # ~10 generations is standard for ASW
    engine = InferenceEngine(emissions, transitions)
    
    get_cm = lambda x: interpolate_genetic_position(x, phys, gen)
    results = engine.run_viterbi(snp_positions, genotypes, get_cm)
    
    # 5. Visualize
    plot_ancestry(snp_positions, results, save_path=f"ancestry_{asw_sample_id}.png")

if __name__ == "__main__":
    main()