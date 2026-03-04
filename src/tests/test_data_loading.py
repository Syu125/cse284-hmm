"""
Test script to validate data loading and basic HMM component functionality.

This script:
1. Checks that the panel file is being parsed correctly
2. Validates frequency calculation across populations
3. Tests emission probability calculations
4. Identifies SNPs with strong population differentiation
"""

from data.data_parser import get_population_dict, get_allele_frequencies
from hmm.emission import EmissionModel


def main():
    # 1. Validate panel parsing
    panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
    populations = get_population_dict(panel_path)
    
    print("[*] Population Summary:")
    print(f"    Populations found: {list(populations.keys())}")
    print(f"    YRI samples: {len(populations.get('YRI', []))}")
    print(f"    CEU samples: {len(populations.get('CEU', []))}")
    
    # 2. Test frequency calculation
    vcf_path = "../data/chr22_slice.vcf"
    print("\n[*] Calculating allele frequencies...")
    print("    YRI frequencies from VCF...")
    yri_freqs = get_allele_frequencies(vcf_path, populations['YRI'])
    
    print("    CEU frequencies from VCF...")
    ceu_freqs = get_allele_frequencies(vcf_path, populations['CEU'])
    
    # 3. Find high-signal SNPs
    common_snps = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
    
    if common_snps:
        print(f"\n[*] Found {len(common_snps)} common SNPs across both populations")
        print("    Searching for high differentiation SNPs (delta > 0.3)...")
        
        found_high_signal = False
        for test_pos in common_snps:
            f_yri = yri_freqs[test_pos]
            f_ceu = ceu_freqs[test_pos]
            delta = abs(f_yri - f_ceu)
            
            if delta > 0.3:
                print(f"\n    [!] HIGH-SIGNAL SNP at Position: {test_pos}")
                print(f"        YRI Frequency: {f_yri:.4f}")
                print(f"        CEU Frequency: {f_ceu:.4f}")
                print(f"        Delta: {delta:.4f}")
                found_high_signal = True
                break 
        
        if not found_high_signal:
            print("    [!] No high-signal SNPs found. Showing first 5 common SNPs:")
            for pos in common_snps[:5]:
                delta = abs(yri_freqs[pos] - ceu_freqs[pos])
                print(f"        Pos {pos}: YRI={yri_freqs[pos]:.4f}, CEU={ceu_freqs[pos]:.4f}, Delta={delta:.4f}")
    
    # 4. Test emission model
    print("\n[*] Testing EmissionModel...")
    model = EmissionModel(yri_freqs, ceu_freqs)
    
    # Test with a high-frequency allele in YRI
    test_snp = list(common_snps)[0] if common_snps else None
    if test_snp:
        test_genotype = (1, 1)  # Homozygous for alternate allele
        
        probs = model.get_emission_probs(test_snp, test_genotype)
        
        print(f"    Test SNP at position {test_snp}:")
        print(f"        Genotype: {test_genotype}")
        print(f"        Likelihood if YRI: {probs['YRI']:.6f}")
        print(f"        Likelihood if CEU: {probs['CEU']:.6f}")
        
        if probs['YRI'] > probs['CEU']:
            print(f"        -> More consistent with YRI ancestry")
        else:
            print(f"        -> More consistent with CEU ancestry")


if __name__ == "__main__":
    main()
