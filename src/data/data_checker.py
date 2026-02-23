from data_parser import get_population_dict, get_allele_frequencies

def main():
    
    # Check sample panel
    panel_path = "../../data/integrated_call_samples_v3.20130502.ALL.panel"
    populations = get_population_dict(panel_path)
    
    print(populations.keys())
    
    print(f"YRI samples: {len(populations.get('YRI', []))}")
    print(f"CEU samples: {len(populations.get('CEU', []))}")
    
    # Checking frequency calculator
    vcf_path = "../../data/chr22_slice.vcf"
    print("[-] Calculating YRI frequencies...")
    yri_freqs = get_allele_frequencies(vcf_path, populations['YRI'])
    
    print("[-] Calculating CEU frequencies...")
    ceu_freqs = get_allele_frequencies(vcf_path, populations['CEU'])
    
    # Find a SNP that exists in both
    common_snps = sorted(set(yri_freqs.keys()) & set(ceu_freqs.keys()))
    
    if common_snps:
        print(f"\nScanning {len(common_snps)} common SNPs for a signal...")
        
        # Let's find a SNP where the difference (Delta) is significant
        for test_pos in common_snps:
            f_yri = yri_freqs[test_pos]
            f_ceu = ceu_freqs[test_pos]
            delta = abs(f_yri - f_ceu)
            
            if delta > 0.3:  # Find a SNP with > 30% difference
                print(f"Found High-Signal SNP at Position: {test_pos}")
                print(f"YRI Frequency: {f_yri:.4f}")
                print(f"CEU Frequency: {f_ceu:.4f}")
                print(f"Delta: {delta:.4f}")
                break 
        else:
            print("No high-signal SNPs found in this small slice. Showing first 5 instead:")
            for pos in common_snps[:5]:
                print(f"Pos {pos}: YRI={yri_freqs[pos]} CEU={ceu_freqs[pos]}")
    
if __name__ == "__main__":
    main()