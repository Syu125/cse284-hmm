"""
This module contains functions for parsing genetic data.
"""

import pandas as pd
import pysam

def get_population_dict(panel_file):
    df = pd.read_csv(panel_file, sep='\t')
    populations = df.groupby('pop')['sample'].apply(list).to_dict()
    return populations

def get_allele_frequencies(vcf_path, sample_list):
    vcf = pysam.VariantFile(vcf_path)
    
    # Check for empty sample list
    if not sample_list:
        print("Warning: sample_list is empty!")
        return {}

    # Subset the VCF to our target population
    vcf.subset_samples(sample_list)
    
    freqs = {}
    first_record = True
    
    for record in vcf:
        # 1. Filter: Only handle SNPs with 1 REF and 1 ALT
        if len(record.alts) != 1:
            continue
        if len(record.ref) != 1 or len(record.alts[0]) != 1:
            continue

        # Debug check on the very first valid SNP we find
        if first_record:
            # Check how many samples pysam actually found in this VCF
            active_samples = list(record.samples.keys())
            print(f"Debug: Processing {len(active_samples)} samples in VCF.")
            first_record = False

        counts = [0, 0] # [Count of 0s, Count of 1s]
        
        for sample in record.samples.values():
            gts = sample.allele_indices
            if gts is not None:
                for allele in gts:
                    if allele is not None:
                        # allele will be 0 or 1
                        counts[allele] += 1
        
        total = sum(counts)
        if total > 0:
            freqs[record.pos] = counts[1] / total
            
    return freqs