"""
Utility functions for caching and optimization of the HMM ancestry inference pipeline.
"""

import pickle
import os
from data.data_parser import get_allele_frequencies


def get_cached_frequencies(vcf_path, population_samples, cache_path):
    """
    Load allele frequencies from cache file, or calculate and cache if not available.
    
    This is useful for large VCF files where frequency calculation is slow.
    
    Parameters:
    -----------
    vcf_path : str
        Path to the VCF file
    population_samples : list
        List of sample IDs for the population
    cache_path : str
        Path where the pickled frequencies will be stored
        
    Returns:
    --------
    dict : Dictionary mapping SNP positions to allele frequencies
    """
    # Load from cache if it exists
    if os.path.exists(cache_path):
        print(f"[+] Loading cached frequencies from {cache_path}...")
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    
    # Otherwise, calculate and save to cache
    print(f"[-] Cache not found. Calculating frequencies (this may take a while)...")
    freqs = get_allele_frequencies(vcf_path, population_samples)
    
    # Ensure cache directory exists
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    
    with open(cache_path, 'wb') as f:
        pickle.dump(freqs, f)
    
    print(f"[+] Cached frequencies to {cache_path}")
    return freqs
