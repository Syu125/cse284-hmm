import pickle
import os

from data.data_parser import get_allele_frequencies, get_population_dict

def get_cached_frequencies(vcf_path, pop_indices, cache_path):
    # If the cache file exists, just load it
    if os.path.exists(cache_path):
        print(f"[-] Loading cached frequencies from {cache_path}...")
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    
    # Otherwise, calculate and save it
    print(f"[-] Cache not found. Calculating frequencies (this may take a while)...")
    freqs = get_allele_frequencies(vcf_path, pop_indices)
    
    with open(cache_path, 'wb') as f:
        pickle.dump(freqs, f)
    
    return freqs


vcf_path = "../data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
panel_path = "../data/integrated_call_samples_v3.20130502.ALL.panel"
map_path = "../data/genetic_map_GRCh37_chr22.txt"
pops = get_population_dict(panel_path)

# In your main() function:
yri_cache = "data/yri_freqs_chr22.pkl"
ceu_cache = "data/ceu_freqs_chr22.pkl"

yri_freqs = get_cached_frequencies(vcf_path, pops['YRI'], yri_cache)
ceu_freqs = get_cached_frequencies(vcf_path, pops['CEU'], ceu_cache)