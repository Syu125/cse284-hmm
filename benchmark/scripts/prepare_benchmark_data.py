#!/usr/bin/env python
"""
Prepare benchmark data by splitting chr22_slice.vcf.gz into:
- Query file: ASW samples
- Reference file: YRI + CEU samples
"""

import gzip
import subprocess
import sys
from pathlib import Path


def read_panel(panel_path):
    """Read panel file and return dict of {sample: population}"""
    samples_by_pop = {}
    
    with open(panel_path, 'r') as f:
        header = f.readline().strip().split('\t')
        sample_idx = header.index('sample')
        pop_idx = header.index('pop')
        
        for line in f:
            fields = line.strip().split('\t')
            sample = fields[sample_idx]
            pop = fields[pop_idx]
            samples_by_pop[sample] = pop
    
    return samples_by_pop


def extract_vcf_samples(vcf_path):
    """Extract sample names from VCF file"""
    samples = []
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                # Header line with samples
                fields = line.strip().split('\t')
                samples = fields[9:]  # Samples start at column 9
                break
    return samples


def filter_vcf(input_vcf, output_vcf, sample_list):
    """Use bcftools to filter VCF to specific samples"""
    samples_str = ','.join(sample_list)
    
    cmd = [
        'bcftools', 'view',
        '-s', samples_str,
        '-O', 'z',  # Output compressed VCF
        '-o', output_vcf,
        input_vcf
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}", file=sys.stderr)
        return False
    
    # Index the output file
    subprocess.run(['bcftools', 'index', '-t', output_vcf], check=True)
    return True


def main():
    # Paths
    workspace_root = Path(__file__).parent.parent.parent
    panel_path = workspace_root / 'data' / 'raw' / 'panels' / 'integrated_call_samples_v3.20130502.ALL.panel'
    input_vcf = workspace_root / 'data' / 'processed' / 'chr22_slice.vcf.gz'
    output_dir = Path(__file__).parent / 'data'
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    query_vcf = output_dir / 'query_asw_chr22.vcf.gz'
    reference_vcf = output_dir / 'reference_yri_ceu_chr22.vcf.gz'
    
    print("Step 1: Reading panel file...")
    samples_by_pop = read_panel(panel_path)
    
    print("Step 2: Extracting samples from VCF...")
    vcf_samples = extract_vcf_samples(input_vcf)
    
    # Categorize samples
    asw_samples = [s for s in vcf_samples if samples_by_pop.get(s) == 'ASW']
    yri_samples = [s for s in vcf_samples if samples_by_pop.get(s) == 'YRI']
    ceu_samples = [s for s in vcf_samples if samples_by_pop.get(s) == 'CEU']
    
    print(f"Found {len(asw_samples)} ASW samples")
    print(f"Found {len(yri_samples)} YRI samples")
    print(f"Found {len(ceu_samples)} CEU samples")
    
    if not asw_samples:
        print("Error: No ASW samples found! Check panel file.", file=sys.stderr)
        return 1
    
    if not (yri_samples or ceu_samples):
        print("Error: No YRI or CEU samples found! Check panel file.", file=sys.stderr)
        return 1
    
    # Extract query VCF (ASW)
    print(f"\nStep 3: Extracting query VCF ({len(asw_samples)} ASW samples)...")
    if not filter_vcf(str(input_vcf), str(query_vcf), asw_samples):
        return 1
    print(f"Created: {query_vcf}")
    
    # Extract reference VCF (YRI + CEU)
    print(f"\nStep 4: Extracting reference VCF ({len(yri_samples)} YRI + {len(ceu_samples)} CEU samples)...")
    reference_samples = yri_samples + ceu_samples
    if not filter_vcf(str(input_vcf), str(reference_vcf), reference_samples):
        return 1
    print(f"Created: {reference_vcf}")
    
    print("\nSuccess! Benchmark data files ready:")
    print(f"  Query:     {query_vcf}")
    print(f"  Reference: {reference_vcf}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
