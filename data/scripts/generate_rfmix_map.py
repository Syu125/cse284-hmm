#!/usr/bin/env python3
"""Generate RFMix sample map from 1000 Genomes panel file."""

panel_path = "integrated_call_samples_v3.20130502.ALL.panel"
output_path = "rfmix_sample_map.txt"

# Create sample map: sample_id -> population (only YRI and CEU for reference)
with open(panel_path) as f, open(output_path, 'w') as out:
    next(f)  # skip header
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        sample_id, pop = parts[0], parts[1]
        
        # Only output reference populations
        if pop in ['YRI', 'CEU']:
            out.write(f"{sample_id}\t{pop}\n")

print(f"Created {output_path}")
