def clean_map_for_rfmix(input_path, output_path):
    last_gen = -1.0
    last_phys = -1
    
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        # RFMix v2 often requires this header
        outfile.write("chromosome\tpos\tgpos\n")
        
        for line in infile:
            parts = line.strip().split()
            if len(parts) < 4: continue
            
            try:
                # Standardize chromosome name to '22' (remove 'chr' if present)
                chrom = parts[0].replace('chr', '')
                phys = int(parts[1])
                gen = float(parts[3]) / 100.0 # Convert cM to Morgans
                
                # STRICTLY increasing check
                if phys > last_phys and gen > last_gen:
                    # Use TABS (\t) for separation
                    outfile.write(f"{chrom}\t{phys}\t{gen}\n")
                    last_phys = phys
                    last_gen = gen
            except ValueError:
                continue

clean_map_for_rfmix("genetic_map_GRCh37_chr22.txt", "genetic_map_chr22_rfmix.txt")