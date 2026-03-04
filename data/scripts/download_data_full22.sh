#!/bin/bash

# 1. Genetic Map (GRCh37)
echo "[-] Downloading Genetic Map..."
wget -nc https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr22.txt.gz
gunzip -f genetic_map_GRCh37_chr22.txt.gz

# 2. Sample Panel (1000 Genomes)
echo "[-] Downloading Sample Panel..."
wget -nc https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# 3. Full Chromosome 22 VCF
# Note: This is a large file (~500MB+). 
VCF_FILE="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$VCF_FILE"

echo "[-] Downloading Full Chromosome 22 VCF..."
wget -nc "$URL"

# 4. Download the existing Index (.tbi) 
# It is much faster to download the index than to recreate it.
echo "[-] Downloading VCF Index..."
wget -nc "$URL.tbi"

# 5. Verification
if [ -f "$VCF_FILE" ] && [ -f "$VCF_FILE.tbi" ]; then
    echo "[+] Success: Full Chromosome 22 and Index are ready."
else
    echo "[!] Error: Download failed. Checking for tabix to re-index locally..."
    # Fallback: if the .tbi download failed, try to index manually
    tabix -p vcf "$VCF_FILE"
fi