#!/bin/bash

# ============================================================================
# Data Download & Preparation Script for HMM Local Ancestry Inference
# ============================================================================
# 
# This script handles downloading and preparing all required data.
# 
# Usage:
#   bash data/prepare_data.sh slice    # Download ~50MB slice (quick)
#   bash data/prepare_data.sh full     # Download ~1GB full chromosome (production)
#   bash data/prepare_data.sh all      # Download everything
#

set -e  # Exit on error

DATASET_TYPE="${1:-slice}"

# Set up directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DIR="$SCRIPT_DIR/raw"
PROCESSED_DIR="$SCRIPT_DIR/processed"

mkdir -p "$RAW_DIR/vcf" "$RAW_DIR/panels" "$RAW_DIR/maps" "$PROCESSED_DIR"

echo "==============================================="
echo "HMM Data Preparation Script"
echo "==============================================="

# ============================================================================
# STAGE 1: Common downloads (needed for all datasets)
# ============================================================================

download_common_files() {
    echo ""
    echo "[*] Downloading common files (genetic map & sample panel)..."
    
    # Download genetic map
    if [ ! -f "$RAW_DIR/maps/genetic_map_GRCh37_chr22.txt" ]; then
        echo "    [-] Downloading genetic map..."
        cd "$RAW_DIR/maps"
        wget -q https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr22.txt.gz
        gunzip -f genetic_map_GRCh37_chr22.txt.gz
        echo "    [+] Genetic map ready"
    else
        echo "    [✓] Genetic map already exists"
    fi
    
    # Download sample panel
    if [ ! -f "$RAW_DIR/panels/integrated_call_samples_v3.20130502.ALL.panel" ]; then
        echo "    [-] Downloading sample panel..."
        cd "$RAW_DIR/panels"
        wget -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        echo "    [+] Sample panel ready"
    else
        echo "    [✓] Sample panel already exists"
    fi
}

# ============================================================================
# STAGE 2: Slice dataset (~50MB, recommended for testing)
# ============================================================================

download_slice() {
    echo ""
    echo "[*] Downloading chromosome 22 slice (50MB, ~5-10 min)..."
    
    if [ -f "$PROCESSED_DIR/chr22_slice.vcf" ]; then
        echo "    [✓] Slice already exists"
        return
    fi
    
    echo "    [-] Downloading VCF index..."
    cd "$RAW_DIR/vcf"
    wget -nc -q https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
    
    echo "    [-] Slicing chromosome 22 (16-17 Mb region)..."
    URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    REGION="22:16000000-17000000"
    
    bcftools view "$URL" "$REGION" -O v -o "$PROCESSED_DIR/chr22_slice.vcf" 2>/dev/null
    
    echo "    [-] Indexing VCF..."
    cd "$PROCESSED_DIR"
    bgzip -c chr22_slice.vcf > chr22_slice.vcf.gz
    tabix -p vcf chr22_slice.vcf.gz
    
    echo "    [+] Slice ready: $PROCESSED_DIR/chr22_slice.vcf"
}

# ============================================================================
# STAGE 3: Full dataset (~1GB, for production use)
# ============================================================================

download_full() {
    echo ""
    echo "[*] Downloading full chromosome 22 (1GB, ~20-30 min)..."
    
    if [ -f "$RAW_DIR/vcf/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ]; then
        echo "    [✓] Full dataset already exists"
        return
    fi
    
    VCF_FILE="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/$VCF_FILE"
    
    echo "    [-] Downloading full VCF (~1GB)..."
    cd "$RAW_DIR/vcf"
    wget -nc "$URL"
    
    echo "    [-] Downloading VCF index..."
    wget -nc "$URL.tbi"
    
    if [ -f "$VCF_FILE" ] && [ -f "$VCF_FILE.tbi" ]; then
        echo "    [+] Full chromosome 22 ready"
    else
        echo "    [!] Download incomplete, attempting to re-index..."
        tabix -p vcf "$VCF_FILE"
    fi
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

case "$DATASET_TYPE" in
    slice)
        download_common_files
        download_slice
        ;;
    full)
        download_common_files
        download_full
        ;;
    all)
        download_common_files
        download_slice
        download_full
        ;;
    *)
        echo ""
        echo "Usage: bash prepare_data.sh [TYPE]"
        echo ""
        echo "TYPE options:"
        echo "  slice   - Download small 50MB slice (quick, recommended for testing)"
        echo "  full    - Download full 1GB chromosome 22 (for production)"
        echo "  all     - Download both slice and full"
        echo ""
        echo "Examples:"
        echo "  bash data/prepare_data.sh slice"
        echo "  bash data/prepare_data.sh full"
        echo ""
        exit 1
        ;;
esac

echo ""
echo "==============================================="
echo "[+] Data preparation complete!"
echo "==============================================="
echo ""
echo "Data location:"
echo "  Raw files: $RAW_DIR/"
echo "  Processed: $PROCESSED_DIR/"
echo ""
echo "Next steps:"
echo "  1. Verify files exist"
echo "  2. Run: python scripts/01_sanity_check.py"
echo ""
