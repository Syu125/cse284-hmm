#!/bin/bash
# Prepare and run RFMix on the smaller chr22 slice (memory-efficient)

cd /mnt/c/Users/yufam/Workspace/cse284-hmm/data/slice_data

echo "[-] Preparing slice VCFs for RFMix..."

# Use sample lists created earlier (in parent dir)
# If they don't exist, create them now
if [ ! -f ../yri_ceu_samples.txt ]; then
    echo "[-] Creating sample lists..."
    awk '$2=="YRI" || $2=="CEU" {print $1}' ../integrated_call_samples_v3.20130502.ALL.panel | tail -n +2 > ../yri_ceu_samples.txt
    awk '$2=="ASW" {print $1}' ../integrated_call_samples_v3.20130502.ALL.panel > ../asw_samples.txt
fi

echo "[-] Creating reference VCF (YRI+CEU from slice)..."
bcftools view -S ../yri_ceu_samples.txt chr22_slice.vcf -Oz -o reference_yri_ceu_slice.vcf.gz
bcftools index -t reference_yri_ceu_slice.vcf.gz

echo "[-] Creating query VCF (ASW from slice)..."
bcftools view -S ../asw_samples.txt chr22_slice.vcf -Oz -o query_asw_slice.vcf.gz
bcftools index -t query_asw_slice.vcf.gz

echo "[-] Running RFMix on slice..."
rfmix \
  -f query_asw_slice.vcf.gz \
  -r reference_yri_ceu_slice.vcf.gz \
  -m ../rfmix_sample_map.txt \
  -g genetic_map_chr22_rfmix.txt \
  -o rfmix_slice_output \
  --chromosome=22 \
  -e 5

echo "[+] RFMix complete! Output in rfmix_slice_output.*"
