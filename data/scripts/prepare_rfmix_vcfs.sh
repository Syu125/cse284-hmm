#!/bin/bash
# Split VCF into reference (YRI+CEU) and query (ASW) files for RFMix

PANEL="integrated_call_samples_v3.20130502.ALL.panel"
VCF="ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Extract YRI and CEU sample IDs for reference
awk '$2=="YRI" || $2=="CEU" {print $1}' $PANEL | tail -n +2 > yri_ceu_samples.txt

# Extract ASW sample IDs for query
awk '$2=="ASW" {print $1}' $PANEL > asw_samples.txt

echo "Reference samples (YRI+CEU): $(wc -l < yri_ceu_samples.txt)"
echo "Query samples (ASW): $(wc -l < asw_samples.txt)"

# Create reference VCF (only YRI and CEU)
bcftools view -S yri_ceu_samples.txt $VCF -Oz -o reference_yri_ceu_chr22.vcf.gz
bcftools index -t reference_yri_ceu_chr22.vcf.gz

# Create query VCF (only ASW)
bcftools view -S asw_samples.txt $VCF -Oz -o query_asw_chr22.vcf.gz
bcftools index -t query_asw_chr22.vcf.gz

echo "Created reference_yri_ceu_chr22.vcf.gz and query_asw_chr22.vcf.gz"
