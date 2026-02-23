# Commands to download the data files for the project.

# Genetic Map
# wget https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr22.txt.gz
# gunzip genetic_map_GRCh37_chr22.txt.gz

# Sample Panel
# wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# # Small VCF Slice
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
REGION="22:16000000-17000000"

bcftools view "$URL" "$REGION" -O v -o chr22_slice.vcf