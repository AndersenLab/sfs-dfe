#!/usr/bin/bash
base_path=`git rev-parse --show-toplevel`
echo ${base_path}

# Fetch Tajima's D
release='20170531'
curl http://storage.googleapis.com/elegansvariation.org/releases/${release}/tajima/WI.${release}.tajima.bed.gz > ${base_path}/data/tajima.bed.gz
curl http://storage.googleapis.com/elegansvariation.org/releases/${release}/tajima/WI.${release}.tajima.bed.gz.tbi > ${base_path}/data/tajima.bed.gz.tbi

# Get orthologs
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/orthologs/c_elegans.current.orthologs.txt.gz -O ${base_path}/data/orthologs.gz

# Pull out operons
zcat ${base_path}/data/WS256.annotations.gff3.gz | awk -v OFS="\t" '$0 ~ "operon" { print $1, $4, $5, "TRUE" }' | bgzip > ${base_path}/data/operons.bed.gz
tabix -p bed ${base_path}/data/operons.bed.gz

# Download expression data and analyze
wget ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/RNASeq_controls_FPKM/c_elegans.canonical_bioproject.current_development.RNASeq_controls_FPKM.dat -O ${base_path}/data/RNASeq_FPKM.txt
Rscript ${base_path}/scripts/process_RNAseq_data.R

# Merge in sp_34
bcftools reheader -s <(echo 'sp34') ${base_path}/data/sp34_results/vcf/merged.vcf.gz | \
bcftools view -O v | sed 's/\//|/g' | \
bcftools view -O z  > ${base_path}/data/sp34_results/vcf/sp34.vcf.gz
bcftools index ${base_path}/data/sp34_results/vcf/sp34.vcf.gz
bcftools merge -O z \
               -m all \
               ${base_path}/data/WI.20170531.impute.vcf.gz \
               ${base_path}/data/sp34_results/vcf/sp34.vcf.gz > ${base_path}/data/WI.20170531.sp34.impute.vcf.gz

# Perform snpeff annotation
cd ${base_path}/data/
vcfanno vcfanno.toml ${base_path}/data/WI.20170531.sp34.impute.vcf.gz | \
bcftools filter --set-GTs . --include 'GT != "0|1"' | \
# Used to determine location of promoter
snpeff eff -upDownStreamLen 2000 WS256 - | \
snpeff eff -upDownStreamLen 2000 -classic -geneId -oicr WS256 - | \
bcftools view -O z >  ${base_path}/data/WI.20170531.impute.snpeff.vcf.gz

