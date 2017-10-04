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

# Perform snpeff annotation
cd ${base_path}/data/
vcfanno vcfanno.toml WI.20170531.impute.vcf.gz
bcftools view ${base_path}/data/WI.20170531.impute.vcf.gz | \
bcftools filter --set-GTs . --include 'GT != "0|1"' | \
# Used to determine location of promoter
snpeff eff -upDownStreamLen 2000 WS256 - | \
snpeff eff -upDownStreamLen 2000 -classic -geneId -oicr WS256 - | \
bcftools view -O z >  ${base_path}/data/WI.20170531.impute.snpeff.vcf.gz

