#!/usr/bin/bash
# @author: Daniel E. Cook
# @date: 2017-10-10

set -e
base_path=`git rev-parse --show-toplevel`
echo ${base_path}


# DEFINE GLOBALS
GFF3_BUILD="WS260"
SNPEFF_BUILD="WS261"
CENDR_RELEASE='20170531'

cd ${base_path}/data/

# DEFINE MAIN FILES
IMPUTE_VCF=original_vcf/WI.${CENDR_RELEASE}.impute.vcf.gz


# Create directories
mkdir -p annotation
mkdir -p expression
mkdir -p gene_set
mkdir -p vcf
mkdir -p original_vcf

# Fetch imputed VCF
if [ ! -s ${IMPUTE_VCF} ];
then
    wget -O ${IMPUTE_VCF} https://storage.googleapis.com/elegansvariation.org/releases/${CENDR_RELEASE}/WI.${CENDR_RELEASE}.impute.vcf.gz
    wget -O ${IMPUTE_VCF}.csi https://storage.googleapis.com/elegansvariation.org/releases/${CENDR_RELEASE}/WI.${CENDR_RELEASE}.impute.vcf.gz.csi
fi

# Fetch Tajima's D
if [ ! -s annotation/tajima.bed.gz ];
then
    curl --fail https://storage.googleapis.com/elegansvariation.org/releases/${CENDR_RELEASE}/tajima/WI.${CENDR_RELEASE}.tajima.bed.gz > annotation/tajima.bed.gz
    curl --fail https://storage.googleapis.com/elegansvariation.org/releases/${CENDR_RELEASE}/tajima/WI.${CENDR_RELEASE}.tajima.bed.gz.tbi > annotation/tajima.bed.gz.tbi
fi

# Download gff3
if [ ! -s annotation/${GFF3_BUILD}.annotations.gff3.gz ];
then
    wget -O annotation/${GFF3_BUILD}.annotations.gff3.gz \
    ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/gff/c_elegans.PRJNA13758.${GFF3_BUILD}.annotations.gff3.gz
fi

# Get orthologs
if [ ! -s annotation/orthologs.gz ];
then
    wget -O annotation/orthologs.gz \
    ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/orthologs/c_elegans.current.orthologs.txt.gz 
fi


# Pull out operons
if [ ! -s annotation/operons.bed.gz.tbi ];
then
    zcat annotation/${GFF3_BUILD}.annotations.gff3.gz | \
    awk -v OFS="\t" '$0 ~ "operon" { print $1, $4, $5, "TRUE" }' | \
    bgzip > annotation/operons.bed.gz
    tabix -p bed annotation/operons.bed.gz
fi

# Download expression data and analyze
if [ ! -s expression/RNASeq_FPKM.txt.gz ];
then
    wget -O expression/RNASeq_FPKM.txt \
    ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/RNASeq_controls_FPKM/c_elegans.canonical_bioproject.current_development.RNASeq_controls_FPKM.dat
    gzip expression/RNASeq_FPKM.txt
    rm expression/RNASeq_FPKM.txt
    Rscript ${base_path}/scripts/process/expression.R
fi

# Dauer Data
if [ ! -s gene_set/dauer_genes.txt ];
then
    #Download dauer data and analyze (doi: 10.1242/dev.00363; http://dev.biologists.org/content/130/8/1621)
    wget http://dev.biologists.org/highwire/filestream/1201065/field_highwire_adjunct_files/6/TableS4.xls -O gene_set/dauer_genes.xls
    Rscript ${base_path}/scripts/process/dauer.R
fi

#===========#
# Sp 34 VCF #
#===========#

# First merge sp34 in...
if [ ! -s vcf/WI.${CENDR_RELEASE}.sp34.impute.vcf.gz.csi ];
then
    # First merge
    bcftools reheader -s <(echo 'sp34') ${base_path}/data/sp34_results/vcf/merged.vcf.gz | \
    bcftools view -O v | sed 's/\//|/g' | \
    bcftools view -O z  > sp34_results/vcf/sp34.vcf.gz
    bcftools index sp34_results/vcf/sp34.vcf.gz
    bcftools merge -O z \
                   -m all \
                   ${IMPUTE_VCF} \
                   sp34_results/vcf/sp34.vcf.gz > vcf/WI.${CENDR_RELEASE}.sp34.impute.vcf.gz
    bcftools index -f vcf/WI.${CENDR_RELEASE}.sp34.impute.vcf.gz

fi

# Annotate sp34 VCF
if [ ! -s vcf/WI.${CENDR_RELEASE}.sp34.impute.snpeff.vcf.gz.csi ]
then
    # Perform snpeff annotation for sp34 VCF
    vcfanno -p 8 vcfanno.toml vcf/WI.${CENDR_RELEASE}.sp34.impute.vcf.gz | \
    bcftools filter --set-GTs . --include 'GT != "0|1"' | \
    # Used to determine location of promoter
    snpeff eff -upDownStreamLen 2000 ${SNPEFF_BUILD} - | \
    snpeff eff -upDownStreamLen 2000 -classic -geneId -oicr ${SNPEFF_BUILD} - | \
    bcftools view -O z >  vcf/WI.${CENDR_RELEASE}.sp34.impute.snpeff.vcf.gz
    bcftools index -f vcf/WI.${CENDR_RELEASE}.sp34.impute.snpeff.vcf.gz
fi

#============#
# SnpEff VCF #
#============#

if [ ! -s vcf/WI.${CENDR_RELEASE}.impute.snpeff.vcf.gz.csi ];
then
    vcfanno -p 8 vcfanno.toml ${IMPUTE_VCF} | \
    bcftools filter --set-GTs . --include 'GT != "0|1"' | \
    # Used to determine location of promoter
    snpeff eff -upDownStreamLen 2000 ${SNPEFF_BUILD} - | \
    snpeff eff -upDownStreamLen 2000 -classic -geneId -oicr ${SNPEFF_BUILD} - | \
    bcftools view -O z > vcf/WI.${CENDR_RELEASE}.impute.snpeff.vcf.gz
    bcftools index -f vcf/WI.${CENDR_RELEASE}.impute.snpeff.vcf.gz
fi
