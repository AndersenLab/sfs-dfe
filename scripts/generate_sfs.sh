#!/usr/bin/bash
set -e
base_path=`git rev-parse --show-toplevel`

mkdir -p ${base_path}/results/sp34
mkdir -p ${base_path}/results/QX1211
mkdir -p ${base_path}/results/XZ1516
mkdir -p ${base_path}/results/reports

# sp 34
bcftools view ${base_path}/data/vcf/WI.20170531.sp34.impute.snpeff.vcf.gz | \
python run_sfs.py sp34

# QX1211
bcftools view ${base_path}/data/vcf/WI.20170531.impute.snpeff.vcf.gz | \
python run_sfs.py QX1211

# XZ1516
bcftools view ${base_path}/data/vcf/WI.20170531.impute.snpeff.vcf.gz | \
python run_sfs.py XZ1516

# Generate reports
for outgroup in sp34; #QX1211 XZ1516;
do
    Rscript -e "rmarkdown::render('${base_path}/scripts/plot_sfs.Rmd', output_file = '${base_path}/results/reports/${outgroup}.html')" ${outgroup}
done;