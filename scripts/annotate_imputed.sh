#!/usr/bin/bash
base_path=`git rev-parse --show-toplevel`
echo ${base_path}
bcftools view ${base_path}/data/WI.20170531.impute.vcf.gz | \
bcftools filter --set-GTs . --include 'GT != "0|1"' | \
# Used to determine location of promoter
snpeff eff -upDownStreamLen 2000 -geneId -oicr WS256 - | \
bcftools view -O z >  ${base_path}/data/WI.20170531.impute.snpeff.vcf.gz

