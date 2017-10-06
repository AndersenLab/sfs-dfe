#!/usr/env/bash

base_path=`git rev-parse --show-toplevel`
echo ${base_path}

zcat ~/.genome/WS245/WS245.fa.gz > ${base_path}/data/sp34/WS245.fa
faidx --split-files ${base_path}/data/sp34/WS245.fa

zcat ${base_path}/data/sp34/Caenorhabditis_sp34_NK74SC_v710.scaffolds.fa.gz > ${base_path}/data/sp34/Caenorhabditis_sp34_NK74SC_v710.scaffolds.fa
faidx --split-files ${base_path}/data/sp34/Caenorhabditis_sp34_NK74SC_v710.scaffolds.fa

lastz ${base_path}/data/sp34/chr1_sp34.fa ${base_path}/data/sp34/I.fa
lastz chr1_sp34.fa I.fa


lastz --step=500 --format=sam  MtDNA.fa[unmask] CSP34.Sp34_mitochondrion.fa | samtools mpileup --fasta-ref ~/.genome/WS245/WS245.fa.gz -g - | bcftools call  -O v -v -m - | head -n 50

lastz --step=50 --format=sam  MtDNA.fa[unmask] CSP34.Sp34_mitochondrion.fa | \
samtools mpileup --fasta-ref ~/.genome/WS245/WS245.fa.gz -g - | \
bcftools call  -O z -v -c - --pval-threshold 1.0 > Sp_34_MtDNA.vcf.gz

lastz --step=50 --format=sam --rdotplot=I.rdotplot I.fa[unmask] CSP34.Sp34_Chr1.fa | \
samtools mpileup --fasta-ref ~/.genome/WS245/WS245.fa.gz -g - | \
bcftools call  -O z -v -c - --pval-threshold 1.0 > Sp_34_I.vcf.gz

for i in I II III IV V X; do
	lastz --step=50 --format=sam  MtDNA.fa[unmask] CSP34.Sp34_${i}.fa | \
	samtools mpileup --fasta-ref ~/.genome/WS245/WS245.fa.gz -g - | \
	bcftools call  -O z -v -c - --pval-threshold 1.0 > Sp_34_${i}.vcf.gz
