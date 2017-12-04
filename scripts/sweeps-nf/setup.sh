bcftools view   WI.20170531.impute.vcf.gz | \
awk '$0 ~ "^#" { print } $0 !~ "^#" && NR % 500 == 0 { print }' | \
bcftools view -O z > test2.vcf.gz && bcftools index test2.vcf.gz