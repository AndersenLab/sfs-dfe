


vcf_in=Channel.fromPath("test.vcf.gz")


r2window = Channel.from([500, 1000, 5000, 10000, 50000])
r2max = Channel.from([0.05, 0.10, 0.15, 0.20, 0.25])

process ibdseq {
    input:
        file(vcf) from vcf_in

    """
    five_maf=$(bcftools query --list-samples ${vcf} | wc -l | awk '{ print \$0*0.05 }' | awk '{printf("%d\n", \$0+=\$0<0?0:0.9)}')
    #for chrom in I II III IV V X do
    for chrom in I do
    java -jar ibdseq.r1206.jar \
        gt=WI.20170531.impute.vcf.gz \
        out=haplotype_out \
        r2max=.8 \
        minalleles=\${five_maf} \
        nthreads=4 \
        chrom=\${chrom}
    done;
    """
}