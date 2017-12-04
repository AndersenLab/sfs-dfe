

ibdseq=file("ibdseq.r1206.jar")
process_ibd=file("process_ibd.R")
vcf_in=Channel.fromPath("test2.vcf.gz")

// Species the minimum number of samples carrying the minor allele. 
minalleles=Channel.from([0.00, 0.05])

// Specifies the number of markers in the sliding window used to detect correlated markers.
r2window = Channel.from([75, 150, 500, 1000])

ibdtrim = Channel.from([0, 3, 5, 10, 20])

r2max = Channel.from([0.8, 1.0])


ibd_set = vcf_in.combine(minalleles).combine(r2window).combine(ibdtrim).combine(r2max)


process ibdseq {

    tag { "${minalleles} ${r2window}" }

    echo true

    input:
        set file(vcf), val(minalleles), val(r2window), val(ibdtrim), val(r2max) from ibd_set

    output:
        set file("haplotype_tsv.ibd")

    """
    echo \$(bcftools query --list-samples ${vcf} | wc -l | awk '{ print \$0*${minalleles} }' | awk '{printf("%d\\n", \$0+=\$0<0?0:0.9)}')
    minalleles=\$(bcftools query --list-samples ${vcf} | wc -l | awk '{ print \$0*${minalleles} }' | awk '{printf("%d\\n", \$0+=\$0<0?0:0.9)}')
    if [[ \${minalleles} -lt 2 ]];
    then
        minalleles=2;
    fi;
    for chrom in I II III IV V X; do
        java -jar ${ibdseq} \\
            gt=${vcf} \\
            out=haplotype_\${chrom} \\
            ibdtrim=${ibdtrim} \\
            minalleles=\${minalleles} \\
            r2max=${r2max} \\
            nthreads=4 \\
            chrom=\${chrom}
        done;
    cat *.ibd | awk '{ print \$0 "\t${minalleles}\t${ibdtrim}\t${r2window}\t${r2max}" }' > haplotype.tsv
    Rscript ${process_ibd}
    """
}

