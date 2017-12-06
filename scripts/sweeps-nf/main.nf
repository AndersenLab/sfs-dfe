
params.vcf=""
ibdseq=file("ibdseq.r1206.jar")
process_ibd=file("process_ibd.R")
vcf_in=Channel.fromPath(params.vcf)

// Species the minimum number of samples carrying the minor allele. 
minalleles=Channel.from([0.0, 0.05])

// Specifies the number of markers in the sliding window used to detect correlated markers.
r2window = Channel.from([150, 750, 1500])

ibdtrim = Channel.from([0, 5, 10])

r2max = Channel.from([0.15, 0.5, 0.8])


ibd_set = vcf_in.combine(minalleles).combine(r2window).combine(ibdtrim).combine(r2max)


process ibdseq {

    echo true

    publishDir "results/", mode: 'copy'

    errorStrategy 'ignore'

    tag { "${minalleles} ${r2window} ${r2max} ${ibdtrim}" }

    echo true

    input:
        set file(vcf), val(minalleles), val(r2window), val(ibdtrim), val(r2max) from ibd_set

    output:
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.ibd.tsv")
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.haplen.png")
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.gw_sort.png")
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.gw.png")
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.haplotype.png")
        file("${minalleles}_${r2window}_${r2max}_${ibdtrim}.sweeps.tsv")

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
    cp haplotype.tsv ${minalleles}_${r2window}_${r2max}_${ibdtrim}.ibd.tsv
    cp max_haplotype_genome_wide.png ${minalleles}_${r2window}_${r2max}_${ibdtrim}.gw.png
    cp max_haplotype_sorted_genome_wide.png ${minalleles}_${r2window}_${r2max}_${ibdtrim}.gw_sort.png
    cp haplotype_length.png ${minalleles}_${r2window}_${r2max}_${ibdtrim}.haplen.png
    cp haplotype.png ${minalleles}_${r2window}_${r2max}_${ibdtrim}.haplotype.png
    cp sweep_summary.tsv ${minalleles}_${r2window}_${r2max}_${ibdtrim}.sweeps.tsv
    """
}

