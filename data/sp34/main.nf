
/*
	Setup Genomes
*/
results_dir = 'results'
target_file = file("data/sitelist.tsv.gz")

chrom_set = Channel.from([
						      ['I', 'AP018151.1', 'sp34'],
                              ['II', 'AP018152.1', 'sp34'],
                              ['III', 'AP018153.1', 'sp34'],
                              ['IV', 'AP018154.1', 'sp34'],
                              ['V', 'AP018155.1', 'sp34'],
                              ['X', 'AP018156.1', 'sp34'],
                              ['I', 'BX284601.5', 'ce'],
                              ['II', 'BX284602.5', 'ce'],
                              ['III', 'BX284603.4', 'ce'],
                              ['IV', 'BX284604.4', 'ce'],
                              ['V', 'BX284605.5', 'ce'],
                              ['X', 'BX284606.5', 'ce']
                         ])

process download_chromosome {
	
	executor 'local'

	tag { "${chrom} ${species}"}

	input:
		set val(chrom), val(accession), val(species) from chrom_set

	output:
		set val(species), val(chrom), file("${chrom}_${species}.fa") into chrom_out

	"""
		wget -qO- 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${accession}&rettype=fasta&retmode=fasta' | \
		tr '[:lower:]' '[:upper:]' | \
		awk 'NR == 1 { print "> ${chrom} [sp34]"; } NR > 1 { print }' > ${chrom}_${species}.fa
	"""
}

mito_chrom = Channel.from([
				['ce', 'MtDNA', file('mito_chrom/ce_MtDNA.fa')],
				['sp34', 'MtDNA', file('mito_chrom/sp34_MtDNA.fa')]
			])

ce_ch = Channel.create()
sp34_ch = Channel.create()

// Split apart by species
mito_chrom
	.concat(chrom_out)
	.choice(ce_ch, sp34_ch) { sp -> sp[0] == 'ce' ? 0 : 1 }

// Combine by chromosome
by_chrom_ch = ce_ch
                  .phase(sp34_ch) { it -> it[1] }
                  .map { it -> [it[0][1], it[0][2], it[1][2]] }



process lastz {

	publishDir results_dir + '/vcf', mode: 'copy', pattern: '*.vcf.gz'

	tag { chrom }

	input:
		set val(chrom), file(ce), file(sp34) from by_chrom_ch

	output:
		set val(chrom), file("${chrom}.rdotplot") into dotplots
		set file("${chrom}.vcf.gz") into vcf_by_chrom

	script:
		if (chrom == 'MtDNA') {
			range = ""
		} else {
			//range = "[34000..4000000]"
			range = ""
		}


	"""
		lastz --step=50 \\
			  --format=sam \\
			  --rdotplot=${chrom}.rdotplot \\
			  --gfextend \\
			  --chain \\
			  --gapped \\
			  ${ce}${range} \\
			  ${sp34} | \\
		samtools mpileup --fasta-ref ~/.genome/WS245/WS245.fa.gz -g - | \\
		bcftools call -O v -c --pval-threshold 1.0 \\
					  --targets-file ${target_file} - | \\
		awk '{ gsub("0/1", "1/1", \$0); print }' | \\
		bcftools filter --exclude 'sum(DP4[0] + DP4[1]) > 0 && sum(DP4[2] + DP4[3]) > 0' | \\
		bcftools view -O z > ${chrom}.vcf.gz
	"""

}

process concatenate_vcf {

	publishDir results_dir + '/vcf', mode: 'copy', pattern: '*.vcf.gz'

	input:
		file(chrom_set) from vcf_by_chrom.toSortedList()

	output:
		file("merged.vcf.gz")

	"""	
		chrom_list="`ls -v -1 | grep -v 'MtDNA'` MtDNA.vcf.gz"
		bcftools concat \${chrom_list} | \\
		bcftools view -O z > merged.vcf.gz
	"""
}


process plot_dotplot {

	publishDir results_dir + '/dotplot', mode: 'copy'

	errorStrategy 'ignore'

	tag { chrom }

	input:
		set val(chrom), file("${chrom}.rdotplot") from dotplots

	output:
		file("${chrom}.dotplot.png")

	"""
	#!/usr/bin/env Rscript --vanilla
	png("${chrom}.dotplot.png",
		width=1200,
		height=1200)
	dots = read.table("${chrom}.rdotplot",header=T)
    plot(dots/1E6, type="l", xlab="C. elegans", ylab="species 34")
    dev.off()
	"""

}