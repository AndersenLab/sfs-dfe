############################################################################### GENERATE ANNOTATION FILES
# put yourself in annotation directory
cd /projects/b1059/projects/Stefan/quick_sfs/

bcftools view --samples AB1,BRC20067,BRC20263,CB4852,CB4854,CB4856,CB4932,CX11262,CX11264,CX11271,CX11276,CX11285,CX11292,CX11307,CX11314,CX11315,DL200,DL226,DL238,ECA189,CB4853,CB4855,CB4857,CB4858,ECA252,PB306,ECA348,ECA349,ECA36,ECA369,ECA372,ECA396,ED3005,ED3011,ED3012,ED3017,ED3040,ED3046,ED3048,ED3049,ED3052,ED3073,ED3077,EG4347,EG4349,EG4724,EG4725,EG4946,GXW1,JT11398,JU1088,JU1172,JU1200,JU1212,JU1213,JU1242,JU1246,JU1249,JU1395,JU1400,JU1409,JU1440,JU1491,JU1530,JU1543,JU1568,JU1580,JU1581,JU1586,JU1652,JU1666,JU1792,JU1793,JU1808,JU1896,JU1934,JU2001,JU2007,JU2016,JU2017,JU2106,JU2131,JU2141,JU2234,JU2250,JU2257,JU2464,JU2466,JU2478,JU2513,JU2519,JU2522,JU2526,JU2534,JU2565,JU2566,JU2570,JU2572,JU2575,JU2576,JU2578,JU258,JU2581,JU2586,JU2587,JU2592,JU2593,JU2600,JU2610,JU2619,JU2800,JU2811,JU2825,JU2829,JU2838,JU2841,JU2853,JU2862,JU2866,JU2878,JU2879,JU2906,JU2907,JU310,JU311,JU3125,JU3127,JU3128,JU3132,JU3134,JU3135,JU3137,JU3140,JU3144,JU323,JU346,JU360,JU367,JU393,JU394,JU397,JU406,JU440,JU561,JU642,JU751,JU774,JU775,JU778,JU782,JU792,JU830,JU847,KR314,LKC34,LSJ1,MY1,MY10,MY16,MY18,MY2147,MY2212,MY23,MY2453,MY2530,MY2535,MY2573,MY2585,MY2693,MY2713,MY2741,MY518,MY679,MY772,MY795,MY920,N2,NIC1,NIC1049,NIC1107,NIC166,NIC195,NIC199,NIC2,NIC207,NIC231,NIC236,NIC242,NIC251,NIC252,NIC255,NIC256,NIC258,NIC259,NIC260,NIC261,NIC262,NIC265,NIC266,NIC267,NIC268,NIC269,NIC271,NIC272,NIC274,NIC275,NIC276,NIC277,NIC3,NIC501,NIC511,NIC513,NIC514,NIC515,NIC522,NIC523,NIC526,NIC527,NIC528,NIC529,PB303,PS2025,PX179,QG2075,QG556,QG557,QW947,QX1212,QX1233,QX1791,QX1792,QX1793,QX1794,RC301,WN2001,WN2033,WN2050,XZ1513,XZ1514,XZ1516 \
/projects/b105/projects/Stefan/CePopGen-nf/input_files/Ce330_GATK4_STRELKA2_Intersect.vcf.gz |\
bcftools file -i N_MISSING=0 -Oz -o WI.GATK_STRELKA.HTAphenotyped.vcf.gz

tabix -p vcf WI.GATK_STRELKA.HTAphenotyped.vcf.gz

mv WI* files/

# make ancestor bed
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t[%TGT]\n' files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > Annotation_Files/ANC.bed.gz

tabix Annotation_Files/ANC.bed.gz

# parse snpeff
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
awk -F"|" '$1=$1' OFS="\t" |\
awk '{print $1, $2, $3, $5}' OFS="\t" |\
awk '{if($4 == "") print $1, $2 = $2 - 1, $3, "intergenic_region";
	  else print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > Annotation_Files/SNPEFF_REGION.bed.gz

tabix Annotation_Files/SNPEFF_REGION.bed.gz

# make 4fold site annotation - 112045 sites
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
grep "synonymous_variant" |\
grep "protein_coding" |\
grep -v "splice_" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' |\
awk '{if(($4 ~ "p.Ser" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Val" || $4 ~ "p.Leu" || $4 ~ "p.Gly" || $4 ~ "p.Arg") && $7 == 0) print $1, $2=$2-1, $3, "4FOLD"}' OFS="\t" |\
bgzip > Annotation_Files/4FOLD_SITES.bed.gz

tabix Annotation_Files/4FOLD_SITES.bed.gz


# make 0fold site annotation - 138537 sites
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
grep "missense" |\
awk -F"||," '{print $1,$2,$3,$4}' OFS="\t" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' OFS="\t" |\
awk '$4 !~ "p.Ser" {print}' |\
awk '{if ($4 ~ "p.Met" || $4 ~ "p.Trp")  print $0, "0FOLD";
	  else if(($4 ~ "p.Phe" || $4 ~ "p.Ile" || $4 ~ "p.Val" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Tyr" || $4 ~ "p.His" || $4 ~ "p.Gln" || $4 ~ "p.Asn" || $4 ~ "p.Lys" || $4 ~ "p.Asp" || $4 ~ "p.Glu" || $4 ~ "p.Cys" || $4 ~ "p.Gly") && ($7 == 1 || $7 == 2)) print $0, "0FOLD";
	  else if(($4 ~ "p.Leu" || $4 ~ "p.Arg") && $7 == 1) print $0, "0FOLD"}' OFS="\t" |\
cut -f1,2,3,8 |\
awk '{print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > Annotation_Files/0FOLD_SITES.bed.gz

tabix Annotation_Files/0FOLD_SITES.bed.gz



# arms and centers
gunzip -c ARMS_CENTERS.bed.gz | awk '{print}' OFS="\t" | bgzip > ARMS_CENTERS.bed.gz
tabix ARMS_CENTERS.bed.gz

# extract intronic regions
# intron gff can be found in andersen lab dropbox worm reagents
gunzip -c Annotation_Files/intron.gff.gz |\
grep WormBase |\
awk '{print $1,$2=$4-1,$5,$3}' OFS="\t" |\
bgzip > Annotation_Files/WS245_INTRONS.bed.gz

tabix Annotation_Files/WS245_INTRONS.bed.gz

# extract exon sequences
gunzip -c Annotation_Files/exon.gff.gz |\
grep WormBase |\
awk '{print $1,$2=$4-1,$5,$3}' OFS="\t" |\
bgzip > Annotation_Files/WS245_EXONS.bed.gz

tabix Annotation_Files/WS245_EXONS.bed.gz

# indel variants
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%REF\t%ALT\n' files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
awk '{if(length($4) != length($5)) print $1, $2 = $2 - 1, $3, "INDEL"}' OFS="\t" |\
bgzip > Annotation_Files/INDEL.bed.gz

tabix Annotation_Files/INDEL.bed.gz


############################################################################### ANNOTATE VCF
# annotate vcf
vcfanno Annotation_Files/ANNOTATION_conf.toml files/WI.GATK_STRELKA.HTAphenotyped.vcf.gz |\
bcftools view -Oz -o files/WI.GATK_STRELKA.HTAphenotyped_ANNOTATED.vcf.gz

tabix -p vcf files/WI.GATK_STRELKA.HTAphenotyped_ANNOTATED.vcf.gz

# generate SFS input - remove MASKED regions from WS266
bcftools query --samples ^XZ1516 -f '[%CHROM\t%POS\t%REF\t%ALT\t%SNPEFF_REGION\t%AA\t%MASKED\t%4FOLD\t%0FOLD\t%GENOMIC_REGION\t%TGT\t%ANN\n]' files/WI.GATK_STRELKA.HTAphenotyped_ANNOTATED.vcf.gz |\
awk '$7 != "Masked" {print $0}' OFS="\t" |\
sed 's/\t\(.*\)\/\(.*\)\t/\t\1\t/g' |\
bgzip > files/SFS_INPUT.tsv.gz

# options
# 3_prime_UTR_variant
# 5_prime_UTR_premature_start_codon_gain_variant
# 5_prime_UTR_variant
# downstream_gene_variant
# initiator_codon_variant
# intergenic_region
# intron_variant
# missense_variant
# missense_variant&splice_region_variant
# non_coding_transcript_exon_variant
# splice_acceptor_variant&intron_variant
# splice_acceptor_variant&splice_donor_variant&intron_variant
# splice_acceptor_variant&splice_region_variant&intron_variant
# splice_donor_variant&intron_variant
# splice_region_variant
# splice_region_variant&intron_variant
# splice_region_variant&non_coding_transcript_exon_variant
# splice_region_variant&stop_retained_variant
# splice_region_variant&synonymous_variant
# start_lost
# start_lost&splice_region_variant
# stop_gained
# stop_gained&splice_region_variant
# stop_lost
# stop_lost&splice_region_variant
# stop_retained_variant
# synonymous_variant
# upstream_gene_variant

############################################################################### generate spectra

######## NEUTRAL SPECTRA - 4FOLD SITE - GENOME - 110531 sites
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t4FOLD\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_NEUTRAL_4FOLD_GENOME.tsv

# 4FOLD SITE - ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t4FOLD\t' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_NEUTRAL_4FOLD_ARMS.tsv

# 4FOLD SITE - CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t4FOLD\t' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_NEUTRAL_4FOLD_CENTERS.tsv



############################################################################### generate spectra

######## SELECTED SPECTRA

################################################################## 0FOLD
# 0 FOLD - GENOME
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t0FOLD\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_GENOME.tsv

# 0 FOLD - ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t0FOLD\t' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_ARMS.tsv

# 0 FOLD - CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep -P '\t0FOLD\t' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_CENTERS.tsv

################################################################## HIGH
# HIGH SITE - GENOME
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_HIGH_GENOME.tsv

# HIGH SITE - ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_HIGH_ARMS.tsv

# HIGH SITE - CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_HIGH_CENTERS.tsv

################################################################## HIGH - OFOLD
# GENOME
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_GENOME.tsv

# ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_ARMS.tsv

# CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_CENTERS.tsv


################################################################## MODERATE HIGH
# MODERATE HIGH SITE - GENOME
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|MODERATE|\||HIGH|' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_MODERATE_HIGH_GENOME.tsv

# MODERATE HIGH SITE - ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|MODERATE|\||HIGH|' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_MODERATE_HIGH_ARMS.tsv

# MODERATE HIGH SITE - CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|MODERATE|\||HIGH|' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_MODERATE_HIGH_CENTERS.tsv

################################################################## HIGH - OFOLD - MODERATE
# GENOME
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD\||MODERATE|' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_MODERATE_GENOME.tsv

# ARMS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD\||MODERATE|' |\
grep -P '\tARM\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_MODERATE_ARMS.tsv

# CENTERS
gunzip -c files/SFS_INPUT.tsv.gz |\
grep '|HIGH|\|0FOLD\||MODERATE|' |\
grep -P '\tCENTER\t' |\
datamash -s -g 1,2,6,11 count 11 |\
awk '{if($3 == $4) print $0, "ANCESTOR"; else print $0, "DERIVED"}' OFS="\t" > spectra_gatk/SFS_SELECTED_0FOLD_HIGH_MODERATE_CENTERS.tsv

