#!/usr/env/python

import sys
import math
import re
import csv
from collections import Counter
from collections import defaultdict


# Folded - How many strains have the rare allele
# 1-124 (cut in half?)

# Unfolded
# Number of strains that share allele with
# outgroup (QX1211)

class sfs:
    def __init__(self):
        self.folded = Counter()
        self.unfolded = Counter()


# Open codon table
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

sfs_out = defaultdict(sfs)

site_types = ['3_prime_UTR_variant',
              '5_prime_UTR_premature_start_codon_gain_variant',
              '5_prime_UTR_variant',
              'downstream_gene_variant',
              'initiator_codon_variant',
              'intergenic_region',
              'intron_variant',
              'missense_variant',
              'non_coding_transcript_exon_variant',
              'splice_region_variant',
              'start_lost',
              'stop_gained',
              'stop_lost',
              'stop_retained_variant',
              'synonymous_variant',
              'upstream_gene_variant',
              'pseudogene']


for line in sys.stdin:
    if line.startswith("#CHROM"):
        samples = line.strip().split("\t")[9:]
        outgroup = "QX1211"
        outgroup_index = samples.index(outgroup)
    if not line.startswith("#"):
        gts = line.strip().split("\t")[9:]
        gt_count = Counter(gts)
        mac = min(gt_count.values())
        oac = gts.count(gts[outgroup_index]) - 1
        
        # Whole sfs
        sfs_out['all_sites'].folded.update([mac])
        sfs_out['all_sites'].unfolded.update([oac])

        for site_type in site_types:
            if site_type in line:
                sfs_out[site_type].folded.update([mac])
                sfs_out[site_type].unfolded.update([oac])


        # match HGVS change and calculate degeneracy
        m = re.match(".*\|([A-Za-z]{3})/([A-Za-z]{3})\|.*", line)
        if m:
            codon = m.group(1)
            pos = [codon.index(x) for x in 'ATCG' if x in codon][0]
            changes = [(codon[:pos] + x + codon[pos+1:]).lower() for x in 'ATGC']
            changes.remove(codon.lower())
            orig_aa = codon_table[codon.lower()]
            changes = [codon_table[x] for x in changes if x != codon]
            degeneracy = changes.count(orig_aa) + 1
            aa = set(changes)
            print(codon, changes, aa, degeneracy)

for k, v in sfs_out.items():
    for sfs_type in ['folded', 'unfolded']:
        with open(f'results/{k}_{sfs_type}.sfs', 'w') as f:
            n = len(samples)
            chrom_count = {'folded': math.floor(n / 2),
                           'unfolded': n - 1}[sfs_type]
            f.write(f"{chrom_count} {sfs_type}\n")
            sfs = getattr(v, sfs_type)
            # Yes, we include freq 0.
            freq = map(str, [sfs[x] for x in range(0, chrom_count)])
            f.write(" ".join(freq) + "\n")
            mask = [1] + ([0] * (chrom_count - 1))
            f.write(' '.join(list(map(str, mask))))

            