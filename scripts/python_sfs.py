#!/usr/env/python

import sys
import math
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


sfs_out = defaultdict(sfs)


site_types = ['upstream_gene_variant',
              'missense',
              'intron_variant',
              'stop_gained',
              'intergenic_region']

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
            mask = [1] + ([0] * chrom_count)
            f.write(' '.join(list(map(str, mask))))

            