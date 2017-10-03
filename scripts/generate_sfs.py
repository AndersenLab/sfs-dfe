#!/usr/env/python

import sys
import math
import re
import csv
from collections import Counter
from collections import defaultdict
from itertools import chain

from os.path import dirname
from subprocess import Popen, PIPE
from pprint import pprint as pp


def repo_path():
    path = Popen(['git', 'rev-parse', '--show-toplevel'],
                 stdout=PIPE).communicate()[0]
    return str(path.strip(), encoding='UTF-8')

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

# chrom : left tip, left arm, center, right arm, right tip
# From Rockman Krugylak 2009
chrom_arm_center = {'I': (527, 3331, 7182, 3835, 197),
                    'II': (306, 4573, 7141, 2589, 670),
                    'III': (494, 3228, 6618, 2877, 567),
                    'IV': (720, 3176, 9074, 3742, 782),
                    'V': (643, 5254, 10653, 3787, 583),
                    'X': (572, 5565, 6343, 3937, 1302)}

for k, v in chrom_arm_center.items():
    chrom_arm_center[k] = (0,
                           sum(v[0:2]) * 1e3,
                           sum(v[0:3]) * 1e3)


def arm_or_center(chrom, pos):
    if chrom == 'MtDNA':
        return None
    ca = chrom_arm_center[chrom]
    if pos > ca[0]:
        c = 'arm'
    if pos > ca[1]:
        c = 'center'
    if pos > ca[2]:
        c = 'arm'
    return c


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

ANN_header = ["allele",
              "effect",
              "impact",
              "gene_name",
              "gene_id",
              "feature_type",
              "feature_id",
              "transcript_biotype",
              "exon_intron_rank",
              "nt_change",
              "aa_change",
              "cDNA_position/cDNA_len",
              "cds_position",
              "protein_position",
              "distance_to_feature",
              "error"]

for line in sys.stdin:
    if line.startswith("#CHROM"):
        samples = line.strip().split("\t")[9:]
        outgroup = "QX1211"
        outgroup_index = samples.index(outgroup)
    if not line.startswith("#"):
        sp_line = line.split("\t")
        gts = line.strip().split("\t")[9:]
        gt_count = Counter(gts)
        mac = min(gt_count.values())
        oac = gts.count(gts[outgroup_index]) - 1
        


        # Whole sfs
        sfs_out['all_sites'].folded.update([mac])
        sfs_out['all_sites'].unfolded.update([oac])

        INFO = [x.split("=")[1] for x in line.split("\t")[7].split(";") if x.split("=")[0] == "ANN"][0]
        ANN_SET = [dict(zip(ANN_header, x.split("|"))) for x in INFO.split(",")]


        for effect in set(sum([x['effect'].split("&") for x in ANN_SET], [])):
            if effect:
                sfs_out['effect_' + effect].folded.update([mac])
                sfs_out['effect_' + effect].unfolded.update([oac])

        # impact
        for impact in {x['impact'] for x in ANN_SET}:
            if impact:
                sfs_out['impact_' + impact].folded.update([mac])
                sfs_out['impact_' + impact].unfolded.update([oac])

        # biotype
        for biotype in {x['transcript_biotype'] for x in ANN_SET}:
            if biotype:
                sfs_out['biotype_' + biotype].folded.update([mac])
                sfs_out['biotype_' + biotype].unfolded.update([oac])

        # chromosome
        chrom = line.split("\t")[0]
        sfs_out['chrom_' + chrom].folded.update([mac])
        sfs_out['chrom_' + chrom].unfolded.update([oac])

        # Arms vs. Centers (defined by Rockman)
        aoc = arm_or_center(sp_line[0], int(sp_line[1]))
        if aoc in ['arm', 'center']:
            sfs_out['chrom_' + aoc].folded.update([mac])
            sfs_out['chrom_' + aoc].unfolded.update([oac])


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
            fold = {}

            # '1-fold' --> 0-fold
            if degeneracy == 1:
                degeneracy = 0

            site_type = f"{degeneracy}-fold"
            sfs_out['fold_' + site_type].folded.update([mac])
            sfs_out['fold_' + site_type].unfolded.update([oac])

for k, v in sfs_out.items():
    for sfs_type in ['folded', 'unfolded']:
        with open(f'{repo_path()}/results/{k}_{sfs_type}.sfs', 'w') as f:
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

            