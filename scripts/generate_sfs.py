#!/usr/env/python

import sys
import math
import re
import csv
import bisect
import sys
from collections import Counter
from collections import defaultdict, OrderedDict
from itertools import chain
from grantham import grantham, three_letter_to_one
from os.path import dirname
from subprocess import Popen, PIPE
from pprint import pprint as pp
from cyvcf2 import VCF

outgroup = sys.argv[1]

def repo_path():
    """
        Retrieves the repo path
    """
    path = Popen(['git', 'rev-parse', '--show-toplevel'],
                 stdout=PIPE).communicate()[0]
    return str(path.strip(), encoding='UTF-8')

class sfs:

    ancestral_allele = None

    def __init__(self):
        self.folded = Counter()
        self.unfolded = Counter()

    def update_sfs(self, fold_freq, unfold_freq):
        self.folded.update(fold_freq)
        if sfs.ancestral_allele != './.':
            self.unfolded.update(unfold_freq)



# Open codon table
bases = ['t', 'c', 'a', 'g']
codons = [a+b+c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))

# Load dauer genes
with open(repo_path() + "/data/gene_set/dauer_genes.txt", 'r') as f:
    dauer_genes = [x.strip() for x in f.readlines()]

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

f = []

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

TAJIMA_BINS = OrderedDict(
                [
                    ['lt_neg_2', -2],
                    ['neg_2_to_0', 0],
                    ['0_to_2', 2],
                    ['gt_2', 10]
                ]
                )


vcf = VCF("-")
samples = vcf.samples
outgroup_index = samples.index(outgroup)
sample_len = len(samples)

for line in vcf:
    line_str = str(line)
    if not line_str.startswith("#") and line_str.count("./.") == 0:
        # Get Allele Counts
        ancestral_allele_count = sum(i.gt_types[outgroup_index] == i.gt_types) - 1 # Subtract one for outgroup
        if sum(i.gt_types == 0) <= sum(i.gt_types == 3):
            minor_allele_count = sum(i.gt_types == 0)
        else:
            minor_allele_count = sum(i.gt_types == 3)

        out = OrderedDict([
                ('chrom', line.CHROM),
                ('pos', line.POS),
                ('allele_frequency', line.aaf),
                ('ancestral_allele_count', ancestral_allele_count),
                ('minor_allele_count', minor_allele_count),
                ('N_GT', sample_len)
              ])

        # Add a column for every site type.
        out.update(list(zip(site_types, len(site_types)*[False])))

        #=====#
        # ANN #
        #=====#

        dict(i.INFO)['ANN'].split("|")
        ANN = i.INFO.get("ANN")
        if ANN:
            ANN_SET = [dict(zip(ANN_header, x.split("|"))) for x in ANN.split(",")]

        for effect in set(sum([x['effect'].split("&") for x in ANN_SET], [])):
            out[effect] = True


        # impact
        for impact in {x['impact'] for x in ANN_SET}:
            if impact:
                sfs_out['impact_' + impact].update_sfs([minor_allele_count], [ancestral_allele_count])

        # biotype
        for biotype in {x['transcript_biotype'] for x in ANN_SET}:
            if biotype:
                sfs_out['biotype_' + biotype].update_sfs([minor_allele_count], [ancestral_allele_count])

        # dauer genes
        for gene in {x['gene_id'] for x in ANN_SET if x['impact'] != 'MODIFIER'}:
            if gene in dauer_genes:
                sfs_out['dauer_yes'].update_sfs([minor_allele_count], [ancestral_allele_count])
            else:
                sfs_out['dauer_no']

        # chromosome
        chrom = line.split("\t")[0]
        sfs_out['chrom_' + chrom].update_sfs([minor_allele_count], [ancestral_allele_count])

        # Arms vs. Centers (defined by Rockman)
        aoc = arm_or_center(sp_line[0], int(sp_line[1]))
        if aoc in ['arm', 'center']:
            sfs_out['chrom_' + aoc].update_sfs([minor_allele_count], [ancestral_allele_count])

        def extract_aa(aa_change):
            aa1, aa2 = re.split('[0-9]+', aa_change[2:])
            try:
                return (three_letter_to_one[aa1], three_letter_to_one[aa2],)
            except:
                pass

        # Gene flags
        if 'operon' in sp_line[7]:
            sfs_out['gene_operon'].update_sfs([minor_allele_count], [ancestral_allele_count])

        # Grantham score (non-synonymous only)
        for x in {extract_aa(x['aa_change']) for x in ANN_SET if x['aa_change']}:
            if x in grantham.keys() and x[0] != x[1]:
                g = grantham[x]
                g_out = str(int(math.floor(g/40.0)*40))
                sfs_out['Grantham_' + g_out].update_sfs([minor_allele_count], [ancestral_allele_count])

        # Tajima's D
        m = re.match(".*tajima=([\-0-9\.]+);.*", line)
        if m:
            tajima = float(m.group(1))
            bi = bisect.bisect_right(sorted(TAJIMA_BINS.values()), tajima)
            tbin = list(TAJIMA_BINS.keys())[bi]
            sfs_out['tajima_' + tbin].update_sfs([minor_allele_count], [ancestral_allele_count])
        

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

            site_type = str(degeneracy)
            sfs_out['fold_' + site_type].folded.update([minor_allele_count])
            sfs_out['fold_' + site_type].unfolded.update([ancestral_allele_count])




for k, v in sfs_out.items():

    # Output format: fitdadi
    for sfs_type in ['folded', 'unfolded']:
        with open(f'{repo_path()}/results/{outgroup}/{k}_{sfs_type}.sfs', 'w') as f:
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

    # Output format: DFE-alpha
    for sfs_type in ['unfolded']:
        with open(f'{repo_path()}/results/fitdadi/{outgroup}/{k}_{sfs_type}.sfs', 'w') as f:
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
    with open(f'{repo_path()}/results/multidfe/{outgroup}/{k}_{sfs_type}.sfs', 'w') as f:
        n = len(samples)
        chrom_count = {'folded': math.floor(n / 2),
                       'unfolded': n - 1}[sfs_type]
        f.write(f"{n}\n")
        sfs = getattr(v, sfs_type)
        # Yes, we include freq 0.
        freq = map(str, [sfs[x] for x in range(0, chrom_count)])
        f.write(" ".join(freq) + "\n")
        mask = [1] + ([0] * (chrom_count - 1))
        f.write(' '.join(list(map(str, mask))))

            