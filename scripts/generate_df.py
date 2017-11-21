#!/usr/env/python

import sys
import math
import re
import bisect
from collections import Counter
from collections import defaultdict, OrderedDict
from grantham import grantham, three_letter_to_one
from subprocess import Popen, PIPE
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

# Load expression genes
expr = {}
for x in range(1,5):
    with open(repo_path() + f"/data/expression/q{x}.expression.tsv", 'r') as f:
        expr[x] = [x.strip() for x in f.readlines()]

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


def extract_aa(aa_change):
    aa1, aa2 = re.split('[0-9]+', aa_change[2:])
    try:
        return (three_letter_to_one[aa1], three_letter_to_one[aa2],)
    except:
        pass


def arm_or_center(chrom, pos):
    if chrom == 'MtDNA':
        return False
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
              'pseudogene',
              'splice_donor_variant']

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


biotypes = ['protein_coding',
            'pseudogene',
            'ncRNA',
            'miRNA',
            'piRNA',
            'tRNA',
            'lincRNA',
            'rRNA',
            'scRNA',
            'snoRNA',
            'snRNA',
            'asRNA',
            'Coding',
            'Noncoding']


TAJIMA_BINS = OrderedDict(
                [
                    ['tajima__lt_neg_2', -2],
                    ['tajima__neg_2_to_0', 0],
                    ['tajima__0_to_2', 2],
                    ['tajima__gt_2', 10]
                ]
                )

g_bins = 40
grantham_set = [str(int(math.floor(g/g_bins*1.0)*g_bins)) for g in range(0,max(grantham.values()),g_bins)]

TF = {0: "F", 1: "T", False: "F", True: "T"}

vcf = VCF("-")
samples = vcf.samples
outgroup_index = samples.index(outgroup)
sample_len = len(samples)

header_printed = False

for line in vcf:
    line_str = str(line)
    if not line_str.startswith("#") and line_str.count("./.") == 0:
        # Get Allele Counts
        ancestral_allele_count = sum(line.gt_types[outgroup_index] == line.gt_types) - 1 # Subtract one for outgroup
        if sum(line.gt_types == 0) <= sum(line.gt_types == 3):
            minor_allele_count = sum(line.gt_types == 0)
        else:
            minor_allele_count = sum(line.gt_types == 3)

        out = OrderedDict([
                ('chrom', line.CHROM),
                ('pos', line.POS),
                ('allele_frequency', line.aaf),
                ('ancestral_allele_count', ancestral_allele_count),
                ('minor_allele_count', minor_allele_count),
                ('N_GT', sample_len),
                ('outgroup', outgroup)
              ])

        #==============#
        # Degeneracy   #
        #==============#
        # match HGVS change and calculate degeneracy
        fold_set = list(map(str, range(0,5)))
        out.update(list(zip(["fold__" + x for x in fold_set], len(fold_set)*[False])))
        m = re.match(".*\|([A-Za-z]{3})/([A-Za-z]{3})\|.*", line_str)
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
            out['fold__' + site_type] = True

        #=====#
        # ANN #
        #=====# 

        ANN = line.INFO.get("ANN")
        if ANN:
            ANN_SET = [dict(zip(ANN_header, x.split("|"))) for x in ANN.split(",")]


        #========#
        # effect #
        #========#

        # Add a column for every effect type.
        out.update(list(zip(["effect__" + x for x in site_types], len(site_types)*[False])))
        for effect in set(sum([x['effect'].split("&") for x in ANN_SET], [])):
            out["effect__" + effect] = True

        #========#
        # impact #
        #========#
        impact_set = ['MODIFIER', 'LOW', 'MODERATE', 'HIGH']
        out.update(list(zip(["impact__" + x for x in impact_set], len(impact_set)*[False])))
        # impact
        out["is__gene"] = False
        if 'transcript' in str(line):
            out['is__gene'] = True
        for impact in {x['impact'] for x in ANN_SET}:
            out["impact__" + impact] = True

        #=========#
        # biotype #
        #=========#
        out.update(list(zip(["biotype__" + x for x in biotypes], len(biotypes)*[False])))
        for biotype in {x['transcript_biotype'] for x in ANN_SET}:
            if biotype:
                out["biotype__" + biotype] = True


        #=======#
        # dauer #
        #=======#
        out['dauer__gene'] = False
        for gene in {x['gene_id'] for x in ANN_SET if x['impact'] != 'MODIFIER'}:
            if gene in dauer_genes:
                out['dauer__gene'] = True
            else:
                out['dauer__gene'] = False

        #============#
        # Expression #
        #============#
        for q in range(1,5):
            out[f'expression__q{q}'] = False
            for gene in {x['gene_id'] for x in ANN_SET if x['impact'] != 'MODIFIER'}:
                if gene in expr[q]:
                    out[f'expression__q{q}'] = True

        #==============#
        # Arm v Center #
        #==============#
        # (defined by Rockman)
        aoc = arm_or_center(line.CHROM, line.POS)
        out['aoc__arm'] = False
        out['aoc__center'] = False
        if aoc:
            out['aoc__' + aoc] = True


        # operon
        out['operon__operon'] = line.INFO.get("operon") or False

        # Grantham score (non-synonymous only)
        for k in grantham_set:
            out.update({f"granthem__{k}": False})
        for x in {extract_aa(x['aa_change']) for x in ANN_SET if x['aa_change']}:
            if x in grantham.keys() and x[0] != x[1]:
                g = grantham[x]
                g_out = str(int(math.floor(g/40.0)*40))
                out['granthem__' + g_out] = True


        # Tajima's D
        for k in TAJIMA_BINS.keys():
            out.update({k: False})

        tajima = line.INFO.get('tajima')
        if tajima:
            bi = bisect.bisect_right(sorted(TAJIMA_BINS.values()), tajima)
            tbin = list(TAJIMA_BINS.keys())[bi]
            out[tbin] = True


        if not header_printed and line.CHROM == "I":
            print('\t'.join(out.keys()))
            header_printed = True
        for k, v in out.items():
            if '__' in k:
                out[k] = TF[v]
        print('\t'.join(out.keys()))
        print('\t'.join(map(str, out.values())))


