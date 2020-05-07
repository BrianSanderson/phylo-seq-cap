#!/usr/bin/env python

"""
Summarize variation in nucleotide sequences between individuals within a
species as Nei's Pi.

    Parameters:
        vcf (str): Path to VCF containing variants. Can be compressed as .gz.
        output (str): Path to write output TSV
        num_vcf (int): Optional, number of variants in the VCF. If
            provided will generate more informative progress bar

    TODO:
        Currently hard coded for specific structure of VCF used for original
        analysis for reproducibility. Needs to be generalized to make more
        useful aside from reproducing my study.
"""

__author__ = "Brian J. Sanderson <brian@biologicallyrelevant.com>"


import vcfpy
import csv
import argparse
from os import path
from itertools import combinations
from tqdm import tqdm


def calculate_pi(vcf, output, individuals, num_vcf):
    """
    Quantify SNPs and indels for each position in the candidate probes
    and writes out table as TSV

    Parameters:
        vcf (str): Path to VCF containing variants. Can be compressed as .gz.
        output (str): Path to write output TSV
        individuals (dict): Dictionary mapping individuals to species names
        num_vcf (int): Optional, number of variants in the VCF. If
            provided will generate more informative progress bar

    Notes:
        This function writes out a tab-delimited text file with the following
        fields
            0: the gene
            1: the site coordinate
            2: the species name
            3: the site diversity (Pi)
    TODO:
        Currently hard coded for specific structure of VCF used for original
        analysis for reproducibility. Needs to be generalized to make more
        useful aside from reproducing my study.
    """
    Ind_keys = list(Individuals.keys())
    Ind_keys_ix = [i for i, _ in enumerate(Ind_keys)]

    with open(output, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(["Gene", "Position", "Species", "Pi"])
        with tqdm(total=num_vcf) as pbar:
            vcf_reader = vcfpy.Reader.from_path(vcf)
            for record in vcf_reader:
                # skipping Idesia because there is only one individual
                it = iter(Ind_keys_ix[1:])
                for x in it:
                    x1, x2 = x, next(it)
                    if all([record.calls[x1].gt_bases != (None, None),
                            record.calls[x2].gt_bases != (None, None)]):
                        if(all([any(x >= 2 for x in record.calls[x1].data['AD']),
                                any(x >= 2 for x in record.calls[x2].data['AD'])])):
                            first, second = record.calls[x1].gt_bases, record.calls[x2].gt_bases
                            both = first + second
                            comp = []
                            for a, b in combinations(both, 2):
                                comp.append(a != b)
                            if len(comp) != 0:
                                pi = sum(comp) / len(comp)
                                writer.writerow([record.CHROM, record.POS,
                                                 Individuals[Ind_keys[x]],
                                                 pi])
                            else:
                                writer.writerow([record.CHROM, record.POS,
                                                 Individuals[Ind_keys[x]],
                                                 0])
                        else:
                            writer.writerow([record.CHROM, record.POS,
                                             Individuals[Ind_keys[x]],
                                             "NA"])
                    elif any([record.calls[x1].gt_bases != (None, None),
                              record.calls[x2].gt_bases != (None, None)]):
                        if record.calls[x1].gt_bases != (None, None):
                            if(x >= 2 for x in record.calls[x1].data['AD']):
                                both = record.calls[x1].gt_bases
                                comp = []
                                for a, b in combinations(both, 2):
                                    comp.append(a != b)
                                if len(comp) != 0:
                                    pi = sum(comp) / len(comp)
                                    writer.writerow([record.CHROM, record.POS,
                                                     Individuals[Ind_keys[x]],
                                                     pi])
                                else:
                                    writer.writerow([record.CHROM, record.POS,
                                                     Individuals[Ind_keys[x]],
                                                     0])
                            elif record.calls[x2].gt_bases != (None, None):
                                if(x >= 2 for x in record.calls[x2].data['AD']):
                                    both = record.calls[x2].gt_bases
                                    comp = []
                                    for a, b in combinations(both, 2):
                                        comp.append(a != b)
                                    if len(comp) != 0:
                                        pi = sum(comp) / len(comp)
                                        writer.writerow([record.CHROM, record.POS,
                                                         Individuals[Ind_keys[x]],
                                                         pi])
                                    else:
                                        writer.writerow([record.CHROM, record.POS,
                                                         Individuals[Ind_keys[x]],
                                                         0])
                            else:
                                writer.writerow([record.CHROM, record.POS,
                                                 Individuals[Ind_keys[x]],
                                                 "NA"])
                    else:
                        writer.writerow([record.CHROM, record.POS,
                                         Individuals[Ind_keys[x]],
                                         "NA"])
                pbar.update(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument('-o', dest='output', metavar='example.tsv',
                               help='tab-delimited output of nucleotide diversity',
                               required=True)
    requiredNamed.add_argument('-v', dest='vcf', metavar='vcf',
                               help='the vcf with SNPs to quantify '
                               'degeneracy of each site', required=True)
    optionalNamed = parser.add_argument_group("optional named arguments")
    optionalNamed.add_argument("--num_vcf", dest="num_vcf", metavar="num_vcf",
                               help="number of sites in VCF; if provided will "
                               "generate an informative progress bar",
                               required=False, default=None)
    args = parser.parse_args()

    if path.exists(args.vcf):
        vcf_file = args.vcf
    else:
        raise FileNotFoundError

    if args.num_vcf is None:
        num_vcf = '-inf'
    else:
        num_vcf = int(args.num_vcf)

    Individuals = {"I_polycarpa_WGS-2": "I_polycarpa",
                   "P_balsamifera_mgr-01": "P_balsamifera",
                   "P_balsamifera_mgr-04": "P_balsamifera",
                   "P_mexicana_pm3": "P_mexicana",
                   "P_mexicana_pm5": "P_mexicana",
                   "P_tremula_r01-01": "P_tremula",
                   "P_tremula_r04-01": "P_tremula",
                   "S_exigua_se002": "S_exigua",
                   "S_exigua_se053": "S_exigua",
                   "S_nigra_sg037": "S_nigra",
                   "S_nigra_sg051": "S_nigra",
                   "S_phlebophylla_sp15m": "S_phlebophylla",
                   "S_phlebophylla_sp7f": "S_phlebophylla"}

    calculate_pi(vcf_file, args.output, Individuals, num_vcf)
