#!/usr/bin/env python

"""
Characterize the degeneracy for each site in a FASTA file using the standard
codon table.

    Parameters:
        fasta (str): Path to FASTA file
        output (str): Path to write output TSV

    TODO:
        This is not a very efficient script. I am new at working with pandas,
        and I need to re-evaluate whether I can accomplish something similar
        without the loop.
"""

__author__ = "Brian J. Sanderson <brian@biologicallyrelevant.com>"

import csv
import pandas as pd
import argparse
from tqdm import tqdm
from os import path
from Bio import SeqIO


def site_degeneracy(targ_dict, codon_dict, output):
    """
    Iterate through FASTA file codon-by-codon, and output the codon base
    degeneracy for each position. Assumes that the FASTA sequence is in the
    correct reading frame.

    Parameters:
        targ_dict (dict): Dictionary of FASTA names (keys) and sequence
            (values)
        codon_dict (dict): Dictionary of codons (keys) and site degeneracy for
            each of the 3 positions (values)
        output (str): Path to write output TSV

    Notes:
        This function writes a tab-delimited text file with the following
        fields
            0: the gene
            1: the position
            2: the site degenracy (0, 2, 3, or 4-fold)
    """
    col_names = ['Gene', 'Position', 'Subs']
    with open(output, 'w', newline='\n') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerow(col_names)
        with tqdm(total=len(targ_dict)) as pbar:
            for key in targ_dict.keys():
                testSeq = targ_dict[key].seq
                i = 0
                while i < len(testSeq):
                    outLine = codon_dict[str(testSeq[i:i+3])]
                    newRows = pd.DataFrame([[key, i, outLine[0]],
                                            [key, i+1, outLine[1]],
                                            [key, i+2, outLine[2]]],
                                           columns=col_names)
                    newRows.to_csv(csvfile, mode='a', header=False,
                                   sep='\t', index=False)
                    i = i + 3
                pbar.update(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='prepare a table of site  '
                                     'degeneracies for all positions in an '
                                     'in-reading-frame FASTA file')
    requiredNamed = parser.add_argument_group("required named arguments")
    requiredNamed.add_argument('-f', dest='fasta', metavar='example.fasta',
                               help='in-reading-frame FASTA file',
                               required=True)
    requiredNamed.add_argument('-o', dest='output', metavar='example.tsv',
                               help='tab-delimited output of site degeracy',
                               required=True)
    args = parser.parse_args()

    if path.exists(args.fasta):
        fasta_file = args.fasta
    else:
        raise FileNotFoundError

    targ_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    codon_dict = {"TTT": [0, 0, 2], "TTC": [0, 0, 2], "TTA": [2, 0, 2], "TTG": [2, 0, 2],
                  "CTT": [0, 0, 4], "CTC": [0, 0, 4], "CTA": [2, 0, 4], "CTG": [2, 0, 4],
                  "ATT": [0, 0, 3], "ATC": [0, 0, 3], "ATA": [0, 0, 3], "ATG": [0, 0, 0],
                  "GTT": [0, 0, 4], "GTC": [0, 0, 4], "GTA": [0, 0, 4], "GTG": [0, 0, 4],
                  "TCT": [0, 0, 4], "TCC": [0, 0, 4], "TCA": [0, 0, 4], "TCG": [0, 0, 4],
                  "CCT": [0, 0, 4], "CCC": [0, 0, 4], "CCA": [0, 0, 4], "CCG": [0, 0, 4],
                  "ACT": [0, 0, 4], "ACC": [0, 0, 4], "ACA": [0, 0, 4], "ACG": [0, 0, 4],
                  "GCT": [0, 0, 4], "GCC": [0, 0, 4], "GCA": [0, 0, 4], "GCG": [0, 0, 4],
                  "TAT": [0, 0, 2], "TAC": [0, 0, 2], "TAA": [0, 0, 2], "TAG": [0, 0, 2],
                  "CAT": [0, 0, 2], "CAC": [0, 0, 2], "CAA": [0, 0, 2], "CAG": [0, 0, 2],
                  "AAT": [0, 0, 2], "AAC": [0, 0, 2], "AAA": [0, 0, 2], "AAG": [0, 0, 2],
                  "GAT": [0, 0, 2], "GAC": [0, 0, 2], "GAA": [0, 0, 2], "GAG": [0, 0, 2],
                  "TGT": [0, 0, 2], "TGC": [0, 0, 2], "TGA": [0, 0, 0], "TGG": [0, 0, 0],
                  "CGT": [0, 0, 4], "CGC": [0, 0, 4], "CGA": [2, 0, 4], "CGG": [2, 0, 4],
                  "AGT": [0, 0, 2], "AGC": [0, 0, 2], "AGA": [2, 0, 2], "AGG": [2, 0, 2],
                  "GGT": [0, 0, 4], "GGC": [0, 0, 4], "GGA": [0, 0, 4], "GGG": [0, 0, 4]}

    site_degeneracy(targ_dict, codon_dict, args.output)
