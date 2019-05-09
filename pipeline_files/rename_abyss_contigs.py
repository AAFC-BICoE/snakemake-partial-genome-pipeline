"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import glob
import argparse

def main():
    parser = argparse.ArgumentParser(description='Renames Abyss contigs to more closely match SPAdes')
    parser.add_argument("input", type=str,
                        help='Input File')
    parser.add_argument('output', type=str,
                        help='Output File')
    args = parser.parse_args()
    print("Renaming Contigs in {}".format(args.input))

    rename_contigs(args.input, args.output)

def rename_contigs(input, output):
    seqs = []
    with open(input, "r") as f:
        for seq in SeqIO.parse(f, 'fasta'):
            seq.name = ""
            split = seq.description.split(" ")
            header = "NODE_{}_length_{}_cov_{}".format(split[0],split[1],split[2])
            seq.id = header
            seq.description = ""
            seqs.append(seq)

    with open(output, "w") as g:
        SeqIO.write(seqs, handle=g, format="fasta")


if __name__ == "__main__":
    main()