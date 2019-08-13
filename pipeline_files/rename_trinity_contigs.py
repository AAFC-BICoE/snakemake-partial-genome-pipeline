"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import glob
import argparse
import re

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
        count = 0
        for seq in SeqIO.parse(f, 'fasta'):
            seq.name = ""
            m = re.search(r"len=(\d+)", str(seq.description))
            length = m.group(1)
            #split = seq.description.split(" ")
            #>TRINITY_DN28260_c0_g1_i1 len = 120 path = [0:0 - 119]
            header = "NODE_{}_length_{}_cov_{}".format(count,length,0)
            seq.id = header
            seq.description = ""
            seqs.append(seq)
            count += 1
    with open(output, "w") as g:
        SeqIO.write(seqs, handle=g, format="fasta")


if __name__ == "__main__":
    main()