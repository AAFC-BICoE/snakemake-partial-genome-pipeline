"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
from Bio import SeqIO
import os
import glob
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description='Combines various log files into a CSV')
    parser.add_argument('-o', type=str,
                        help='Output File', required=True)
    parser.add_argument('-s', type=str,
                        help='SPAdes assembly metrics file', required=True)
    parser.add_argument('-r', type=str,
                        help='rnaSPAdes assembly metrics file', required=True)
    parser.add_argument('-su', type=str,
                        help='SPAdes uce metrics file', required=True)
    parser.add_argument('-ru', type=str,
                        help='rnaSPAdes uce metrics file', required=True)
    parser.add_argument('-sp', type=str,
                        help='SPAdes phylue log', required=True)
    parser.add_argument('-rp', type=str,
                        help='rnaSPAdes phyluce log', required=True)
    args = parser.parse_args()
    print("Merging SPAdes and rnaSPAdes UCEs together into {} directory".format(args.o))
    #print(args)
    extract_data(args)

def extract_data(args):
    rna_json = args.r
    spades_json = args.s
    with open(args.o) as g:
        g.write("")
        files = [args.s, args.r, args.su, args.ru]
        for file in files:
            with open(file) as f:
                the_dict =  csv.DictReader(f, delimiter="\t")
                for row in the_dict:
                    print(row["filename"], row["n_contigs"], row["ctg_N50"])


if __name__ == "__main__":
    main()
