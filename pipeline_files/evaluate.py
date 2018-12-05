"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""

import os
import glob
import argparse
import csv


def main():
    parser = argparse.ArgumentParser(description='Combines various log files into a CSV')
    parser.add_argument('-s1', type=str,
                        help='SPAdes UCE Log Input', required=True)
    parser.add_argument('-r1', type=str,
                        help='rnaSPAdes UCE Log Input', required=True)
    parser.add_argument('-s2', type=str,
                        help='SPAdes UCE Output', required=True)
    parser.add_argument('-r2', type=str,
                        help='rnaSPAdes UCE Output', required=True)

    args = parser.parse_args()
    #print(args)
    summarize_uces(args.r1, args.r2)
    summarize_uces(args.s1, args.s2)

def summarize_uces(input_path, output_path):
    with open(output_path, "w") as g:
        with open (input_path) as f:

            index = 0
            index_start = 0
            index_end = 0
            lines = f.readlines()
            for line in lines:
                if "INFO - ---" in line:
                    if index_start > 0:
                        index_end = index
                    else:
                        index_start = index
                index += 1

            specimen_lines = lines[index_start+1: index_end]
            g.write("Species, UCEs, Contigs, Dupes, UCEs Filtered, Contigs Filtered\n")
            for line in specimen_lines:
                if "Writing" in line:
                    continue
                slice = line[76:]
                split = slice.split(" ")
                species = split[0].replace(":", "")
                uniques = split[1]
                contigs = split[5]
                dupes = split[7]
                removed = split[11]
                match = split[19]

                g.write("{},{},{},{},{},{}\n".format(species, uniques, contigs, dupes, removed, match))

if __name__ == "__main__":
    main()
