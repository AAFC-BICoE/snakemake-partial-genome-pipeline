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
    parser = argparse.ArgumentParser(description='Merges Phyluce UCEs from SPAdes and rnaSPAdes')
    parser.add_argument('-o', type=str,
                        help='Output Folder', required=True)
    parser.add_argument('-i', type=str,
                        help='Input folder of merged fastas', required=True)
    args = parser.parse_args()
    print("Counts merged_uces into a summary file in {} directory".format(args.o))

    count_uces(args.o, args.i)


def count_uces(output_directory, input_directory):
    # Gather each specimen file produced from the Phyluce
    merged_fastas = glob.glob(os.path.join(input_directory, "*_merged.fasta"))

    # Put all the contigs into a single dictionary
    specimen_dict = {}
    for fasta in merged_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("_merged.fasta", "").replace("-","_")
        with open(fasta) as f:
            count = 0
            abyss_count = 0
            spades_count = 0
            rnaspades_count = 0
            for seq in SeqIO.parse(fasta, 'fasta'):
                if "_A" in seq.id[-2:]:
                    abyss_count += 1
                if "_R" in seq.id[-2:]:
                    rnaspades_count += 1
                if "_S" in seq.id[-2:]:
                    spades_count += 1
                count += 1
        if specimen_name in specimen_dict:
            specimen_dict[specimen_name] = [count, abyss_count, spades_count, rnaspades_count]
        else:
            specimen_dict[specimen_name] = [count, abyss_count, spades_count, rnaspades_count]

    output_file = os.path.join(output_directory, "merged_uce_summary.csv")
    with open(output_file, "w") as g:
        g.write("Specimen, Merged Targets, Abyss Contribution, SPAdes Contribution, rnaSPAdes Contribution\n")
        for key, value in specimen_dict.items():
            g.write("{},{},{},{},{}\n".format(key, value[0],value[1],value[2],value[3]))


if __name__ == "__main__":
    main()
