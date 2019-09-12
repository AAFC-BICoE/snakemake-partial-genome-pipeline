"""
Creates a summary report using Phyluce logs and fastq metrics from BBTools
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""
import os
import argparse


def main():
    parser = argparse.ArgumentParser(description='Combines various log files into a CSV')
    parser.add_argument('-i', type=str,
                        help='UCE Log Input', required=True)
    parser.add_argument('-f', type=str,
                        help='Fastq Metrics from statswrapper.sh', required=True)
    parser.add_argument('-o', type=str,
                        help='UCE Output', required=True)

    args = parser.parse_args()
    summarize_uces(args.i, args.f, args.o)


def summarize_uces(input_path, fastq_metrics, output_path):
    with open(output_path, "w") as g:
        reads = {}

        with open(fastq_metrics) as f:
            lines = f.readlines()
            lines.pop(0)
            for line in lines:
                split = line.rstrip().split("\t")
                read_count = split[0]
                file_name = split[-1]
                sample_name = os.path.basename(file_name).\
                    replace("_L001_R1_001.fastq.gz", "").replace("_L001_R2_001.fastq.gz", "")
                reads[sample_name] = read_count

        with open(input_path) as f:

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
            g.write("Species, Reads, Targets, Contigs, Dupes, Targets Filtered, Contigs Filtered\n")
            for line in specimen_lines:
                if "Writing" in line:
                    continue
                sliced = line[76:]
                split = sliced.split(" ")
                species = split[0].replace(":", "")
                species_name = split[0].replace("_A:", "").replace("_S:", "").replace("_R:", "").replace("_AU:", "")
                read_count = 0
                if species_name in reads:
                    read_count = reads[species_name]
                uniques = split[1]
                contigs = split[5]
                dupes = split[7]
                removed = split[11]
                match = split[19]

                g.write("{},{},{},{},{},{},{}\n".format(species, read_count, uniques, contigs, dupes, removed, match))


if __name__ == "__main__":
    main()
