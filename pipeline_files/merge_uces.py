"""
Merges the UCE targets found in Abyss, SPAdes and rnaSPAdes assemblies into a single file suitable for further
processing with Phyluce phylogeny tools.
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
    parser.add_argument('-s', type=str,
                        help='SPAdes exploded-fastas folder', required=True)
    parser.add_argument('-r', type=str,
                        help='rnaSPAdes exploded-fastas folder', required=True)
    parser.add_argument('-a', type=str,
                        help='Abyss exploded-fastas folder', required=True)
    parser.add_argument('-u', type=str,
                        help='Abyss Unmerged exploded-fastas folder', required=True)
    args = parser.parse_args()
    print("Merging SPAdes and rnaSPAdes UCEs together into {} directory".format(args.o))

    combine_uces(args.o, args.s, args.r, args.a, args.u)


def combine_uces(output_directory, spades_directory, rnaspades_directory, abyss_directory, abyss_u_directory):
    """
    Takes the UCES from various assembly runs and creates a seperate file taking only the best sequence per UCE
    :return:
    """

    # Verify folders exist
    if os.path.isdir(spades_directory) and os.path.isdir(rnaspades_directory) and os.path.isdir(abyss_directory):
        pass
    else:
        print("Missing either {} or {} or {}".format(spades_directory, rnaspades_directory, abyss_directory))
        return

    # Gather each specimen file produced from the Phyluce
    spades_fastas = glob.glob(os.path.join(spades_directory, "*.fasta"))
    rnaspades_fastas = glob.glob(os.path.join(rnaspades_directory, "*.fasta"))
    abyss_fastas = glob.glob(os.path.join(abyss_directory, "*.fasta"))
    abyss_u_fastas = glob.glob(os.path.join(abyss_u_directory, "*.fasta"))
    # Put all the contigs into a single dictionary
    specimen_dict = {}
    for fasta in spades_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-S.unaligned.fasta", "")
        specimen_dict[specimen_name] = [fasta]

    for fasta in rnaspades_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-R.unaligned.fasta", "")
        if specimen_name in specimen_dict:
            specimen_dict[specimen_name].append(fasta)

    for fasta in abyss_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-A.unaligned.fasta", "")
        if specimen_name in specimen_dict:
            specimen_dict[specimen_name].append(fasta)

    for fasta in abyss_u_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-AU.unaligned.fasta", "")
        if specimen_name in specimen_dict:
            specimen_dict[specimen_name].append(fasta)

    # For each specimen, add all the UCES to a single dictionary from every file, then examine each UCE sequence and
    # choose the one with the greatest length. Write all filtered UCEs to both a merged file, and monolithic file
    for key, value in specimen_dict.items():
        all_uces = {}
        for fasta in value:
            for seq in SeqIO.parse(fasta, 'fasta'):
                uce = seq.description.split("|")[-1]
                if uce in all_uces:
                    all_uces[uce].append(seq)
                else:
                    all_uces[uce] = [seq]
        print(key, len(all_uces))

        final_uces = []
        for k, v in all_uces.items():
            max_uce = None
            max_length = 0
            for seq in v:
                if len(seq.seq) > max_length:
                    max_uce = seq
            final_uces.append(max_uce)

        # Write Final UCES to merged file
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        file_name = str(key) + "_merged.fasta"
        file_path = os.path.join(output_directory, file_name)
        with open(file_path, "w") as f:
            for seq in final_uces:
                SeqIO.write(seq, handle=f, format="fasta")

        file_name = "all-taxa-incomplete-merged-renamed.fasta"
        file_path = os.path.join(output_directory, file_name)
        with open(file_path, "a") as f:
            for seq in final_uces:
                uce = str(seq.id).split("_")[0]
                specimen = key
                seq.description = "|" + uce
                seq.id = uce + "_" + specimen
                SeqIO.write(seq, handle=f, format="fasta")

        # # Log all the changes made to the SPAdes UCE file to create the merged file
        # file_name = "UCE_Change_Log.txt"
        # file_path = os.path.join(new_directory, file_name)
        # with open(file_path, "a") as f:
        #     f.writelines(uce_change_log)


if __name__ == "__main__":
    main()

