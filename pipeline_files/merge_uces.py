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
    parser.add_argument('-s', type=str,
                        help='SPAdes exploded-fastas folder', required=True)
    parser.add_argument('-r', type=str,
                        help='rnaSPAdes exploded-fastas folder', required=True)
    args = parser.parse_args()
    print("Merging SPAdes and rnaSPAdes UCEs together into {} directory".format(args.o))

    combine_uces(args.o, args.s, args.r)

def combine_uces(output_directory, spades_directory, rnaspades_directory):

    """
    Takes the UCES from an rnaSPAdes and SPAdes run and creates a seperate file taking only the best sequence per UCE
    :return:
    """

    # Verify folders exist
    if os.path.isdir(spades_directory) and os.path.isdir(rnaspades_directory):
        pass
    else:
        print("Missing either {} or {}".format(spades_directory, rnaspades_directory))
        return

    spades_fastas = glob.glob(os.path.join(spades_directory, "*.fasta"))
    rnaspades_fastas = glob.glob(os.path.join(rnaspades_directory, "*.fasta"))

    # Put all the contigs into a single dictionary
    final_uces = {}
    uce_change_log = []
    for fasta in spades_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-S.unaligned.fasta", "")
        final_uces[specimen_name] = {}
        for seq in SeqIO.parse(fasta, 'fasta'):
            uce = seq.description.split("|")[-1]
            final_uces[specimen_name][uce] = seq

    # Go through each rnaSPAdes UCE contig, compare with the SPAdes ones in dictionary and pick the best
    # The best is currently the longest
    for fasta in rnaspades_fastas:
        specimen = os.path.basename(fasta)
        specimen_name = specimen.replace("-R.unaligned.fasta", "")

        for seq in SeqIO.parse(fasta, 'fasta'):
            uce = seq.description.split("|")[-1]

            if uce in final_uces[specimen_name]:
                # Look for SPAdes uce
                spades_uce = final_uces[specimen_name][uce]
                if len(spades_uce.seq) > len(seq.seq):
                    # SPAdes UCE is longer so change nothing:
                    pass
                elif len(spades_uce.seq) < len(seq.seq):
                    # rnaSPAdes UCE is longer so use that one
                    final_uces[specimen_name][uce] = seq
                    change = "Specimen: {}, UCE: {}, Change: Replaced with rnaSPAdes Contig\n".format(specimen_name,
                                                                                                      uce)
                    uce_change_log.append(change)
                else:
                    # Lengths are equal between both UCEs, look for differences
                    if spades_uce.seq == seq.seq:
                        pass
                    else:
                        final_uces[specimen_name][uce] = seq
                        change = "Specimen: {}, UCE: {}, Change: Replaced with rnaSPAdes Contig\n".format(specimen_name,
                                                                                                          uce)
                        uce_change_log.append(change)
            else:
                # Add RNA UCE
                change = "Specimen: {}, UCE: {}, Change: Added rnaSPAdes Contig\n".format(specimen_name, uce)
                uce_change_log.append(change)
                final_uces[specimen_name][uce] = seq

    # Write Final UCES to merged file
    new_directory = os.path.join(output_directory, "merged_uces")
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)

    for key, value in final_uces.items():
        file_name = str(key) + "_merged.fasta"
        file_path = os.path.join(new_directory, file_name)
        with open(file_path, "w") as f:
            for k, seq in value.items():
                SeqIO.write(seq, handle=f, format="fasta")

    # Log all the changes made to the SPAdes UCE file to create the merged file
    file_name = "UCE_Change_Log.txt"
    file_path = os.path.join(new_directory, file_name)
    with open(file_path, "w") as f:
        f.writelines(uce_change_log)

if __name__ == "__main__":
    main()

