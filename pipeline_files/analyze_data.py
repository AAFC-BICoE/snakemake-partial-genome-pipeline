from Bio import SeqIO

uces = {}
specimens = {}
with open ("/home/eyresj/Desktop/staphylinidae_analysis/all-taxa-merged.fasta") as f:
    for seq in SeqIO.parse(f, "fasta"):
        uce = str(seq.id.split("_")[0])
        specimen = str(seq.id.split("_")[1]).replace("-","_")
        if specimen in specimens:
            specimens[specimen].add(uce)
        else:
            specimens[specimen] = set(uce)
        # if specimen in uces:
        #     uces[specimen].append(uce)
        # else:
        #     uces[specimen] = [uce]

# for k, v1 in uces.items():
#     for k2, v2 in uces.items():
#         print(k, k2, len(set(v2) & set(v1)))
# for k, v in uces.items():
#     print(k, len(v))
import collections
# print(specimens)

a1 = "/home/eyresj/Desktop/staphylinidae_analysis/abyss_duplicates.txt"
r1 = "/home/eyresj/Desktop/staphylinidae_analysis/rnaspades_duplicates.txt"
s1 = "/home/eyresj/Desktop/staphylinidae_analysis/spades_duplicates.txt"



def count_mismatches(file_name):
    samples = {}
    with open(file_name) as f:
        lines = f.readlines()
        current_sample = lines[0].rstrip().replace("[","").split(" -")[0]
        for line in lines:
            if line.startswith("["):
                sample = line.rstrip().replace("[","").replace("_A ","").replace("_R ","").replace("_S ","").split("-")[0]
                if current_sample == sample:
                    if sample in samples:
                        pass
                    else:
                        samples[sample] = set()
                else:
                    current_sample = sample
                    samples[sample] = set()
            elif line.startswith("uce"):
                samples[current_sample].add(line.split(":")[0])
    return samples

abyss_samples = count_mismatches(a1)
rnaspades_samples = count_mismatches(r1)
spades_samples = count_mismatches(s1)



def merge_dicts(new_dict, merged_dict):
    for k, v in new_dict.items():
        if k in merged_dict:
            merged_dict[k] = merged_dict[k].union(new_dict[k])
        else:
            merged_dict[k] = new_dict[k]
    return merged_dict

merged_dict = {}
merged_dict = merge_dicts(abyss_samples,merged_dict)
merged_dict = merge_dicts(rnaspades_samples, merged_dict)
merged_dict = merge_dicts(spades_samples, merged_dict)

with open("/home/eyresj/Desktop/staphylinidae_analysis/compiled_results.csv", "w") as g:
    g.write("Sample, Targets Merged, Targets Filtered Merged, Unique Merged Filtered Targets, Unique Targets All\n")
    for k, v in merged_dict.items():
        target_sum = len(set(specimens[k]).union(set(v)))
        target_difference = len(set(v).difference(set(specimens[k])))
        print(k,len(v), target_difference)
        if target_sum > 800:
            print(k,target_sum)
        #print(k, len(specimens[k]), len(v))
        g.write("{},{},{},{},{}\n".format(k, len(specimens[k]), len(v), target_difference, target_sum))
