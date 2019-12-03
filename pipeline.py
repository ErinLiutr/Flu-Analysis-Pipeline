import os
import dendropy
import numpy
import seaborn
from matplotlib import pyplot as plt
from dendropy.calculate import treecompare


## Step 1 read mapping
## run read_mapping.py



## Step 2 classify .fasta files into 8 segments
input_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/UpperB/"
output_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments/"
concate_output_path =  "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/"
reference_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/reference/"

map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}
# segments = {}
#
# for file in os.listdir(input_path):
#     name = file.split(".fasta")[0]
#     lines = ""
#     with open(input_path + file, "r") as f:
#         for line in f.readlines():
#             lines += line.replace("\n","")
#
#         f.close()
#
#     entries = lines.split(">")
#     if "" in entries:
#         entries.remove("")
#
#     for seg in entries:
#         s = seg.split("segment ")
#         seg_num = int (s[1][0])
#         sequence = s[1][1:]
#         id = s[0]
#         newID = id.split()[-1]
#         id = " ".join(id.split()[:-1])
#         header = ">" + name + " " + newID + " " + id + "\n" + sequence
#         if seg_num in segments.keys():
#             segments[seg_num].append(header)
#         else:
#             segments[seg_num] = [header]
#
# ## add reference segments
# ref_file = reference_path + "/B_ref.fasta"
# with open(ref_file, "r") as f:
#     entries = f.readlines()
# f.close()
#
# new_entries = []
# lines = ""
# for line in entries:
#     if line[0] == ">":
#         if lines != "":
#             new_entries.append(lines+"\n")
#         header = ">" + "reference " + line.split(">")[1]
#         new_entries.append(header)
#         lines = ""
#     else:
#         nl = line.replace("\n","")
#         lines += nl
# new_entries.append(lines)
# #print new_entries
# for idx in range(0, len(new_entries), 2):
#     segments[idx/2 + 1].append(new_entries[idx] + new_entries[idx+1])
#
#
# os.mkdir(output_path)
# ## write into files
# for num in range(1, 9):
#     segment = map[num]
#     with open(output_path + "/%s.fasta"%segment, "wb") as f:
#         for lines in segments[num]:
#             f.write(lines + "\n")
#     f.close()


# classify fragments
input_path2 = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments"
output_path2 = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/separate_frags"
os.mkdir(output_path2)
for file in os.listdir(input_path2):
    id = file.split(".fasta")[0]
    if len(id) in [2,3]:
        os.mkdir(output_path2 + "/" + id)
        with open(input_path2 + "/" + file, "r") as f:
            entries = f.readlines()

        if "\n" in entries:
            entries.remove("\n")
        for idx in range(0, len(entries), 2):
            filename = entries[idx].split(" ")[0].replace(">","") + ".fasta"
            with open(output_path2 + "/" + id + "/%s"%filename, "wb") as f2:
                f2.write(entries[idx])
                f2.write(entries[idx+1])
                f2.close()
        f.close()
    else:
        os.mkdir(output_path2 + "/" + "reference")
        with open(input_path2 + "/" + file, "r") as f:
            entries = f.readlines()
        new_entries = []
        lines = ""
        for line in entries:
            if line[0] == ">":
                if lines != "":
                    new_entries.append(lines+"\n")
                new_entries.append(line)
                lines = ""
            else:
                nl = line.replace("\n","")
                lines += nl
        new_entries.append(lines)
        for idx in range(0, len(new_entries), 2):
            filename = "reference_" + map[idx/2 + 1] + ".fasta"
            with open(output_path2 + "/" + "reference" + "/%s"%filename, "wb") as f2:
                f2.write(new_entries[idx])
                f2.write(new_entries[idx+1])
                f2.close()
        f.close()

## genrate the concatenated.fasta file
# lines = []
#
# for file in os.listdir(input_path):
#     with open(input_path + file, "r") as f:
#         read_lines = f.readlines()
#         subject_info = file.split(".fasta")[0]
#         id = read_lines[0].split(" ")[-3]
#         lines.append("\n" + ">" + subject_info + " " + id + "\n")
#         for l in read_lines:
#             if not l[0] == ">":
#                 lines.append(l.replace("\n",""))
#     f.close()
#
# ## write the reference genome (Optional)
# lines.append("\n>reference B/District Of Columbia/03/2018\n")
# with open(reference_path + "B_ref.fasta") as f:
#     for line in f.readlines():
#         if not line[0] == ">":
#             lines.append(line.replace("\n",""))
#
# with open(concate_output_path + "concatenated.fasta", "wb") as f:
#     for idx in range(len(lines)):
#         if idx == 0:
#             f.write(lines[idx][1:])
#         else:
#             f.write(lines[idx])
#     f.close()

## Step 3 Phylogenetic analysis
## run analysis.py

## Visualize all the phylogenetic trees

# def treeImage(self, newick, rooted=False, outgroup=False):
#     """
#         Given a newick string, creates an image of the tree.
#         Used in L Statistic GUI.
#     """
#     plt.figure(figsize=(8, 4))
#     plt.axis('off')
#     ax = plt.subplot(1, 1, 1)
#     ax.axis('off')
#
#     # Create the tree object
#     tree = Phylo.read(newick, "newick")
#     tree.rooted = rooted
#
#     if rooted:
#         tree.root_with_outgroup(outgroup)
#
#     # Create the tree image
#     Phylo.draw(tree, axes=ax, do_show=False)
#     plt.savefig('test.png')

