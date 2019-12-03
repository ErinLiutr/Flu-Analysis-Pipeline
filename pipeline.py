import os
import dendropy
import numpy
import seaborn
from biopython import Phylo
from matplotlib import pyplot as plt
from dendropy.calculate import treecompare


## Step 1 read mapping
## run read_mapping.py



## Step 2 classify .fasta files into 8 segments
# input_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/UpperA"
# output_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments"
# reference_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/reference"
#
# map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}
# segments = {}
#
# for file in os.listdir(input_path):
#     name = file.split(".fasta")[0]
#     lines = ""
#     with open(input_path + "/" + file, "r") as f:
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
# ref_file = reference_path + "/A_ref.fasta"
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
#     print new_entries[idx] + new_entries[idx+1]
#     segments[idx/2 + 1].append(new_entries[idx] + new_entries[idx+1])
#
#
# ## write into files
# for num in range(1, 9):
#     segment = map[num]
#     with open(output_path + "/%s.fasta"%segment, "wb") as f:
#         for lines in segments[num]:
#             f.write(lines + "\n")
#     f.close()
#

## classify fragments
# input_path2 = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments"
# output_path2 = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/separate_frags"
# os.mkdir(output_path2)
# for file in os.listdir(input_path2):
#     id = file.split(".fasta")[0]
#     if len(id) in [2,3]:
#         os.mkdir(output_path2 + "/" + id)
#         with open(input_path2 + "/" + file, "r") as f:
#             entries = f.readlines()
#         for idx in range(0, len(entries), 2):
#             filename = entries[idx].split(" ")[0].replace(">","") + ".fasta"
#             with open(output_path2 + "/" + id + "/%s"%filename, "wb") as f2:
#                 f2.write(entries[idx])
#                 f2.write(entries[idx+1])
#                 f2.close()
#         f.close()
#     else:
#         os.mkdir(output_path2 + "/" + "reference")
#         with open(input_path2 + "/" + file, "r") as f:
#             entries = f.readlines()
#         new_entries = []
#         lines = ""
#         for line in entries:
#             if line[0] == ">":
#                 if lines != "":
#                     new_entries.append(lines+"\n")
#                 new_entries.append(line)
#                 lines = ""
#             else:
#                 nl = line.replace("\n","")
#                 lines += nl
#         new_entries.append(lines)
#         for idx in range(0, len(new_entries), 2):
#             filename = "reference_" + map[idx/2 + 1] + ".fasta"
#             with open(output_path2 + "/" + "reference" + "/%s"%filename, "wb") as f2:
#                 f2.write(new_entries[idx])
#                 f2.write(new_entries[idx+1])
#                 f2.close()
#         f.close()

## Step 3 Phylogenetic analysis
## run analysis.py

## Visualize all the phylogenetic trees

def treeImage(self, newick, rooted=False, outgroup=False):
    """
        Given a newick string, creates an image of the tree.
        Used in L Statistic GUI.
    """
    plt.figure(figsize=(8, 4))
    plt.axis('off')
    ax = plt.subplot(1, 1, 1)
    ax.axis('off')

    # Create the tree object
    tree = Phylo.read(newick, "newick")
    tree.rooted = rooted

    if rooted:
        tree.root_with_outgroup(outgroup)

    # Create the tree image
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('test.png')


## Tree comparison
def treeComparison(type):
    """
    
    :param type:
    :return: 
    """
    groundTruthFile = ""
    estimationFile = ""

    tns = dendropy.TaxonNamespace()
    gtTree = dendropy.Tree.get(file=open(groundTruthFile, 'r'), schema='newick', taxon_namespace=tns)
    estimateTree = dendropy.Tree.get(file=open(estimationFile, 'r'), schema='newick', taxon_namespace=tns)

    # metrics, weighted RF is unsymmetric, unweighted RF is symmetric distance
    weightedRF = treecompare.weighted_robinson_foulds_distance(gtTree, estimateTree)
    unweightedRF = treecompare.unweighted_robinson_foulds_distance(gtTree, estimateTree)
    euclideanDist = treecompare.euclidean_distance(gtTree, estimateTree)

    print weightedRF
    print unweightedRF
    print euclideanDist

## Step 4 Transmission tree analysis

## Step 5 Transmission bottleneck analysis