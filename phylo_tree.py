"""
This is the script for phylogenetic analysis

Input:

What it does:

Part 1:
"""
import os
import sys
import shutil

path = "fragments"
full_path = "08-fasta"
threads = 8


threads = sys.argv[1] # 8


def msa(path):
    """
    This function performs multiple sequence alignment
    It creates two directories:
    1. fragments_aligned: 8 fasta files for segment alignment results
    2. full_aligned: fasta file for concatenated sequence alignment result
    :param path: path to the 8 fragment fasta files
    :return: NA
    """
    # multiple sequence alignment of the fragment files
    os.system("mkdir fragments_aligned")

    for file in os.listdir(path):
        # only want to align fragments, not the reference genome that's
        # also under this directory
        id = file.split(".fasta")[0]
        if len(id) in [2, 3]:
            input = path + "/" + file
            os.system("mafft --anysymbol %s > fragments_aligned/%s_aligned.fasta" % (input, file))

    # multiple sequence alignment of the concatenated file
    os.system("mkdir full_aligned")
    os.system("mafft --anysymbol concatenated/concatenated.fasta > full_aligned/concatenated_aligned.fasta")


def raxml_concatenated():
    """
    This function runs RAxML to generate the concatenated tree.
    It creates a direcotry under full_aligned called raxml_output that contains
    the output files. raxml_output/bootstrapping contains the output after
    performing bootstrapping with RAxML

    :return: NA
    """
    os.chdir("full_aligned")
    # Run raxml on concatenated fasta file
    # We add extension (.nexus) to view the trees directly in dendroscope
    for file in os.listdir("."):
        os.system("raxmlHPC -m GTRCAT -p 12345 -s %s -n %s.nexus"%(file, file.split(".fasta")[0]))
    os.chdir("../")

    for file in os.listdir("full_aligned"):
        if file.endswith(".fasta"):
            ## simply the output file name
            ## only want the segment information!!! + aligned
            new_name = file.split(".fasta_aligned")[0]
            os.system("raxmlHPC -m GTRCAT -p 12345 -# 100 -s full_aligned/%s -n %sT1"%(file, new_name))
            os.system("raxmlHPC -m GTRCAT -p 12345 -b 12345 -# 100 -s full_aligned/%s -n %sT2"%(file, new_name))
            os.system("raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.%sT1 -z RAxML_bootstrap.%sT2 -n %sT3"%(new_name, new_name, new_name))

            # strict consensus tree
            os.system("raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.%sT2 -n %sT4"%(new_name, new_name))

    # Remove useless files
    for file in os.listdir("full_aligned"):
        if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
            os.remove("full_aligned/%s"%file)

    for file in os.listdir("."):
        if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
            os.remove("%s"%file)

    os.mkdir("full_aligned/raxml_output")
    os.mkdir("full_aligned/raxml_output/bootstrapping")

    # Move output files to corresponding directories
    for file in os.listdir("full_aligned"):
        if file.startswith("RAxML"):
            shutil.move("full_aligned/%s"%file, "full_aligned/raxml_output/bootstrapping")

    for file in os.listdir("."):
        if file.startswith("RAxML"):
            os.rename(file, "%s.nexus"%file)
            shutil.move("%s.nexus"%file, "full_aligned/raxml_output/bootstrapping")


def raxml_segments():
    """
    This function runs RAxML to generate segment trees.
    It creates a directory called raxml_output that contains the output files.
    raxml_output/bootstrapping contains the output after performing bootstrapping with RAxML

    :return: NA
    """
    os.chdir("fragments_aligned")
    for file in os.listdir("."):
        os.system("raxmlHPC -m GTRCAT -p 12345 -s %s -n %s.nexus"%(file, file.split(".fasta")[0]))
    os.chdir("../")

    os.mkdir("raxml_output")

    # move raxml output files to raxml_output
    for file in os.listdir("fragments_aligned"):
        if file.endswith(".nexus"):
            shutil.move("fragments_aligned/%s"%file, 'raxml_output')

    # remove useless files
    for file in os.listdir("fragments_aligned"):
        if file.endswith(".reduced"):
            os.remove("fragments_aligned/%s"%file)

    for file in os.listdir("raxml_output"):
        if file.endswith(".reduced") or file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree"):
            os.remove("raxml_output/%s"%file)


    # tree bootstrapping with raxml
    for file in os.listdir("fragments_aligned"):
        if file.endswith(".fasta"):
            # simply the output file name
            # only want the segment information!!! + aligned
            new_name = file.split(".fasta_aligned")[0]
            os.system("raxmlHPC -m GTRCAT -p 12345 -# 100 -s fragments_aligned/%s -n %sT1"%(file, new_name))
            os.system("raxmlHPC -m GTRCAT -p 12345 -b 12345 -# 100 -s fragments_aligned/%s -n %sT2"%(file, new_name))
            os.system("raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.%sT1 -z RAxML_bootstrap.%sT2 -n %sT3"%(new_name, new_name, new_name))

            # strict consensus tree
            os.system("raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.%sT2 -n %sT4"%(new_name, new_name))

    # remove useless files
    for file in os.listdir("fragments_aligned"):
        if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
            os.remove("fragments_aligned/%s"%file)

    for file in os.listdir("."):
        if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
            os.remove("%s"%file)

    # Move output files to corresponding directories
    # and add extension to view in dendroscope
    os.mkdir("raxml_output/bootstrapping")
    for file in os.listdir("."):
        if file.startswith("RAxML"):
            os.rename(file, "%s.nexus"%file)
            shutil.move("%s.nexus"%file, "raxml_output/bootstrapping")


def parsnp(reference):
    """
    This function runs Parsnp to generate segment trees and the concatenated tree.
    It creates a directory called parsnp_output that contains the output files.
    :param reference:
    :return: NA
    """
    # Parsnp on concatenated tree
    os.system("parsnp -v -x -c -r reference/B_ref.fasta -d UpperB -C 1000 -o parsnp_output/concatenated")

    #3. run Parsnp (fast2Tree) on 8 segments with reference
    # view tree (parsnp.ggr) with Gingr
    os.mkdir("parsnp_output")
    for dir in os.listdir("separate_frags"):
        if dir != "reference":
            os.mkdir("parsnp_output/%s"%dir)
            reference = "reference_%s.fasta"%dir
            os.system("parsnp -v -x -c -r separate_frags/reference/%s -d separate_frags/%s -C 1000 -o parsnp_output/%s"%(reference, dir, dir))



# Part 1
msa(path)

# Part 2
raxml_concatenated()
raxml_segments()

# Part 3
parsnp()

# Part 4
# To be added: BEAST 2
