"""
This is the script for checking homologoous recombination
signal using 3SEQ

Input:
1. output_dir: the path to where the output 3 seq reports will be stored

NOTE:
    1. This script uses the fullrun mode
    2. This script runs 3seq on a 700 p value table
"""
import os
import sys
import shutil


output_dir = sys.argv[1]

for file in os.listdir("fragments_aligned"):
    # create a separate directory for each segment
    segment = file.split("."[0])
    os.mkdir(output_dir + "/%s"%segment)

    # first generate a p value table
    os.system("./3seq -g my3seqTable700 700")

    # run 3 seq
    os.system("./3seq -f fragments_aligned/%s"%file)

    # put files into output directory
    for file in os.listdir("."):
        if file.startswith("3s."):
            shutil.move(file, output_dir + "/%s"%segment)

# remove the p value table after running 3 seq
os.remove("my3seqTable700")