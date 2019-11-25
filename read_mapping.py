# -*- coding: utf-8 -*-
import os
import sys
import pandas as pd

type = sys.argv[1]
# UpperA UpperB

reference = sys.argv[2]
# a_ref = "A_Maryland_90_2017.fasta"
# b_ref = "B_DOC_03_2018.fasta"

reads = sys.argv[3]
# readsA = "Upper_resp-fluA"
# readsB = "Upper_resp-fluB"

# metadata file 1
m1 = sys.argv[4]
# A = "fluA.csv"
# B = "fluB.csv"

m2 = sys.argv[5]
#metadataA = "metadataA.csv"
#metadataB = "metadataB.csv"

#
#
# print "---------------Read Mapping Starting for %s %s---------------"%(reads, type)
#
# # Setting up working directories
# dirs = ["references","01-trimmomatic", "02-bowtie2-sam", "03-bam", "04-sorted-bam", "05-merged-bam", "06-vcf", "07-coverage", "08-fasta"]
#
# for d in dirs:
#     if not os.path.exists(d):
#         os.mkdir(d)
#
# print "---------------Directories all set!---------------"
#
# # STEP 1 clean up reads with Trimmomatic
#
# # STEP 2 read mapping with bowtie2
# # build reference
# build_name = type[-1].lower() + "_build"
# os.system("bowtie2-build %s %s"%(reference, build_name))
# os.system("mv %s* references"%build_name)
#
# # for id in os.listdir(readsA):
# #     print "---------------Current id: %s---------------"%id
# #     for file in os.listdir(readsA + "/" + id):
# #         if "R1" in file:
# #             read1 = readsA + "/" + id + "/" + file
# #         if "R2" in file:
# #             read2 = readsA + "/" + id + "/" + file
# #
# #     output = "02-bowtie2-sam/UpperA/%s"%id
# #
# #     print "read 1: %s \nread 2: %s\noutput: %s\n"%(read1, read2, output)
# #     os.system("bowtie2 -x references/a_build --quiet -1 %s -2 %s -S %s.sam"%(read1, read2, output))
#
# for id in os.listdir(reads):
#     print "---------------Current id: %s---------------"%id
#     for file in os.listdir(reads + "/" + id):
#         if "R1" in file:
#             read1 = reads + "/" + id + "/" + file
#         if "R2" in file:
#             read2 = reads + "/" + id + "/" + file
#
#     output = "02-bowtie2-sam/%s/%s"%(type, id)
#
#     print "read 1: %s \nread 2: %s\noutput: %s\n"%(read1, read2, output)
#     os.system("bowtie2 -x references/%s --quiet -1 %s -2 %s -S %s.sam"%(build_name, read1, read2, output))
#
# print "---------------Read Mapping Finished---------------"
#
# STEP 3 convert to .bam & .bam.bai files with samtools
for file in os.listdir("02-bowtie2-sam/%s"%type):
     id = file.split(".")[0]
     print id
     os.system("samtools view -Sb 02-bowtie2-sam/%s/%s > 03-bam/%s/%s.bam"%(type, file, type, id))

     os.system("samtools sort 03-bam/%s/%s.bam -o 04-sorted-bam/%s/%s.sorted.bam"%(type, id, type, id))

print "---------------Conversion to Bam Files Finished---------------"

# STEP 4 merge the reads from same subject

df = pd.read_csv(m1)

# format: keys: sample name; values: a list of two tuples, each represents one replicate
# inside each tuple, the first one comes from Q the second one from F
samples = {}
for idx, row in df.iterrows():
    name = row["Sample Name"].split("-")[0]
    replicate1 = (row["VDBPM ID"], row["VDBPM ID.2"])
    replicate2 = (row["VDBPM ID.1"], row["VDBPM ID.3"])
    samples[name] = [replicate1, replicate2]

# merge the the Q and F for each replicate

for name, ids in samples.items():
    for id in ids:
        merged = str(id[0]) + "-" + str(id[1])
        os.system("samtools merge 05-merged-bam/%s/%s-merged.bam 04-sorted-bam/%s/%s.sorted.bam 04-sorted-bam/%s/%s.sorted.bam"
                  %(type, merged, type, id[0], type, id[1]))

for file in os.listdir("05-merged-bam/%s"%type):
    id = file.split(".bam")[0]
    print "samtools index 05-merged-bam/%s/%s.bam"%(type, id)
    os.system("s 05-merged-bam/%s/%s.bam"%(type, id))

print "---------------Merge Bam Files finished---------------"

# STEP 5 Variant calling
# Build fasta index
os.system("lofreq faidx %s"%reference)

for file in os.listdir("05-merged-bam/%s"%type):
    if file.endswith("bam"):
        id = file.split(".")[0]
        os.system("lofreq call-parallel --pp-threads 2 -f %s -o 06-vcf/%s/%s.vcf 05-merged-bam/%s/%s"%(reference, type, id, type, file))

print "---------------Variant Calling Finished---------------"

# STEP 6 Get the coverage information (write into .coverage files)
for file in os.listdir("05-merged-bam/%s"%type):
    if file.endswith("bam"):
        print file
        id = file.split(".")[0]
        os.system("samtools mpileup 05-merged-bam/%s/%s -o 07-coverage/%s/%s.coverage"%(type, file, type, id))

print "---------------Generating .coverage Files Finished---------------"

# STEP 7 Generate assembled .fasta files
for file in os.listdir("05-merged-bam/%s"%type):
    if file.endswith("bam"):
        id = file.split("-merged")[0]
        print id

        os.system("python fasta.py %s %s %s %s"%(id, m2, type, reference))

print "---------------Conversion to Assembled .fasta Files Finished---------------"
