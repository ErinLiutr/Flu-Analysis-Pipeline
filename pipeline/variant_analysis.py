"""
This is the script for generating low frequency variants

Input:
1. metadata: path to metadata type2
2. vcf_path: path to the directory storing all vcf files
3. output_path: path to which the variant summary will be stored
4. reference_path: path to the reference genome
5. threshod: a list containing two floating point numbers (between
0 and 1). The first is the lower bound and the second the upper bound
to be identified as a low frequency variants

Output: it'll create a directory called variant analysis that
stores two csv files
1. Low Frequency Variants.csv:
This sheet contains all the low frequency variants found within the threshold.
x-axis: subject ID
y-axis: position of the variants

2. Low Frequency Variants Comparison.csv
This sheet contains the number of shared low frequency variants found between
any two subject
x-axis: subject 1 ID
y-axis: subject 2 ID

NOTE:
Notice since we're analyzing output files from lofreq, all the entries in vcf
have passed the filter. This method doesn't consider "QUAL" (quality of each variation)
"""
import allel
import os
import pandas as pd
import numpy as np
import sys


metadata = sys.argv[1]
vcf_path = sys.argv[2]
output_path = sys.argv[3]
reference_path = sys.argv[4]
threshold = sys.argv[5]

map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}

def get_sample_map(reference):
    """
    This function gets a mapping between segment identifier and segment number
    :param reference: path to the reference genome
    :return: a dictionary where the keys are the segment identifier (the first part of each segment header in
    reference fasta file) and values are the corresponding segment number (1-8)

    NOTE:
    example sample map:
    Prometheus flu A = {"CY260950": 1, "CY260951": 2, "CY260952": 3, "CY260953": 4, "CY260954": 5, "CY260955": 6,
                  "CY260956": 7, "CY260957": 8}
    Prometheus flu B{"MH584291": 1, "MH584292": 2, "MH584290":3, "MH584289":4, "MH584288":5, "MH584287":6, "MH584286":7, "MH584285":8}
    EMIT = {"NC_007373": 1, "NC_007372": 2, "NC_007371":3, "NC_007366":4, "NC_007369":5, "NC_007368":6, "NC_007367":7, "NC_007370":8}
    """
    sample_map = {}

    ## extract the headers. Should have exactly 8 of them
    headers = []
    with open(reference, "r") as f:
        for line in f.readlines():
            if line[0] == ">":
                headers.append(line)

    for header in headers:
        seg_id = header.replace(">", "").split(" ")[0]
        seg_num = int(header.split("segment")[1][1])
        sample_map[seg_id] = seg_num

    return sample_map

def read_file(filename):
    """
    This function reads a vcf file and summarizes its low frequency variants
    :param filename: vcf file name
    :return: a dictionary with four keys
    1. Percentage: value is the the percentage of low frequency variants among all variants
    2. Position: value is a list containing the positions of variants
    (eg. Segment 1 PB2 32)
    3. Frequency: value is a list containing the allele frequencies of variants
    4. Nts: value is a list containing a list in the form [a, b]
    where a is the nt at the position in the reference genome
    b is the nt at the position in the reads (alternate allele)
    """
    sample_map = get_sample_map(reference_path)
    callset = allel.read_vcf(filename,
                             fields=['variants/CHROM', "variants/POS", "variants/REF", "variants/ALT",
                                     'variants/AF']).values()

    ## getting a list from each of the following
    chrome = [str(val) for val in callset[0]]
    reference_nt = [str(val) for val in callset[1]]
    position = [int(val) for val in callset[3]]

    ## these two have a slightly different format
    afreq = [list(val) for val in callset[2]]
    allel_freq = []
    for num in range(len(afreq)):
        entry = []
        for item in afreq[num]:
            if not np.isnan(item):
                entry.append(item)
        allel_freq.append(entry)

    ALnt = [list(val) for val in callset[4]]
    alternate_nt = []
    for num in range(len(ALnt)):
        entry = []
        for item in ALnt[num]:
            if not item == "":
                entry.append(str(item))
        alternate_nt.append(entry)

    ## get low frequence variants (allel frequence below 50%)
    low_freq_pos = []
    for num in range(len(chrome)):
        for idx in range(len(alternate_nt[num])):
            if allel_freq[num][idx] < threshold[1] and allel_freq[num][idx] > threshold[0]:
                low_freq_pos.append(num)

    ## sort low frequence variant positions based on Allel frequency (smallest first!)
    frequencies = [min(allel_freq[pos]) for pos in low_freq_pos]
    combined = zip(low_freq_pos, frequencies)
    combined.sort(key=lambda elem: elem[1])
    low_freq_pos = [elem[0] for elem in combined]


    ## calculate the percentage of low frequency variants among all variants
    percentage = len(low_freq_pos) / float(len(chrome))


    ## indexing in low_freq_pos starts from 0!
    ## initial setting for the output

    result = {"Percentage": percentage, "Position":[], "Frequency":[], "Nts":[]}
    for pos in low_freq_pos:
        ## NOTE: we want to extract the segment where this low frequency variant
        ## is located. chrome_name is the segment identifier
        chrome_name = chrome[pos]
        segment_num = sample_map[chrome_name]

        position_name = "Segment " + str(segment_num) + " " + map[segment_num] + " " + str(position[pos])
        result["Position"].append(position_name)
        result["Frequency"].append(allel_freq[pos])
        result["Nts"].append((reference_nt[pos], alternate_nt[pos]))

    return result

def convertFileName(metadata, filename):
    """
    This function converts the filename (from the name of the vcf files) to a string
    containing subject information
    :param metadata: path to metadata type 2
    :param filename: original filename
    :return: a new filename that contains the subject ID and sampling position;
    empty string if the filename is not found in metadata sheet
    """
    df = pd.read_csv(metadata)
    names = df.set_index('ID').T.to_dict()
    dict = {}
    for key, val in names.items():
        new_key = key.split("_")[-1]
        dict[new_key] = [str(val["Subject"]), str(val["Sample Group"]), val["Sample Type"]]

    for key in dict.keys():
        if isinstance(key, str) and key in filename:
            return " ".join(reversed(dict[key]))

    return ""

def getSummaryCSV(metadata, path, output_path):
    """
    This function gets the positions of all variants.
    Besides what it returns, it also generates a sheet that contains
    all the low frequency variants found within the threshold with format:
    x-axis: subject ID
    y-axis: position of the variants

    :param metadata: path to metadata type 2
    :param path: path to the directory storing all vcf files that'll be summarized
    :param output_path: path to which the variant summary sheet will be stored
    :return: a dictionary where keys are subject information and values are lists
    that contain variants information
    """

    pos = set()
    for file in os.listdir(path):
        result = read_file(path + file)
        pos = pos.union(set(result["Position"]))

    summary = {"Position":list(pos)}
    for file in os.listdir(path):
        result = read_file(path + file)
        subject_info = convertFileName(metadata, file)
        column = []
        for p in pos:
            if p in result["Position"]:
                idx = result["Position"].index(p)
                frequency = result["Frequency"][idx]
                nts = result["Nts"][idx]
                column.append("".join([str(freq) for freq in frequency]) + " " + nts[0] + " " + "".join(nts[1]))
            else:
                column.append("")
        summary[subject_info] = column

    df = pd.DataFrame(data=summary)
    df.to_csv(output_path + "Low Frequency Variants.csv")

    return summary

def getComparison(metadata, path, output_path):
    """
    This function gets all shared low frequency variants
    It generates a sheet that contains the number of shared low frequency variants
    found between any two subjects
    x-axis: subject 1 ID
    y-axis: subject 2 ID

    :param metadata: path to metadata type 2
    :param path: path to the directory storing all vcf files that'll be summarized
    :param output_path: path to which the variant summary sheet will be stored
    :return: NA
    """
    summary = getSummaryCSV(metadata, path)
    comparison = {"Sample": [key for key in summary.keys() if key!= "Position"]}
    for key1 in summary.keys():
        column = []
        for key2 in summary.keys():
            if key1 != "Position" and key2 != "Position":
                ## check how many variants the two samples share
                shared_count = 0
                for idx in range(len(summary[key1])):
                    if summary[key1][idx] != "" and summary[key2][idx] != "":
                        shared_count += 1
                column.append(shared_count)
        comparison[key1] = column

    comparison.pop("Position")
    df = pd.DataFrame(data=comparison)
    df.to_csv(output_path + "Low Frequency Variants Comparison.csv")


getSummaryCSV(metadata, vcf_path, output_path)
getComparison(metadata, vcf_path, output_path)