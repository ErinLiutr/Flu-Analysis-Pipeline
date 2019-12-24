import allel
import os
import pandas as pd
import numpy as np
import csv

#path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/UpperB/vcf/"
path = "/Users/liutianrui/Desktop/vcf/"
output_path = "/Users/liutianrui/Desktop/"
map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}
## getting a list of list

## notice since we're analyzing output files from lofreq, all the entries in vcf
## have passed the filter. This method doesn't consider "QUAL" (quality of each SNP)
def read_file(filename):
    callset = allel.read_vcf(filename,
                             fields=['variants/CHROM', "variants/POS", "variants/REF", "variants/ALT",
                                     'variants/AF']).values()
    ## getting a list fromg each of the following
    ## reformat
    chrome = [str(val) for val in callset[0]]
    reference_nt = [str(val) for val in callset[1]]
    position = [int(val) for val in callset[3]]


    ## getting a list of list fromg each of the following
    ## reformat
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
            if allel_freq[num][idx] < 0.5 and allel_freq[num][idx] > 0.05:
                low_freq_pos.append(num)

    ## sort low frequence variant positions based on Allel frequency (smallest first!)
    frequencies = [min(allel_freq[pos]) for pos in low_freq_pos]
    combined = zip(low_freq_pos, frequencies)
    combined.sort(key=lambda elem: elem[1])
    low_freq_pos = [elem[0] for elem in combined]


    ## calculate the percentage of low frequency variants among all vairants
    percentage = len(low_freq_pos) / float(len(chrome))


    ## indexing in low_freq_pos starts from 0!
    ## initial setting for the output

    #sample_map = {"MH584291": 1, "MH584292": 2, "MH584290":3, "MH584289":4, "MH584288":5, "MH584287":6, "MH584286":7, "MH584285":8}
    #sample_map = {"CY260950": 1, "CY260951": 2, "CY260952": 3, "CY260953": 4, "CY260954": 5, "CY260955": 6,
                  #"CY260956": 7, "CY260957": 8}


    sample_map = {"NC_007373": 1, "NC_007372": 2, "NC_007371":3, "NC_007366":4, "NC_007369":5, "NC_007368":6, "NC_007367":7, "NC_007370":8}

    result = {"Percentage": percentage, "Position":[], "Frequency":[], "Nts":[]}
    for pos in low_freq_pos:
        #start = chrome[pos].index("MH")
        #start = chrome[pos].index("CY")
        start = chrome[pos].index("NC")

        #chrome_name = chrome[pos][start:start+8]

        chrome_name = chrome[pos][start:start + 9]
        segment_num = sample_map[chrome_name]

        position_name = "Segment " + str(segment_num) + " " + map[segment_num] + " " + str(position[pos])
        result["Position"].append(position_name)
        result["Frequency"].append(allel_freq[pos])
        result["Nts"].append((reference_nt[pos], alternate_nt[pos]))

    return result

def convertFileName(metadata, filename):
    ## convert the filename to a string containing subject information
    df = pd.read_csv(metadata)
    names = df.set_index('ID').T.to_dict()
    dict = {}
    for key, val in names.items():

        #dict[key] = [str(val["Subject"]), str(val["Sample Group"]), val["Sample Type"]]
        new_key = key.split("_")[-1]
        dict[new_key] = [str(val["Subject"]), str(val["Sample Group"]), val["Sample Type"]]

    for key in dict.keys():
        if isinstance(key, str) and key in filename:
            print "t"
            return " ".join(reversed(dict[key]))

    return ""

def getSummaryCSV(metadata, path):
    ## get all possible variant positions
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
        #print subject_info, column
    #print summary
    df = pd.DataFrame(data=summary)
    df.to_csv(output_path + "Low Frequency Variants.csv")

    return summary

def getComparison(metadata, path):
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
                        #print idx, key1, key2
                column.append(shared_count)
        comparison[key1] = column

    comparison.pop("Position")
    df = pd.DataFrame(data=comparison)
    df.to_csv(output_path + "Low Frequency Variants Comparison.csv")

getSummaryCSV("/Users/liutianrui/Desktop/metadataEMIT.csv", path)
getComparison("/Users/liutianrui/Desktop/metadataEMIT.csv", path)