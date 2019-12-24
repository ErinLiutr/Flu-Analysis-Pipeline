import sys
import pandas as pd

id = sys.argv[1]
metadata = sys.argv[2]
type = sys.argv[3]
reference = sys.argv[4]

# match blinded sample number with subject & sample position
# Usage note!!
# metadata should be in .csv format, with name "metadata.csv"
# sample number with column name: "ID"
# sample group with column name: "Sample Group"
# subject ID with column name: "Subject"
# sampling position with column name: "Sample Type"

def creat_dict():
    df = pd.read_csv("%s"%metadata)
    names = df.set_index('ID').T.to_dict()
    dict = {}
    for key, val in names.items():
        dict[key] = [str(val["Subject"]), str(val["Sample Group"]), val["Sample Type"]]

    return dict

creat_dict()

## read .fasta reference files into a dictionary
## key: segment information
## value: [segment number, sequences]

def create_ref(id):
    ref = {}

    f = open(reference, "r")

    f1 = f.read()
    segments = f1.split(">")

    for i in range(0, len(segments)):
        info = segments[i].split()
        #print info
        length = len(info)
        if length > 0:
            ## concatenate sequences
            for j in range(5, length):
                #print info[j]
                info[5] += info[j]

            ## sequence information
            chromosome = info[0][0:11]

            ## segment number
            info[0] += (" sample " + str(id) + " segment " + str(i))
            intro = ">" + info[0]
            ref[chromosome] = [intro, info[5]]

    return (ref)


## Generate a counts dictionary, mapping the segment name and number of nts of that segment in reference genome
counts = create_ref(id)
for key, val in counts.items():
    counts[key] = len(val[1])

def check_counts(ref):
    remove = False
    for segment in ref.keys():
        count = 0
        for char in ref[segment][1]:
            if char=='N':
                count+=1

        ## delete the assembled the result if it has over 50% nts missing
        if count > counts[segment]*0.5:
            remove = True

    return remove

def coverage(id):
    # return the reference dictionary without the variance
    file = "07-coverage/%s/%s-merged.coverage"%(type, id)

    f = open(file, "r")
    lines = f.readlines()

    ## segments is a dictionary summarizing the positions
    ## key: reference segment name (NC_007371.1) (gb:CY260950|Organism:Influenza)
    ## value: list of tuples (column 2nd and 4th)

    segments = {}
    for line in lines:
        split = line.split()
        if split[0] in segments.keys():
            segments[split[0]].append((split[1], split[3]))
        else:
            segments[split[0]] = []
            segments[split[0]].append((split[1], split[3]))

    ref = create_ref(id)

    ## if segment missing in the coverage file, assembled sequences
    ## only have N's

    for segment in ref.keys():
        ## the string used here (|Organism:Influenza) is specific to the .coverage files
        check = segment + "|Organism:Influenza"

        if check not in segments.keys():
            num = len(ref[segment][1])
            new_seg = ''
            for i in range(0, num):
                new_seg += 'N'
            ref[segment][1] = new_seg


    for segment in segments.keys():

        ## sequences of reference segment
        ref_seg = ref[segment[0:11]][1]

        new_seg = ''
        for i in range(0, len(ref_seg)):
            new_seg += 'N'

        ref_list = list(ref_seg)
        new_list = list(new_seg)

        for pair in segments[segment]:
            if int(pair[1]) > 1:
                new_list[int(pair[0])-1] = ref_list[int(pair[0])-1]

        string = ''.join(new_list)
        ref[segment[0:11]][1] = string

    if not check_counts(ref):
        return ref
    else:
        return {}


def summarize(file):
    summary = {}
    f = open(file, "r")
    f1 = f.readlines()
    for line in f1:
        ## the string used here (gb) is specific to the .coverage files
        if line[0:3] == "gb:":
            info = line.split()
            chromosome = info[0]
            if chromosome in summary:
                # take POS & ALT in .vcf files!!
                summary[chromosome].append((info[1], info[4]))
            else:
                summary[chromosome] = []
                summary[chromosome].append((info[1], info[4]))
    return (summary)

def fasta(id):

    names = creat_dict()

    ref = coverage(id)
    if ref == {}:
        print ("removing..." + id)
        return

    if id in names.keys():
        new_name = "_".join(names[id])
    else:
        new_name = id

    vcf = "06-vcf/%s"%type + id + "-merged.vcf"

    summary = summarize(vcf)

    name = "08-fasta/%s"%type + new_name + ".fasta"
    
    if " " in name:
        name.replace(" ", "_")

    fasta = open(name,"w+")

    # Sort the sequence according to the segment ordering
    sorted_key = sorted(ref.keys(), key = lambda v: int(ref[v][0][-1]))

    # Update indels from the updated .vcf files
    for segment in sorted_key:
        if segment in summary.keys():
            sequence = list(ref[segment][1])
            for update in summary[segment]:
                index = int(update[0])
                # .vcf files start indexing from 1!!
                # minus one to account for this difference
                sequence[index-1] = update[1]
            string = ''.join(sequence)
            header = ref[segment][0].replace(":"," ")
            fasta.write(header)
            fasta.write("\n")
            fasta.write(string)
            fasta.write("\n")
        else:
            fasta.write(ref[segment][0])
            fasta.write("\n")
            fasta.write(ref[segment][1])
            fasta.write("\n")


if __name__== "__main__":
    fasta(id)
