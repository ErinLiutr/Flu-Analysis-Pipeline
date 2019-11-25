import os


## Step 1 read mapping
## run read_mapping.py



## Step 2 classify .fasta files into 8 segments
input_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/UpperA"
output_path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments2"

os.mkdir("fragments2")
map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}
segments = {}

for file in os.listdir(input_path):
    name = file.split(".fasta")[0]
    lines = ""
    with open(input_path + "/" + file, "r") as f:
        for line in f.readlines():
            lines += line.replace("\n","")

        f.close()

    entries = lines.split(">")
    if "" in entries:
        entries.remove("")

    for seg in entries:
        s = seg.split("segment ")
        seg_num = int (s[1][0])
        sequence = s[1][1:]
        id = s[0]
        newID = id.split()[-1]
        id = " ".join(id.split()[:-1])
        header = ">" + name + " " + newID + " " + id + "\n" + sequence
        if seg_num in segments.keys():
            segments[seg_num].append(header)
        else:
            segments[seg_num] = [header]


for num in range(1, 9):
    segment = map[num]
    with open(output_path + "/%s.fasta"%segment, "wb") as f:
        for lines in segments[num]:
            f.write(lines + "\n")
    f.close()


## Step 3 Phylogenetic analysis
## run analysis.py



## Step 4 Transmission tree analysis

## Step 5 Transmission bottleneck analysis