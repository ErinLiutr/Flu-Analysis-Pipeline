"""
This is the script for reformatting the assembled fasta files

Input:
1. input_path: the directory under which assembled fasta files are located
(one fasta file per sample)
2. reference_path: the path to the reference genome fasta file
3. output_path: the directory to store the output fragment files
4. output_path2: the directory to store the output fragment direcotries
5. concate_output_path: the path to store the concatenated fasta file

What it does:

Part 1:
Get 8 fasta files, each containing only sequences of
one of the 8 segments. Write the 8 files and put
under output_path

Part 2:
Generate a separate directory for each of the 8 segment that stores
, for each sample, the segment sequence in separate fasta files
Also generate a directory called reference that stores 8 fasta files,
each storing sequences from one segment in the reference genome.
Put the 9 directories under output_path2

Part 3:
Generate the concatenated fasta file.
Put it under concate_output_path

"""
import os
import sys
import shutil

input_path = sys.argv[1]
reference_path = sys.argv[2]
output_path = sys.argv[3]
output_path2 = sys.argv[4]
concate_output_path = sys.argv[5]

map = {1: "PB2", 2: "PB1", 3:"PA", 4: "HA", 5:"NP", 6:"NA", 7:"MP", 8:"NS"}

def get_segments(input_path):
    """
    This function reads all the fasta files in the input path,
    separates sequences of different segments and then
    stores them in a dictionary

    :param input_path: the directory under which all fasta files are stored

    :return: segments, a dictionary where the keys are segment numbers (int)
    and values are the corresponding DNA sequences (string)
    """
    segments = {}
    for file in os.listdir(input_path):
        name = file.split(".fasta")[0]
        lines = ""
        with open(input_path + file, "r") as f:
            for line in f.readlines():
                lines += line.replace("\n","")
            f.close()

        entries = lines.split(">")
        if "" in entries:
            entries.remove("")

        for seg in entries:
            ## Note: the header of assembled fasta files contain
            ## segment information
            s = seg.split("segment ")
            seg_num = int(s[1][0])
            sequence = s[1][1:]
            id = s[0]
            newID = id.split()[-1]
            id = " ".join(id.split()[:-1])
            header = ">" + name + " " + newID + " " + id + "\n" + sequence
            if seg_num in segments.keys():
                segments[seg_num].append(header)
            else:
                segments[seg_num] = [header]
    return segments

def add_ref(reference_path, segments):
    """
    This function adds the segment information of the reference genome

    :param reference_path: the path to the reference genome

    :param segments: output from get_segments

    :return: a new dictionary with the segments of reference genome added
    """
    with open(reference_path, "r") as f:
        entries = f.readlines()
    f.close()

    new_entries = []
    lines = ""
    for line in entries:
        if line[0] == ">":
            if lines != "":
                new_entries.append(lines+"\n")
            header = ">" + "reference " + line.split(">")[1]
            new_entries.append(header)
            lines = ""
        else:
            nl = line.replace("\n","")
            lines += nl
    new_entries.append(lines)

    for idx in range(0, len(new_entries), 2):
        segments[idx/2 + 1].append(new_entries[idx] + new_entries[idx+1])

    return segments

def write(segments, output_path):
    """
    This function writes the segments dictionary into 8 files and
    stores them under the output_path

    :param segments: dictionary where keys are segments numbers and values
    are segments DNA sequences

    :param output_path: the path to store the 8 segment files

    :return: NA
    """
    ## write into files
    for num in range(1, 9):
        segment = map[num]
        with open(output_path + "/%s.fasta"%segment, "wb") as f:
            for lines in segments[num]:
                f.write(lines + "\n")
        f.close()

def separate_frags(input_path, output_path):
    """
    This function generates a separate directory for each of the 8 segment that stores
    , for each sample, the segment sequence in separate fasta files
    It also generates a directory called reference that stores 8 fasta files,
    each storing sequences from one segment in the reference genome.
    Put the 9 directories under output_path2

    :param input_path: path to the assembled fasta files

    :param output_path: the path to store the separate fragments directories

    :return: NA
    USAGE NOTE:
    1. This function must be called after outputs from part 1 have been
    generated and put under output_path
    """
    os.mkdir(output_path)
    for file in os.listdir(input_path):

        ## get the segment name
        id = file.split(".fasta")[0]

        ## case 1: segment files
        if len(id) in [2,3]:
            ## make a separate directory for each of the segment
            os.mkdir(output_path + "/" + id)
            with open(input_path + "/" + file, "r") as f:
                entries = f.readlines()

            ## write a segment file for each sample and put it under
            ## the segment directory
            if "\n" in entries:
                entries.remove("\n")
            for idx in range(0, len(entries), 2):
                filename = entries[idx].split(" ")[0].replace(">","") + ".fasta"
                with open(output_path + "/" + id + "/%s"%filename, "wb") as f2:
                    f2.write(entries[idx])
                    f2.write(entries[idx+1])
                    f2.close()
            f.close()

        ## case 2: reference file
        else:
            ## make a separate directory for the reference genome
            os.mkdir(output_path + "/" + "reference")
            with open(input_path + "/" + file, "r") as f:
                entries = f.readlines()

            ## new_entries store sequences in the following format:
            ## line 1: segment 1 header
            ## line 2: segment 1 sequence
            ## line 3: segment 2 header
            ## ...
            ## NOTE:
            # we added this extra step here because the reference genomes
            # sequences for prometheus are fragmented into multiple lines
            # eg.
            # (line 1)> header here, segment 1
            # (line 2)AACG...
            # (line 3)TCC...

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

            ## write a segment file for each of the eight segments in the reference genome
            for idx in range(0, len(new_entries), 2):
                filename = "reference_" + map[idx/2 + 1] + ".fasta"
                with open(output_path + "/" + "reference" + "/%s"%filename, "wb") as f2:
                    f2.write(new_entries[idx])
                    f2.write(new_entries[idx+1])
                    f2.close()
            f.close()

def concatenate (input_path, reference_path, reference_name, concate_output_path):
    """
    This function writes a fasta file that concatenates all the sample sequences and
    the reference genome
    Each line in the fasta represents either a sample or the reference

    :param input_path: path to the assembled fasta files

    :param reference_path: path to the reference genome

    :param reference_name: a string to indicate the reference genome
    that will be stored as a header in the concatenated fasta file
    eg. B/District Of Columbia/03/2018

    :param concate_output_path: the path to store the concatenated fasta
    :return: NA
    """
    ## generate the concatenated.fasta file
    ## from reading assembled fasta files
    lines = []
    for file in os.listdir(input_path):
        with open(input_path + file, "r") as f:
            read_lines = f.readlines()
            subject_info = file.split(".fasta")[0]
            id = read_lines[0].split(" ")[-3]
            lines.append("\n" + ">" + subject_info + " " + id + "\n")
            for l in read_lines:
                if not l[0] == ">":
                    lines.append(l.replace("\n",""))
        f.close()

    ## write the reference genome
    lines.append("\n>reference %s\n"%reference_name)
    with open(reference_path) as f:
        for line in f.readlines():
            if not line[0] == ">":
                lines.append(line.replace("\n",""))

    with open(concate_output_path + "concatenated.fasta", "wb") as f:
        for idx in range(len(lines)):
            if idx == 0:
                f.write(lines[idx][1:])
            else:
                f.write(lines[idx])
        f.close()


"""
Part 1: 
"""
segments = get_segments(input_path)
segments = add_ref(reference_path, segments)

## first create the directory to store outputs
os.mkdir(output_path)

## write output
write(segments, output_path)

## copy reference genome to the output_path for future use
shutil.copy(reference_path, output_path + "/reference.fasta")

"""
Part 2: 
"""
separate_frags(output_path, output_path2)


"""
Part 3: 
"""
## Use the reference genome file name as the reference name
reference_name = reference_path.split("/")[-1]
concatenate (input_path, reference_path, reference_name, concate_output_path)