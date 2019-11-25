import os
import shutil

genomes_dir = "../full_genomes"

subject_dict = {}
sample_file = {}

for filename in os.listdir(genomes_dir):
    id = filename.split("_")[0]
    sample_type = filename.split("_")[1]
    # classify files based on id
    if id in subject_dict.keys():
        subject_dict[id].append(sample_type)
    else:
        subject_dict[id] = [sample_type]
    # classify files based on sample location
    if sample_type in sample_file.keys():
        sample_file[sample_type].append(filename)
    else:
        sample_file[sample_type] = [filename]

# print len(subject_dict)
# print subject_dict
# print len(sample_file)
# print sample_file
for s_type in sample_file.keys():
    os.mkdir("../%s"%s_type)
    for f in sample_file[s_type]:
        shutil.copy("../full_genomes/%s"%f, "../%s"%s_type)
