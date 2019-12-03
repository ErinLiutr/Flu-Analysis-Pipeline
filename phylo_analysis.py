import os
import shutil

path = "fragments"
full_path = "UpperB"
threads = 8

#1. multiple sequence alignment with MAFFT
os.system("mkdir fragments_aligned")

for file in os.listdir(path):
	input = path + "/" + file
	os.system("mafft --anysymbol %s > fragments_aligned/%s_aligned.fasta"%(input, file))

os.system("mkdir full_aligned")

## need to automatically generate the concatenated.fasta file
os.system("mafft --anysymbol concatenated.fasta > full_aligned/concatenated_aligned.fasta")

os.chdir("full_aligned")
for file in os.listdir("."):
	print file
	os.system("raxmlHPC -m GTRCAT -p 12345 -s %s -n %s.nexus"%(file, file.split(".fasta")[0]))
os.chdir("../")

# manually removed useless files
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

for file in os.listdir("full_aligned"):
	if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
		os.remove("full_aligned/%s"%file)

for file in os.listdir("."):
	if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
		os.remove("%s"%file)

#add extension to view in dendroscope
os.mkdir("full_aligned/raxml_output")
os.mkdir("full_aligned/raxml_output/bootstrapping")
for file in os.listdir("."):
	if file.startswith("RAxML"):
		os.rename(file, "%s.nexus"%file)
		shutil.move("%s.nexus"%file, "full_aligned/raxml_output/bootstrapping")



# 2. run raxml on the aligned results (8 segments)
os.chdir("fragments_aligned")
for file in os.listdir("."):
	print file
	os.system("raxmlHPC -m GTRCAT -p 12345 -s %s -n %s.nexus"%(file, file.split(".fasta")[0]))
os.chdir("../")

os.mkdir("raxml_output")
for file in os.listdir("fragments_aligned"):
	if file.endswith(".nexus"):
		shutil.move("fragments_aligned/%s"%file, 'raxml_output')

for file in os.listdir("fragments_aligned"):
	if file.endswith(".reduced"):
		os.remove("fragments_aligned/%s"%file)

for file in os.listdir("raxml_output"):
	if file.endswith(".reduced") or file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree"):
		os.remove("raxml_output/%s"%file)


#tree bootstrapping with raxml
for file in os.listdir("fragments_aligned"):
	if file.endswith(".fasta"):
		## simply the output file name
		## only want the segment information!!! + aligned
		new_name = file.split(".fasta_aligned")[0]
		os.system("raxmlHPC -m GTRCAT -p 12345 -# 100 -s fragments_aligned/%s -n %sT1"%(file, new_name))
		os.system("raxmlHPC -m GTRCAT -p 12345 -b 12345 -# 100 -s fragments_aligned/%s -n %sT2"%(file, new_name))
		os.system("raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.%sT1 -z RAxML_bootstrap.%sT2 -n %sT3"%(new_name, new_name, new_name))

		# strict consensus tree
		os.system("raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.%sT2 -n %sT4"%(new_name, new_name))

for file in os.listdir("fragments_aligned"):
	if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
		os.remove("fragments_aligned/%s"%file)

for file in os.listdir("."):
	if file.startswith("RAxML_info") or file.startswith("RAxML_log") or file.startswith("RAxML_result") or file.startswith("RAxML_parsimonyTree") or file.endswith(".reduced"):
		os.remove("%s"%file)

#add extension to view in dendroscope
os.mkdir("raxml_output/bootstrapping")
for file in os.listdir("."):
	if file.startswith("RAxML"):
		os.rename(file, "%s.nexus"%file)
		shutil.move("%s.nexus"%file, "raxml_output/bootstrapping")



#3. run Parsnp (fast2Tree) on 8 segments with reference
# view tree (parsnp.ggr) with Gingr
os.mkdir("parsnp_output")
for dir in os.listdir("separate_frags"):
	if dir != "reference":
		os.mkdir("parsnp_output/%s"%dir)
		reference = "reference_%s.fasta"%dir
		os.system("parsnp -v -x -c -r separate_frags/reference/%s -d separate_frags/%s -C 1000 -o parsnp_output/%s"%(reference, dir, dir))

# Parsnp on concatenated tree
os.system("parsnp -v -x -c -r reference/B_ref.fasta -d UpperB -C 1000 -o parsnp_output/concatenated")






