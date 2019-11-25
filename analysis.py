import os

path = "/Users/liutianrui/Desktop/lab/flu_transmission/data/prometheus/fragments"
threads = 8

os.system("cd ../../Desktop/lab/flu_transmission/data/prometheus")
#1. multiple sequence alignment with MAFFT
print "Multiple sequence alignment"
os.system("mkdir fragments_aligned")

for file in os.listdir(path):
	os.system("mafft --anysymbol %s > fragments_aligned/%s_aligned.fasta"%(file, file))

