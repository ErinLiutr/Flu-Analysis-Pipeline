"""
This is the script for the analysis pipeline
Part 1: Reference guided genome assembly
Part 2: Classification: reformatting FASTA files
Part 3: Recombination Signal Check
Part 4: Variant analysis
Part 5: Multiple sequence alignment with MAFFT
Part 6: Phylogenetic analysis
Part 7: Transmission analysis

NOTE:
1. currently Part1 hasn't been adapted for EMIT
2. hasn't integrated BEAST2 tree generation
3. hasn't integrated part 6 properly

USAGE NOTE:
1. Put all input files under the same directory of this script
2. Softwares used in the pipeline (need to be installed beforehand)
- Bowtie2
- Samtools
- Lofreq
- 3SEQ
- RAxML
- Parsnp
3. Reference genome sequences
"""
import os
import sys

"""
Input to run the pipeline (all strings):
1. The type of the reads (eg. In Prometheus we have UpperA and UpperB)
2. Reference genome filename
    - For Prometheus we used A_Maryland_90_2017.fasta and B_DOC_03_2018.fasta
3. The directory of the reads
4. Metadata file type1 filename
5. Metadata file type2 filename

Note: 
- Metadata file type1 is for pairing up the two sequencing results for each 
specimen replicate (example file: fluA.csv and fluB.csv) 
- Metadata file type2 is for for pairing up samples with subject ID
(example file: metadataA.csv and metadataB.csv)
"""

type = sys.argv[1]
reference = sys.argv[2]
reads = sys.argv[3]
m1 = sys.argv[4]
m2 = sys.argv[5]


"""
Part 1: Reference guided genome assembly
Run read_mapping.py 

Input: the 5 inputs of this pipeline

Output: it'll create 9 directories
references: reference build files for Bowtie2
01-trimmomatic: NA
02-bowtie2-sam: output of Bowtie2 (sam files)
03-bam: bam files (samtools view -Sb)
04-sorted-bam: sorted bam files (samtools sort)
05-merged-bam: merged bam files (samtools merge)
06-vcf: vcf files (lofreq)
07-coverage: coverage files (samtools mpileup)
08-fasta: fasta.py
"""

# run assembly script
os.system("python assembly.py %s %s %s %s"% (type, reference, reads, m1, m2))


"""
Part 2: Classifying .fasta files into 8 segments
For building phylogenetic trees on segments

Input:
1. input_path: the directory under which assembled fasta files are located
(one fasta file per sample)
2. reference_path: the path to the reference genome fasta file
3. output_path: the directory to store the output fragment files
4. output_path2: the directory to store the output fragment direcotries
5. concate_output_path: the path to store the concatenated fasta file

Output: it'll create 3 directories
1. fragments: has 8 fasta files, each containing only sequences of
one of the 8 segments. 

2. separate_frags: has 9 directories.
It has a separate directory for each of the 8 segment that stores
, for each sample, the segment sequence in separate fasta files
It also has a directory called reference that stores 8 fasta files,
each storing sequences from one segment in the reference genome.

3. concatenated: has one .fasta file storing the concatenated sequence
"""

os.mkdir("fragments")
os.mkdir("separate_frags")
os.mkdir("concatenated")
os.system("python classification.py %s %s %s %s %s"%("08-fasta", reference, "fragments", "separate_frags", "concatenated"))



"""
Part 3: Recombination signal check
Input:
1. output_dir: the path to where the output 3 seq reports will be stored

Output: it'll create one directory
1. homo_recombination:
It has a separate directory for each of the 8 segment that stores the 3 seq 
result of that segment
"""
os.mkdir("homo_recombination")
os.system("python homo_recombination.py homo_recombination")


"""
Part 4: Variant analysis

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
any two subjects 
x-axis: subject 1 ID
y-axis: subject 2 ID
"""
os.mkdir("variant_analysis")
os.system("python variant_analysis.py %s %s %s %s"%(m2, "06-vcf", "variant_analysis", reference, [0.05, 0.5]))

"""
Part 5 & 6: Multiple sequence alignment & Phylogenetic tree analysis

Input:
1. threads: the number of threads to be used for reconstructing phylogenetic trees
2. reference_path: path to the reference genome
3. fasta_path: path to the sample sequences (fasta files)
What it does:

Part 1: Multiple sequence alignment
Create two directories:
    -fragments_aligned: 8 fasta files for segment alignment results
    -full_aligned: fasta file for concatenated sequence alignment result

Part 2: Building phylogenetic trees with RAxML
Create two directories:
    -raxml_output: segment trees
    -full_aligned/raxml_output: concatenated tree

Part 3: Building phylogenetic trees with Parsnp
Create one directory:
    -parsnp_output: it has nine folders, one for each of the eight segment
    and one for the reference. These folders all contain a file "parsnp.tree", 
    which is the phylogenetic tree in newick format; the file "parsnp.ggr" can 
    can be viewed using Gingr
    
Part 4: Compare phylogenetic trees
Create one directory 
    -tree_comparison: it contains the tree comparison figures using three metrics
NOTE:
    1. BEAST2 trees haven't been integrated by this script because
    generation of BEAST2 trees haven't been automated yet
"""

# run phylo_tree.py for Part 1 - Part 3
os.system("python phylo_tree.py %s %s %s"%(8, reference, "08-fasta"))

# run tree_comparison.py for Part 4
os.system("python tree_comparison.py")


"""
Part 7: Transmission tree analysis

Input:

Output:
"""

# to be added