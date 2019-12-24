"""
This is the script for the analysis pipeline
Part 1: Reference guided genome assembly
Part 2: Classification: reformatting FASTA files
Part 3: Recombination Signal Check
Part 4: Variant analysis
Part 5: Phylogenetic analysis
Part 6: Transmission analysis

NOTE:
1. currently Part1 hasn't been adapted for EMIT
2. hasn't integrated BEAST2 tree generation
3. hasn't integrated part 6 properly

USAGE NOTE:
1. Put all input files under the same directory of this script
2. softwares used in the pipeline (need to be installed beforehand)
- Bowtie2
- Samtools
- Lofreq
- RAxML
- Parsnp
"""

import os
import sys

"""
Input to run the pipeline (all strings):
1. The type of the reads (eg. In Prometheus we have UpperA and UpperB)
2. Reference genome filename
    - For Prometheus we used A_Maryland_90_2017.fasta and B_DOC_03_2018.fasta
3. The directory of the reads
4. Metadata file 1 filename
5. Metadata file 2 filename

Note: 
- Metadata file 1 is for pairing up the two sequencing results for each 
specimen replicate (example file: fluA.csv and fluB.csv) 
- Metadata file 2 is for for pairing up samples with subject ID
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

## run assembly script
os.system("python assembly.py %s %s %s %s"% (type, reference, reads, m1, m2))


"""
Part 2: Classifying .fasta files into 8 segments
For building phylogenetic trees on segments

Input:

Output:
"""




"""
Part 3: Recombination signal check

Input:

Output:
"""



"""
Part 4: Variant analysis

Input:

Output:
"""


"""
Part 5: Phylogenetic tree analysis

Input:

Output:
"""


"""
Part 6: Transmission tree analysis

Input:

Output:
"""

