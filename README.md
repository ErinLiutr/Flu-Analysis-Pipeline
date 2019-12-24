**Influenza Transmission Analysis Pipeline**

This is a pipeline built for influenza virus transmission analysis. 
Analysis dataset:
Prometheus & EMIT

**What it does:**
1. Reference guided genome assembly
2. Multiple sequence alignment (MAFFT)
3. Check for homologous recombination signal 
4. Variant analysis
5. Phylogenetic analysis (supports RAxML and Parsnp)
6. Transmission tree analysis (in process)

**Downloads before usage:**
- Bowtie2
- Samtools
- Lofreq
- 3SEQ
- RAxML
- Parsnp

**Working directory structure:**

For genome assembly:
1. 01-trimmomatic (in process)
2. 02-bowtie2-sam
3. 03-bam
4. 04-sorted-bam
5. 05-merged-bam
6. 06-vcf
7. 07-coverage
8. 08-fasta

For multiple sequence alignment:
9. fragments
10. separate_frags
11. concatenated

For homologous recombination signal check:
12. homo_recombination

For variant analysis:
13. variant_analysis

For phylogenetic analysis:
14. raxml_output
15. parsnp_output
16. tree_comparison

For transmission tree analysis: (in process)

**Metadata Example**
Example files are all under data
1. reference genome:
Prometheus: 
- A_Maryland_90_2017.fasta 
- B_DOC_03_2018.fasta
EMIT:
- fluA_ny_h3n2.fna

2. metadata file type 1:
For pairing up the two sequencing results for each  specimen replicate

Prometheus:
- fluA.csv
- fluB.csv
EMIT:
- NA

3. metadata file type 2:
For pairing up samples with subject ID

Prometheus:
- metadataA.csv
- metadataB.csv
EMIT:
- metadataEMIT.csv
