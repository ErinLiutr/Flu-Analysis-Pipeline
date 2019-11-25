#!/bin/sh

FRAGMENTS = "/Users/liutianrui/Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/analysis/fragments"
THREADS = 8

#1. multiple sequence alignment with MAFFT
echo "Multiple sequence alignment"
mkdir fragments_aligned
for file in ${FRAGMENTS}*.fasta; do
	mafft --anysymbol $file > fragments_aligned/$file"_aligned.fasta"
done

#2. run raxml on the aligned results
cd fragments_aligned
for file in *.fasta; do
	raxmlHPC -m GTRCAT -p 12345 -s $file -n $file.nexus
done
cd ..
mkdir raxml_output
mv fragments_aligned/*.nexus raxml_output
rm fragments_aligned/*.reduced

tree bootstrapping with raxml
cd fragments_aligned

for file in *.fasta; do
	raxmlHPC -m GTRCAT -p 12345 -# 100 -s $file -n ${file%_aligned.fasta}T1
	raxmlHPC -m GTRCAT -p 12345 -b 12345 -# 100 -s $file -n ${file%_aligned.fasta}T2
	raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.${file%_aligned.fasta}T1 -z RAxML_bootstrap.${file%_aligned.fasta}T2 -n ${file%_aligned.fasta}T3

	# strict consensus tree
	raxmlHPC -m GTRCAT -J STRICT -z RAxML_bootstrap.${file%_aligned.fasta}T2 -n ${file%_aligned.fasta}T4
	rm RAxML_info*
	rm RAxML_log*
	rm RAxML_result*
	rm RAxML_parsimonyTree*
	rm *.reduced
done

#add extension to view in dendroscope
for file in RAxML*; do
	mv $file $file.nexus
done

# move to folder
cd ../raxml_output
mkdir bootstrapping
cd ..
mv fragments_aligned/*.nexus raxml_output/bootstrapping

#3. run Parsnp (fast2Tree) on 8 segments with reference
view tree (parsnp.ggr) with Gingr
mkdir parsnp_output
cd separateFrags
for d in */ ; do
	parsnp -v -x -c -r ../reference/reference_${d%/}.fasta -d ${d%/} -C 1000 -o ../parsnp_output/${d%/};
done

# Parsnp on concatenated tree
parsnp -v -x -c -r reference/fluA_ny_h3n2.fna -d full_genomes -C 1000 -o parsnp_output/concatenated;

