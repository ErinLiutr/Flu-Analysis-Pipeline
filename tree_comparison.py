import dendropy
from dendropy.calculate import treecompare

# groundTruthFile = '/Users/liutianrui/Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/analysis/parsnp_output/PB1/parsnp.tree'#path to ground truth tree, assume it is newick format
# estimationFile = '/Users/liutianrui/Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/analysis/parsnp_output/PB2/parsnp.tree'#path to where your tree is
groundTruthFile = "/Users/liutianrui/Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/Data_DA09252019/segement_genomes_output/parsnp_output_segments/parsnpout_PB2/parsnp.tree"
estimationFile = "/Users/liutianrui/Desktop/lab/flu_transmission/data/Reanalysis_of_EMIT/Data_DA09252019/segement_genomes_output/parsnp_output_segments/parsnpout_PB1/parsnp.tree"

tns = dendropy.TaxonNamespace()
gtTree = dendropy.Tree.get(file=open(groundTruthFile,'r'), schema='newick', taxon_namespace=tns)
estimateTree = dendropy.Tree.get(file=open(estimationFile,'r'), schema='newick', taxon_namespace=tns)

#metrics, weighted RF is unsymmetric, unweighted RF is symmetric distance
weightedRF = treecompare.weighted_robinson_foulds_distance(gtTree, estimateTree)
unweightedRF = treecompare.unweighted_robinson_foulds_distance(gtTree, estimateTree)
euclideanDist = treecompare.euclidean_distance(gtTree, estimateTree)


print weightedRF
print unweightedRF
print euclideanDist