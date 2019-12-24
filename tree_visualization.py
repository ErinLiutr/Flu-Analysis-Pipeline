# Visualize all the phylogenetic trees

import matplotlib as plt
import Phylo

def treeImage(self, newick, rooted=False, outgroup=False):
    """
        Given a newick string, creates an image of the tree.
        Used in L Statistic GUI.
    """
    plt.figure(figsize=(8, 4))
    plt.axis('off')
    ax = plt.subplot(1, 1, 1)
    ax.axis('off')

    # Create the tree object
    tree = Phylo.read(newick, "newick")
    tree.rooted = rooted

    if rooted:
        tree.root_with_outgroup(outgroup)

    # Create the tree image
    Phylo.draw(tree, axes=ax, do_show=False)
    plt.savefig('test.png')
