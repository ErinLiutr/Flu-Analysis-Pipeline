"""
This script compare the various phylogenetic trees by calculating
their unweighted Robinson Foulds, weighted Robinson Foulds, and
Euclidean distances.

Input: NA

What it does:
    Create a directory called tree_comparison
    Generate tree comparison figures using three metrics and save
    the figures under tree_comparison
NOTE:
    1. The first two functions are based on internet resources
    2. BEAST2 trees haven't been integrated by this script because
    generation of BEAST2 trees haven't been automated yet
"""

import dendropy
from dendropy.calculate import treecompare
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import os


def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def get_file(software, segment):
    """
    This function outputs the path to a segment tree generated by a certain software
    :param software: software name
    :param segment: segment name (can also be concatenated)
    :return: path to the tree
    """
    if software == "parsnp":
        return "parsnp_output/%s/parsnp.tree"%segment

    # note we store parsnp output of segment trees and concatenated
    # trees under two different directories
    if software == "raxml":
        if segment != "concatenated":
            # we use the bootstrapped best tree!
            return "raxml_output/bootstrapping/RAxML_bestTree.%sT1.nexus"% segment
        else:
            return "full_aligned/raxml_output/bootstrapping/RAxML_bestTree.concatenated_aligned.fastaT1.nexus"

def get_figure(software1, software2, output_path):
    """
    This function generates a comparison figure of trees generated from software 1 and software 2
    :param software1: name of the software
    :param software2: name of the software
    :param output_path: path to where the output figure will be saved
    :return: NA
    """
    map = {1: "PB2", 2: "PB1", 3: "PA", 4: "HA", 5: "NP", 6: "NA", 7: "MP", 8: "NS", 9: "concatenated"}

    wRF = []
    uwRF = []
    eD = []


    for num1 in range(1,10):
        w = []
        u = []
        e = []
        for num2 in range(1, 10):
            segment1 = map[num1]
            segment2 = map[num2]

            groundTruthFile = get_file(software1, segment1)
            estimationFile = get_file(software2, segment2)

            tns = dendropy.TaxonNamespace()
            gtTree = dendropy.Tree.get(file=open(groundTruthFile, 'r'), schema='newick', taxon_namespace=tns)
            estimateTree = dendropy.Tree.get(file=open(estimationFile, 'r'), schema='newick', taxon_namespace=tns)

            # metrics, weighted RF is unsymmetric, unweighted RF is symmetric distance
            weightedRF = treecompare.weighted_robinson_foulds_distance(gtTree, estimateTree)
            unweightedRF = treecompare.unweighted_robinson_foulds_distance(gtTree, estimateTree)
            euclideanDist = treecompare.euclidean_distance(gtTree, estimateTree)
            w.append(weightedRF)
            u.append(unweightedRF)
            e.append(euclideanDist)

        wRF.append(w)
        uwRF.append(u)
        eD.append(e)

    wRF = np.array(wRF)
    uwRF = np.array(uwRF)
    eD = np.array(eD)

    metric_map = {"Weighted Robinson Foulds": wRF, "Unweighted Robinson Foulds": uwRF, "Euclidean Distances": eD}
    for metric in ["Weighted Robinson Foulds", "Unweighted Robinson Foulds", "Euclidean Distances"]:
        fig, ax = plt.subplots()

        im, cbar = heatmap(metric_map[metric], software1, software2, ax=ax,
                           cmap="YlGn", cbarlabel="Distance")

        texts = annotate_heatmap(im, valfmt="{x:.2f}")

        title = "%s on %s and %s Tree"%(metric, software1.capitalize(), software2.capitalize())
        ax.set_title(title, pad = -330)

        fig.tight_layout()

        # save figure to output path
        plt.savefig(output_path)

        #plt.show()


os.mkdir("tree_comparison")

get_figure("parsnp", "parsnp", "tree_comparison/parsnp_trees.png")
get_figure("raxml", "raxml", "tree_comparison/raxml_trees.png")
get_figure("parsnp", "raxml", "tree_comparison/parsnp_raxml_trees.png")