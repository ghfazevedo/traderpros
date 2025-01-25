#!/usr/bin/env python3

import os
import argparse
import math
import numpy as np
import dendropy # type: ignore
import seaborn as sns # type: ignore
import matplotlib.pyplot as plt

# Function to confirm with the user if directory exists
def confirm_proceed(message="Directory already exists. Do you want to proceed? (y/n): "):
    while True:
        response = input(message).strip().lower()
        if response == 'y':
            return True
        elif response == 'n':
            print("Exiting program.")
            return False
        else:
            print("Please enter 'y' or 'n'.")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate conspecific probabilities and create heatmap")
    parser.add_argument("-pt", "--protracted_tree", required=True, help="Path to protracted tree (Nexus format)")
    parser.add_argument("-od", "--out_dir", default=".", help="Output directory")
    parser.add_argument("-pre", "--prefix", default="conspecific_probs", help="Prefix for output files")
    parser.add_argument("-tax", "--keep_only_taxa", help="Comma-separated list of focal taxa to plot. At least two should be provided")
    parser.add_argument("-fgs", "--figsize", type=str, default="40,32", help="Figure size as width,height (default: '12,10')")
    parser.add_argument("-fts", "--fontsize", type=int, default=10, help="Font size for labels and ticks (default: 10)")

    return parser.parse_args()

def read_tree(file_path):
    """Reads a Nexus tree file using DendroPy."""
    return dendropy.Tree.get(path=file_path, schema="nexus", preserve_underscores=True)

def clone_tree(tree):
    """Clones a tree object to create a new copy."""
    return tree.clone(depth=1)

def transform_branch_lengths(tree):
    """Transform branch lengths to the natural log of the p_no_speciation annotation."""
    for node in tree.preorder_node_iter():
        if node.annotations.get_value("p_no_speciation") is not None:
            p_no_speciation = float(node.annotations.get_value("p_no_speciation"))
            node.annotations.add_new("br_time", node.edge.length)
            node.edge.length = math.log(p_no_speciation) if p_no_speciation > 0 else math.log(1E-20)

def prune_tree_to_keep(tree, taxa_to_retain_list):
    """Prune tree to keep only the focal taxa provided."""
    # Get the labels of taxa present in the tree
    tree_taxa_labels = [taxon.label for taxon in tree.taxon_set]
    
    # Check if all taxa to retain are present in the tree
    taxa_to_retain_set = set(taxa_to_retain_list)
    missing_taxa = taxa_to_retain_set - set(tree_taxa_labels)

    if missing_taxa:
        print(f"Warning: The following taxa are not present in the tree and will be ignored: {', '.join(missing_taxa)}")
    
    # Prune the tree to keep only the taxa that are present
    taxa_to_retain_list = [taxon for taxon in taxa_to_retain_list if taxon in tree_taxa_labels]
    tree.retain_taxa_with_labels(taxa_to_retain_list)

def calculate_distance_matrix(tree):
    """Calculates a pairwise distance matrix for all tips using DendroPy's phylogenetic_distance_matrix()."""
    # Ensure tree is rooted for distance calculations
    tree.is_rooted = True
    
    tips = tree.leaf_nodes()
    num_tips = len(tips)
    distance_matrix = np.zeros((num_tips, num_tips))
    pdm = tree.phylogenetic_distance_matrix()

    # Get Taxon objects from tips and calculate distances
    taxon_list = [tip.taxon for tip in tips]
    
    for i, taxon1 in enumerate(taxon_list):
        for j, taxon2 in enumerate(taxon_list):
            if i != j:
                try:
                    distance_matrix[i, j] = pdm.distance(taxon1, taxon2)
                except KeyError as e:
                    print(f"Warning: Missing distance for pair ({taxon1}, {taxon2}): {e}")
                    distance_matrix[i, j] = np.nan  # or assign a default value, e.g., 0

    return distance_matrix, [taxon.label for taxon in taxon_list]

def save_distance_matrix(matrix, labels, out_dir, prefix):
    """Saves the distance matrix as a CSV file."""
    out_path = os.path.join(out_dir, f"{prefix}_probs.matrix.csv")
    np.savetxt(out_path, matrix, delimiter=",", header=",".join(labels), comments="")
    print(f"Probability matrix saved to {out_path}")

def plot_heatmap(matrix, labels, out_dir, prefix, figsize=(40, 32), fontsize=10):
    """Creates and saves a heatmap with color-coded annotations."""
    plt.figure(figsize=figsize) 
    
    # Define minimum and maximum values for the log scale
    log_min =  math.log(0.5) #matrix.min() 
    log_max =  math.log(1)    #matrix.max()

    # Create the heatmap with the specified log scale limits
    ax = sns.heatmap(
        matrix,
        xticklabels=labels,
        yticklabels=labels,
        cmap="Blues",
        cbar=True,
        square=True,
        linewidths=.05,
        linecolor='gray',
        vmin=log_min,  # Minimum value for color scale
        vmax=log_max   # Maximum value for color scale
    )
    
    # Set the tick labels and increase font size
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)
    
    # Customize colorbar to show the scale from 0 to 1
    colorbar = ax.collections[0].colorbar
    colorbar_ticks = np.linspace(log_min, log_max, 5)  # 5 tick marks from log_min to log_max
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels([f"{np.exp(tick):.2f}" for tick in colorbar_ticks])  # Convert back to [0,1] scale
    colorbar.ax.tick_params(labelsize=fontsize)
    colorbar.set_label('Prob(conspecificity)', fontsize=fontsize)

    # Optional: Tight layout for better spacing
    plt.tight_layout()
    
    # Save the figure
    pdf_path = os.path.join(out_dir, f"{prefix}_heatmap.probs.pdf")
    png_path = os.path.join(out_dir, f"{prefix}_heatmap.probs.png")
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Heatmap saved to {pdf_path} and {png_path}")


def main():
    args = parse_arguments()

    # Convert out_dir to absolute path
    args.out_dir = os.path.abspath(args.out_dir)
    args.protracted_tree   = os.path.abspath(args.protracted_tree)

    # Check if the directory exists and confirm with the user
    if os.path.exists(args.out_dir):
        print(f"Warning: Directory '{args.out_dir}' already exists. Proceeding may erase previous outputs.")
        if not confirm_proceed():
            exit()

    # Create output directory
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # Convert figsize argument to tuple of floats
    figsize = tuple(map(float, args.figsize.split(',')))

    # Read the trees
    protracted_tree = read_tree(args.protracted_tree)

    # Clone and transform the protracted tree
    conspec_prob_tree = clone_tree(protracted_tree)
    transform_branch_lengths(conspec_prob_tree)
    
    # Prune to keep only focal taxa if provided
    if args.keep_only_taxa:
        taxa_to_retain = args.keep_only_taxa.split(',')
        prune_tree_to_keep(conspec_prob_tree, taxa_to_retain)
        conspec_prob_tree.update_bipartitions()

    # Calculate distance matrix
    distance_matrix, labels = calculate_distance_matrix(conspec_prob_tree)

    # Save distance matrix
    save_distance_matrix(distance_matrix, labels, args.out_dir, args.prefix)

    # Plot and save heatmap
    plot_heatmap(distance_matrix, labels, args.out_dir, args.prefix, figsize=figsize, fontsize=args.fontsize)

if __name__ == "__main__":
    main()