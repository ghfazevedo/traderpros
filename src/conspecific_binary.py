#!/usr/bin/env python3

import os
import argparse
import math
import numpy as np
import dendropy
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    parser = argparse.ArgumentParser(description="Create plot of conspecificity.")
    parser.add_argument("-ant", "--annotated_tree", required=False, help="Path to tree with speciation events as branch/node annotations (in Nexus format). The annotation should be name branch_speciation_events.")
    parser.add_argument("-brt", "--branch_length_tree", required=False, help="Path to tree with speciation events as branch lengths (Nexus format)")
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
    """Transform branch lengths to binary speciation."""
    for node in tree.preorder_node_iter():
        if node.annotations.get_value("branch_speciation_events") is not None:
            n_speciation = node.annotations.get_value("branch_speciation_events")[0]
            node.annotations.add_new("br_time", node.edge.length)
            node.edge.length = 1.0 if n_speciation != "0" else 0.0

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
    out_path = os.path.join(out_dir, f"{prefix}_binsp.matrix.csv")
    np.savetxt(out_path, matrix, delimiter=",", header=",".join(labels), comments="")
    print(f"Speciation matrix saved to {out_path}")

def plot_heatmap(matrix, labels, out_dir, prefix, figsize=(40, 32), fontsize=10):
    """Creates and saves a binary heatmap with 0 as 'Same Species' (blue) and 1 as 'Different Species' (white)."""
    plt.figure(figsize=figsize)
    
    # Define the binary colormap for 0 (blue) and 1 (white)
    ax = sns.heatmap(
        matrix,
        xticklabels=labels,
        yticklabels=labels,
        cmap=sns.color_palette(["#0D47A1", "#FFFFFF"]),  # Blue for 0, White for 1
        cbar=False,
        square=True,
        linewidths=.05,
        linecolor='gray',
        vmin=0,
        vmax=1
    )
    
    # Set tick labels and font sizes
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=fontsize)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=fontsize)

    
    # Create custom legend with two squares
    same_species_patch = mpatches.Patch(facecolor="#0D47A1", edgecolor="black", label="Same Species")
    diff_species_patch = mpatches.Patch(facecolor="#FFFFFF", edgecolor="black", label="Different Species")
    
   # Add the custom legend to the right of the plot
    plt.legend(
        handles=[same_species_patch, diff_species_patch],
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),  # Position to the right of the heatmap
        borderaxespad=0,
        title="",
        frameon=True
    )

    # Adjust layout to prevent clipping
    plt.tight_layout()

    # Save the figure
    pdf_path = os.path.join(out_dir, f"{prefix}_heatmap.binsp.pdf")
    png_path = os.path.join(out_dir, f"{prefix}_heatmap.binsp.png")
    plt.savefig(pdf_path, format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(png_path, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Binary heatmap saved to {pdf_path} and {png_path}")


def main():
    args = parse_arguments()

    # Validate that at least one of the two arguments is provided
    if not (args.annotated_tree or args.branch_length_tree):
        print("One of -ant or -brt must be provided. These arguments take nexus trees with information on speciation events.")
        exit()

    if (args.annotated_tree and args.branch_length_tree):
        print("Only one of -ant or -brt must be provided. These arguments take nexus trees with information on speciation events.")
        exit()

    # Convert paths to absolute path
    args.out_dir = os.path.abspath(args.out_dir)
    if (args.annotated_tree):
        args.annotated_tree   = os.path.abspath(args.annotated_tree)
    if (args.branch_length_tree):
        args.branch_length_tree   = os.path.abspath(args.branch_length_tree)

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
    if args.annotated_tree:
        sp_events_tree = read_tree(args.annotated_tree)
    else:
        sp_events_tree = read_tree(args.branch_length_tree)

    # Clone and transform the protracted tree if tree is provided as annotations
    conspec_tree = clone_tree(sp_events_tree)
    if args.annotated_tree:
        transform_branch_lengths(conspec_tree)
    
    # Prune to keep only focal taxa if provided
    if args.keep_only_taxa:
        taxa_to_retain = args.keep_only_taxa.split(',')
        prune_tree_to_keep(conspec_tree, taxa_to_retain)
        conspec_tree.update_bipartitions()

    # Calculate distance matrix
    distance_matrix, labels = calculate_distance_matrix(conspec_tree)

    # Save distance matrix
    save_distance_matrix(distance_matrix, labels, args.out_dir, args.prefix)

    # Plot and save heatmap
    plot_heatmap(distance_matrix, labels, args.out_dir, args.prefix, figsize=figsize, fontsize=args.fontsize)

if __name__ == "__main__":
    main()