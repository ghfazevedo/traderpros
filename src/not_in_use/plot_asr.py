#!/usr/bin/env python3
"""
Plot ancestral state reconstruction from traderpros run.
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import toytree

def rename_state(state, n_obs, n_hid):
    if isinstance(state, str):
        try:
            state = float(state)
        except:
            return state
    state = int(state)
    idx = state
    obs = idx % n_obs
    hid = idx // n_obs
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return f"{obs}{letters[hid]}"

def main():
    parser = argparse.ArgumentParser(description="Plot ancestral state reconstruction from traderpros run.")
    parser.add_argument('-i', '--input', required=True, help="Folder with Traderpros results")
    parser.add_argument('-o', '--output', required=True, help="Output folder name")
    parser.add_argument('-s', '--size', type=float, default=29.7, help="Size of plot in cm (A4 height)")
    parser.add_argument('-n', '--n_observed_states', type=int, default=2, help="Number of observed states")
    parser.add_argument('-H', '--n_hidden_states', type=int, default=2, help="Number of hidden states")
    parser.add_argument('-t', '--type', choices=['pie', 'dot'], default='pie', help="Type of plot")
    parser.add_argument('-c', '--colors', default='colorblind', help="Seaborn color palette")
    parser.add_argument('-l', '--layout', choices=['right', 'r', 'left', 'l', 'down', 'd', 'up', 'u', 'unrooted', 'unr', 'circular', 'c'], default='right', help="Tree layout")

    args = parser.parse_args()

    # Check output folder
    if os.path.exists(args.output):
        response = input(f"Output folder '{args.output}' exists. Continue? (y/n): ")
        if response.lower() != 'y':
            return
    else:
        os.makedirs(args.output)

    # Find files
    cond_file = None
    marg_file = None
    for file in os.listdir(args.input):
        if file.endswith('.traits.MAP.cond.tree'):
            cond_file = os.path.join(args.input, file)
        elif file.endswith('.traits.MAP.marg.tree'):
            marg_file = os.path.join(args.input, file)

    if not cond_file or not marg_file:
        raise FileNotFoundError("Missing .traits.MAP.cond.tree or .traits.MAP.marg.tree")

    # Read trees with toytree, which parses Nexus annotations into node attributes
    cond_tree = toytree.tree(cond_file)
    marg_tree = toytree.tree(marg_file)

    # Collect node data from toytree node attributes
    def get_node_data(tree):
        data = {}
        for node in tree.treenode.traverse():
            node_data = {}
            for key in node.__dict__:
                if key.startswith('anc_state'):
                    value = getattr(node, key)
                    try:
                        value = float(value)
                    except:
                        pass
                    node_data[key] = value
            data[node] = node_data
        return data

    cond_node_data = get_node_data(cond_tree)
    marg_node_data = get_node_data(marg_tree)

    # Prepare data for plotting
    # States
    total_states = args.n_observed_states * args.n_hidden_states
    states = list(range(total_states))
    renamed_states = [rename_state(s, args.n_observed_states, args.n_hidden_states) for s in states]

    # Colors
    palette = sns.color_palette(args.colors, total_states)
    color_dict = dict(zip(renamed_states, palette))

    # Layout
    layout_dict = {
        'right': 'r', 'r': 'r',
        'left': 'l', 'l': 'l',
        'down': 'd', 'd': 'd',
        'up': 'u', 'u': 'u',
        'unrooted': 'unr', 'unr': 'unr',
        'circular': 'c', 'c': 'c'
    }
    layout = layout_dict[args.layout]

    # Plot
    fig, axes = plt.subplots(1, 2, figsize=(args.size/2.54 * 2, args.size/2.54))

    # Conditional
    ax1 = axes[0]
    if args.type == 'pie':
        # Prepare pie data: dict of node idx to {state: prob}
        pie_data = {}
        for node, data in cond_node_data.items():
            probs = {}
            for k, v in data.items():
                if k.startswith('anc_state_') and not k.endswith('_pp'):
                    pp_key = k + '_pp'
                    if pp_key in data:
                        state = v
                        renamed = rename_state(state, args.n_observed_states, args.n_hidden_states)
                        probs[renamed] = data[pp_key]
            if probs:
                pie_data[node.idx] = probs
        # Draw tree with pies
        cond_tree.draw(layout=layout, axes=ax1, node_pie=pie_data, pie_colors=color_dict)
    elif args.type == 'dot':
        # Most probable state
        node_colors = {}
        node_sizes = {}
        for node, data in cond_node_data.items():
            if 'anc_state_1' in data:
                state = data['anc_state_1']
                pp = data.get('anc_state_1_pp', 1)
                renamed = rename_state(state, args.n_observed_states, args.n_hidden_states)
                node_colors[node.idx] = color_dict.get(renamed, 'black')
                node_sizes[node.idx] = pp * 20  # scale
        cond_tree.draw(layout=layout, axes=ax1, node_colors=node_colors, node_sizes=node_sizes, node_markers='o')

    ax1.set_title('Conditional MAP')

    # Marginal
    ax2 = axes[1]
    if args.type == 'pie':
        pie_data = {}
        for node, data in marg_node_data.items():
            probs = {}
            for k, v in data.items():
                if k.startswith('anc_state_') and not k.endswith('_pp'):
                    pp_key = k + '_pp'
                    if pp_key in data:
                        state = v
                        renamed = rename_state(state, args.n_observed_states, args.n_hidden_states)
                        probs[renamed] = data[pp_key]
            if probs:
                pie_data[node.idx] = probs
        marg_tree.draw(layout=layout, axes=ax2, node_pie=pie_data, pie_colors=color_dict)
    elif args.type == 'dot':
        node_colors = {}
        node_sizes = {}
        for node, data in marg_node_data.items():
            if 'anc_state_1' in data:
                state = data['anc_state_1']
                pp = data.get('anc_state_1_pp', 1)
                renamed = rename_state(state, args.n_observed_states, args.n_hidden_states)
                node_colors[node.idx] = color_dict.get(renamed, 'black')
                node_sizes[node.idx] = pp * 20
        marg_tree.draw(layout=layout, axes=ax2, node_colors=node_colors, node_sizes=node_sizes, node_markers='o')

    ax2.set_title('Marginal MAP')

    # Legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[s], markersize=10, label=s) for s in renamed_states]
    fig.legend(handles=handles, loc='center right')

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(args.output, 'Traderpros_ASR_Plot.pdf')
    tiff_path = os.path.join(args.output, 'Traderpros_ASR_Plot.tiff')
    plt.savefig(pdf_path)
    plt.savefig(tiff_path)
    print(f"Plots saved to {pdf_path} and {tiff_path}")

if __name__ == '__main__':
    main()