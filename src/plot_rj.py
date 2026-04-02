#!/usr/bin/env python3
"""
Plot RJ indicators from traderpros run.
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Plot RJ indicators from traderpros run.")
    parser.add_argument('-i', '--input', required=True, help="Folder with results or full path to .RJ.log file")
    parser.add_argument('-o', '--output', required=True, help="Output folder name")
    parser.add_argument('-s', '--size', type=float, default=29.7, help="Size of plot in cm (A4 height)")
    parser.add_argument('-b', '--burn_in', type=float, default=0.10, help="Burn-in percentage")
    parser.add_argument('-f','--force', choices=['yes', 'no'], default='no', help="If yes, it forces to run without checking the existence of an output folder. It may erase previous results")

    args = parser.parse_args()

    # Check output folder
    if args.force == "no":
        if os.path.exists(args.output):
            response = input(f"Output folder '{args.output}' exists. Continue? (y/n): ")
            if response.lower() != 'y':
                return
        else:
            os.makedirs(args.output)
    else:
        print("Not checking for existence of previous outputs. Results might be erased.")

    # Find input file
    if os.path.isfile(args.input):
        input_file = args.input
    elif os.path.isdir(args.input):
        for file in os.listdir(args.input):
            if file.endswith('.RJ.log'):
                input_file = os.path.join(args.input, file)
                break
        else:
            raise FileNotFoundError("No .RJ.log file found in the input folder")
    else:
        raise FileNotFoundError("Input is neither a file nor a directory")

    # Read data
    df = pd.read_csv(input_file, sep='\t')

    # Discard burn-in
    n_samples = len(df)
    burn_in_samples = int(n_samples * args.burn_in)
    df = df.iloc[burn_in_samples:]

    # Variables of interest
    variables = [
        'is_extinction_state_dependent',
        'is_extinction_hidden',
        'is_pop_formation_state_dependent',
        'is_pop_formation_hidden',
        'is_spCompletion_state_dependent',
        'is_spCompletion_hidden',
        'is_reversible'
    ]

    # Plot: 2 columns x 4 rows
    fig, axes = plt.subplots(4, 2, figsize=(21/2.54, args.size/2.54))  # A4 width x height in inches
    axes = axes.flatten()

    # Titles
    title_map = {
        'is_extinction_state_dependent': 'Extinction State Dependent',
        'is_extinction_hidden': 'Extinction Hidden',
        'is_pop_formation_state_dependent': 'Pop Formation State Dependent',
        'is_pop_formation_hidden': 'Pop Formation Hidden',
        'is_spCompletion_state_dependent': 'SpCompletion State Dependent',
        'is_spCompletion_hidden': 'SpCompletion Hidden',
        'is_reversible': 'Reversible'
    }

    for idx, var in enumerate(variables):
        ax = axes[idx]
        if var not in df.columns:
            ax.text(0.5, 0.5, f'No data for {var}', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title_map.get(var, var))
            continue

        # Compute probabilities
        counts = df[var].value_counts(normalize=True)
        prob_no = counts.get(0, 0)
        prob_yes = counts.get(1, 0)

        bar_labels = ['No', 'Yes']
        bar_values = [prob_no, prob_yes]

        bars = ax.bar(bar_labels, bar_values, color=['lightgray', 'darkgray'], edgecolor='black', alpha=0.8)

        # Evidence lines
        BF_vals = np.array([3.2, 10, 100])
        prior = 0.5
        strength = BF_vals / (BF_vals + 1)
        y_lines = np.concatenate(([prior], strength))
        labels = ["no support", "weak", "substantial", "strong"]

        ax.axhline(prior, linestyle='solid', color='grey', linewidth=1)
        linestyles = ['dotted', 'dashed', (0, (5, 10))]
        for y, ls in zip(strength, linestyles):
            ax.axhline(y, linestyle=ls, color='grey', linewidth=1)

        # Annotate evidence labels
        for y, label in zip(y_lines, labels):
            ax.text(1.05, y, label, ha='left', va='center', color='grey', fontsize=8, transform=ax.get_yaxis_transform())

        # Annotate bar values
        for bar in bars:
            height = bar.get_height()
            ax.annotate(f'{height:.3f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 5),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=10, fontweight='bold')

        ax.set_ylim(0, 1.05)
        ax.set_ylabel('Probability')
        ax.set_title(title_map.get(var, var))

    # Hide the last subplot if 7 variables
    if len(variables) < 8:
        axes[-1].axis('off')

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(args.output, 'Traderpros_RJ_Plot.pdf')
    tiff_path = os.path.join(args.output, 'Traderpros_RJ_Plot.tiff')
    plt.savefig(pdf_path)
    plt.savefig(tiff_path)
    print(f"Plots saved to {pdf_path} and {tiff_path}")

if __name__ == '__main__':
    main()