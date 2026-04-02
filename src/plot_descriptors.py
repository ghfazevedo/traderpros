#!/usr/bin/env python3
"""
Plot descriptors from traderpros run.
"""

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
try:
    from joypy import joyplot
    JOYPY_AVAILABLE = True
except ImportError:
    JOYPY_AVAILABLE = False
    print("joypy not available, ridge plots will use violin instead.")

def rename_descriptor(col, n_obs, n_hid):
    idx = int(col.split('[')[1].rstrip(']')) - 1
    obs = idx % n_obs
    hid = idx // n_obs
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return f"{obs}{letters[hid]}"

def main():
    parser = argparse.ArgumentParser(description="Plot descriptors from traderpros run.")
    parser.add_argument('-i', '--input', required=True, help="Folder with results or full path to .Descriptors.log file")
    parser.add_argument('-o', '--output', required=True, help="Output folder name")
    parser.add_argument('-g', '--graph', choices=['violin', 'ridge'], default='violin', help="Type of plot")
    parser.add_argument('-s', '--size', type=float, default=29.7, help="Size of plot in cm (A4 height)")
    parser.add_argument('-n', '--n_observed_states', type=int, default=2, help="Number of observed states")
    parser.add_argument('-H', '--n_hidden_states', type=int, default=2, help="Number of hidden states")
    parser.add_argument('-b', '--burn_in', type=float, default=0.10, help="Burn-in percentage")
    parser.add_argument('--outliers', choices=['include', 'exclude'], default='include', help="Include or exclude outliers")
    parser.add_argument('-f','--force', choices=['yes', 'no'], default='no', help="If yes, it forces to run without checking the existence of an output folder. It may erase previous results")

    args = parser.parse_args()

    # Check output folder
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
            if file.endswith('.Descriptors.log'):
                input_file = os.path.join(args.input, file)
                break
        else:
            raise FileNotFoundError("No .Descriptors.log file found in the input folder")
    else:
        raise FileNotFoundError("Input is neither a file nor a directory")

    # Read data
    df = pd.read_csv(input_file, sep='\t')

    # Discard burn-in
    n_samples = len(df)
    burn_in_samples = int(n_samples * args.burn_in)
    df = df.iloc[burn_in_samples:]

    # Extract columns
    ephemerality_cols = [col for col in df.columns if col.startswith('ephemerality[')]
    evolvability_cols = [col for col in df.columns if col.startswith('evolvability[')]
    net_pop_formation_rate_cols = [col for col in df.columns if col.startswith('net_pop_formation_rate[')]
    pop_turnover_cols = [col for col in df.columns if col.startswith('pop_turnover[')]

    # Rename
    ephemerality_rename = {col: rename_descriptor(col, args.n_observed_states, args.n_hidden_states) for col in ephemerality_cols}
    evolvability_rename = {col: rename_descriptor(col, args.n_observed_states, args.n_hidden_states) for col in evolvability_cols}
    net_pop_formation_rate_rename = {col: rename_descriptor(col, args.n_observed_states, args.n_hidden_states) for col in net_pop_formation_rate_cols}
    pop_turnover_rename = {col: rename_descriptor(col, args.n_observed_states, args.n_hidden_states) for col in pop_turnover_cols}

    # Prepare data
    def prepare_data(cols, rename_dict):
        if not cols:
            return pd.DataFrame(columns=['Descriptor', 'Value'])
        data = df[cols].rename(columns=rename_dict)
        melted = data.melt(var_name='Descriptor', value_name='Value')
        if args.outliers == 'exclude':
            q1 = melted['Value'].quantile(0.25)
            q3 = melted['Value'].quantile(0.75)
            iqr = q3 - q1
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            melted = melted[(melted['Value'] >= lower) & (melted['Value'] <= upper)]
        return melted

    ephem_data = prepare_data(ephemerality_cols, ephemerality_rename)
    evol_data = prepare_data(evolvability_cols, evolvability_rename)
    net_pop_data = prepare_data(net_pop_formation_rate_cols, net_pop_formation_rate_rename)
    pop_turn_data = prepare_data(pop_turnover_cols, pop_turnover_rename)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(21/2.54, args.size/2.54))  # A4 width x height in inches

    plot_order = [
        ('Ephemerality', ephem_data, axes[0,0]),
        ('Evolvability', evol_data, axes[0,1]),
        ('Net Population Formation Rate', net_pop_data, axes[1,0]),
        ('Population Turnover', pop_turn_data, axes[1,1])
    ]

    for title, data, ax in plot_order:
        if data.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue
        if args.graph == 'violin':
            sns.violinplot(data=data, x='Descriptor', y='Value', ax=ax, color='lightgray')
        elif JOYPY_AVAILABLE:
            # Ridge plot using joypy
            joyplot(data, by='Descriptor', ax=ax, colormap=sns.color_palette("Greys", as_cmap=True), alpha=0.7)
        else:
            # Fallback to violin
            sns.violinplot(data=data, x='Descriptor', y='Value', ax=ax, color='lightgray')
        ax.set_title(title)
        ax.grid(False)

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(args.output, 'Traderpros_Descriptors_Plot.pdf')
    tiff_path = os.path.join(args.output, 'Traderpros_Descriptors_Plot.tiff')
    plt.savefig(pdf_path)
    plt.savefig(tiff_path)
    print(f"Plots saved to {pdf_path} and {tiff_path}")

if __name__ == '__main__':
    main()