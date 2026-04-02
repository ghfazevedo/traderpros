#!/usr/bin/env python3
"""
Plot rates from traderpros run.
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

def rename_extinction(col, n_obs, n_hid):
    idx = int(col.split('[')[1].rstrip(']')) - 1
    obs = idx % n_obs
    hid = idx // n_obs
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return f"{obs}{letters[hid]}"

def rename_transition(col, n_obs):
    idx = int(col.split('[')[1].rstrip(']')) - 1
    for i in range(n_obs):
        count = n_obs - 1
        if idx < count:
            js = [j for j in range(n_obs) if j != i]
            j = js[idx]
            return f"q{i}{j}"
        idx -= count
    return col  # fallback

def main():
    parser = argparse.ArgumentParser(description="Plot rates from traderpros run.")
    parser.add_argument('-i', '--input', required=True, help="Folder with results or full path to .Rates.log file")
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
            if file.endswith('.Rates.log'):
                input_file = os.path.join(args.input, file)
                break
        else:
            raise FileNotFoundError("No .Rates.log file found in the input folder")
    else:
        raise FileNotFoundError("Input is neither a file nor a directory")

    # Read data
    df = pd.read_csv(input_file, sep='\t')

    # Discard burn-in
    n_samples = len(df)
    burn_in_samples = int(n_samples * args.burn_in)
    df = df.iloc[burn_in_samples:]

    # Extract columns
    extinction_cols = [col for col in df.columns if col.startswith('extinction[')]
    pop_formation_cols = [col for col in df.columns if col.startswith('pop_formation[')]
    spCompletion_cols = [col for col in df.columns if col.startswith('spCompletion[')]
    transition_rates_cols = [col for col in df.columns if col.startswith('transition_rates[')]

    # Rename
    extinction_rename = {col: rename_extinction(col, args.n_observed_states, args.n_hidden_states) for col in extinction_cols}
    pop_formation_rename = {col: rename_extinction(col, args.n_observed_states, args.n_hidden_states) for col in pop_formation_cols}
    spCompletion_rename = {col: rename_extinction(col, args.n_observed_states, args.n_hidden_states) for col in spCompletion_cols}
    transition_rename = {col: rename_transition(col, args.n_observed_states) for col in transition_rates_cols}

    # Prepare data
    def prepare_data(cols, rename_dict):
        if not cols:
            return pd.DataFrame(columns=['Rate', 'Value'])
        data = df[cols].rename(columns=rename_dict)
        melted = data.melt(var_name='Rate', value_name='Value')
        if args.outliers == 'exclude':
            q1 = melted['Value'].quantile(0.25)
            q3 = melted['Value'].quantile(0.75)
            iqr = q3 - q1
            lower = q1 - 1.5 * iqr
            upper = q3 + 1.5 * iqr
            melted = melted[(melted['Value'] >= lower) & (melted['Value'] <= upper)]
        return melted

    trans_data = prepare_data(transition_rates_cols, transition_rename)
    spComp_data = prepare_data(spCompletion_cols, spCompletion_rename)
    popForm_data = prepare_data(pop_formation_cols, pop_formation_rename)
    ext_data = prepare_data(extinction_cols, extinction_rename)

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(21/2.54, args.size/2.54))  # A4 width x height in inches

    plot_order = [
        ('Transition Rates', trans_data, axes[0,0]),
        ('Speciation Completion Rates', spComp_data, axes[0,1]),
        ('Population Formation Rates', popForm_data, axes[1,0]),
        ('Population Extinction Rates', ext_data, axes[1,1])
    ]

    for title, data, ax in plot_order:
        if data.empty:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue
        if args.graph == 'violin':
            sns.violinplot(data=data, x='Rate', y='Value', ax=ax, color='lightgray')
        elif JOYPY_AVAILABLE:
            # Ridge plot using joypy
            joyplot(data, by='Rate', ax=ax, colormap=sns.color_palette("Greys", as_cmap=True), alpha=0.7)
        else:
            # Fallback to violin
            sns.violinplot(data=data, x='Rate', y='Value', ax=ax, color='lightgray')
        ax.set_title(title)
        ax.grid(False)

    plt.tight_layout()

    # Save
    pdf_path = os.path.join(args.output, 'Traderpros_Rates_Plot.pdf')
    tiff_path = os.path.join(args.output, 'Traderpros_Rates_Plot.tiff')
    plt.savefig(pdf_path)
    plt.savefig(tiff_path)
    print(f"Plots saved to {pdf_path} and {tiff_path}")

if __name__ == '__main__':
    main()