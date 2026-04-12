#!/usr/bin/env python3
"""
Summary statistics for Rates and Descriptors from traderpros run.
"""

import argparse
import os
import pandas as pd
import numpy as np
from scipy import stats

def compute_mode(series):
    # For continuous data, use the mode from histogram
    hist, bin_edges = np.histogram(series, bins=50)
    mode_idx = np.argmax(hist)
    mode_val = (bin_edges[mode_idx] + bin_edges[mode_idx + 1]) / 2
    return mode_val

def main():
    parser = argparse.ArgumentParser(description="Compute summary statistics for Rates and Descriptors.")
    parser.add_argument('-i', '--input', required=True, help="Folder with results")
    parser.add_argument('-o', '--output', required=True, help="Output CSV file")
    parser.add_argument('-b', '--burn_in', type=float, default=0.10, help="Burn-in percentage")

    args = parser.parse_args()

    # Find files
    rates_file = None
    desc_file = None
    for file in os.listdir(args.input):
        if file.endswith('.Rates.log'):
            rates_file = os.path.join(args.input, file)
        elif file.endswith('.Descriptors.log'):
            desc_file = os.path.join(args.input, file)

    if not rates_file:
        raise FileNotFoundError("No .Rates.log file found")
    if not desc_file:
        raise FileNotFoundError("No .Descriptors.log file found")

    # Read data
    rates_df = pd.read_csv(rates_file, sep='\t')
    desc_df = pd.read_csv(desc_file, sep='\t')

    # Discard burn-in
    n_samples_rates = len(rates_df)
    burn_in_samples_rates = int(n_samples_rates * args.burn_in)
    rates_df = rates_df.iloc[burn_in_samples_rates:]

    n_samples_desc = len(desc_df)
    burn_in_samples_desc = int(n_samples_desc * args.burn_in)
    desc_df = desc_df.iloc[burn_in_samples_desc:]

    # Get columns, exclude Iteration, Posterior, Likelihood, Prior
    exclude_cols = ['Iteration', 'Posterior', 'Likelihood', 'Prior']
    rates_cols = [col for col in rates_df.columns if col not in exclude_cols]
    desc_cols = [col for col in desc_df.columns if col not in exclude_cols]

    # Unique variables, prefer rates if overlap
    all_cols = list(set(rates_cols + desc_cols))
    variables = []
    for col in all_cols:
        if col in rates_cols:
            variables.append((col, rates_df))
        elif col in desc_cols:
            variables.append((col, desc_df))

    # Compute stats
    stats_dict = {}
    for var, df in variables:
        series = df[var]
        mean_val = series.mean()
        median_val = series.median()
        mode_val = compute_mode(series)
        # MAP: value with max posterior
        max_post_idx = df['Posterior'].idxmax()
        map_val = df.loc[max_post_idx, var]
        q025 = series.quantile(0.025)
        q975 = series.quantile(0.975)
        min_val = series.min()
        max_val = series.max()
        stats_dict[var] = {
            'mean': mean_val,
            'median': median_val,
            'mode': mode_val,
            'MAP': map_val,
            '2.5%': q025,
            '97.5%': q975,
            'min': min_val,
            'max': max_val
        }

    # To DataFrame
    summary_df = pd.DataFrame.from_dict(stats_dict, orient='index')

    # Save
    summary_df.to_csv(args.output)
    print(f"Summary statistics saved to {args.output}")

if __name__ == '__main__':
    main()