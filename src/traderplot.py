#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Run all TraderPros plotting scripts")
    parser.add_argument('-i', '--input', required=True, help="Folder with Traderpros results")
    parser.add_argument('-o', '--output', required=True, help="Output folder name")
    parser.add_argument('-n', '--n_observed_states', type=int, default=2, help="Number of observed states")
    parser.add_argument('-H', '--n_hidden_states', type=int, default=2, help="Number of hidden states")
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


    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)

    # Commands to run
    commands = [
        ["plot_rates", "-i", str(args.input), "-o", str(args.output), "-n", str(args.n_observed_states), "-H", str(args.n_hidden_states), "--outliers", "exclude", "-f", "yes"],
        ["plot_rates_observed", "-i", str(args.input), "-o", str(args.output), "-n", str(args.n_observed_states), "--outliers", "exclude", "-f", "yes"],
        ["plot_descriptors", "-i", str(args.input), "-o", str(args.output), "-n", str(args.n_observed_states), "-H", str(args.n_hidden_states), "--outliers", "exclude", "-f", "yes"],
        ["plot_descriptors_obs", "-i", str(args.input), "-o", str(args.output), "-n", str(args.n_observed_states), "--outliers", "exclude", "-f", "yes"],
        ["plot_rj", "-i", str(args.input), "-o", str(args.output), "-f", "yes"],
        ["plot_trees", "-i", str(args.input), "-o", str(args.output), "-f", "yes"],
    ]
    
    for cmd in commands:
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, cwd=os.path.dirname(__file__))
        if result.returncode != 0:
            print(f"Command failed: {' '.join(cmd)}")
            sys.exit(1)
    
    print("All plotting scripts completed successfully.")

if __name__ == "__main__":
    main()