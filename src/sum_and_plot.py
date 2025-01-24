#!/usr/bin/env python3

import argparse
import subprocess
import os

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

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Python wrapper for sum_and_plot_traderpros_result in R.")
    parser.add_argument('-od', '--out_dir', type=str, default="./outputstrader", help='Output directory')
    parser.add_argument('-in', '--in_dir', type=str, required=True, help='Input directory')
    parser.add_argument('-pre', '--prefix_used_in_traderpros', type=str, required=True, help='Prefix used in TraderPros output')
    parser.add_argument('-b', '--burn', type=float, default=0.1, help='Burn-in percentage')
    parser.add_argument('--color_completion', type=str, default="#A6CEE3,#1F78B4", help='Comma-separated colors for completion rate')
    parser.add_argument('--color_n_speciation', type=str, default="#BDBDBD,#000000", help='Comma-separated colors for number of speciation events')
    parser.add_argument('--color_birth', type=str, default="#DADAEB,#3F007D", help='Comma-separated colors for birth rates')
    parser.add_argument('--color_death', type=str, default="#FB9A99,#E31A1C", help='Comma-separated colors for death rates')
    parser.add_argument('--color_net', type=str, default="#FFFF99,#B15928", help='Comma-separated colors for net diversification')
    parser.add_argument('--color_transition', type=str, default="#FDBF6F,#FF7F00", help='Comma-separated colors for transitions')
    parser.add_argument('--anc_state_labels', type=str, default="State 0 Hidden A,State 1 Hidden A,State 0 Hidden B,State 1 Hidden B", help='Comma-separated ancestral state labels')
    parser.add_argument('--anc_state_colors', type=str, default="#D73027,#A50026,#4575B4,#313695", help='Comma-separated colors for ancestral states')
    parser.add_argument('--path_to_custom_functions', type=str, default="./", help='Path to custom functions')
    parser.add_argument('--out_images_format', type=str, choices=['pdf', 'png', 'both'], default='both', help='Output image format')

    args = parser.parse_args()

    # Convert out_dir to absolute path
    args.out_dir = os.path.abspath(args.out_dir)
    args.in_dir   = os.path.abspath(args.in_dir)
    args.path_to_custom_functions  = os.path.abspath(args.path_to_custom_functions)
 
    # Check if the directory exists and confirm with the user
    if os.path.exists(args.out_dir):
        print(f"Warning: Directory '{args.out_dir}' already exists. Proceeding may erase previous simulations.")
        if not confirm_proceed():
            exit()

    # Format the arguments for R
    #anc_state_labels_r = ','.join([f'"{k}={v}"' for k, v in eval(f"{{{args.anc_state_labels}}}").items()])
    
    r_args = [
        args.out_dir,
        args.in_dir,
        args.prefix_used_in_traderpros,
        str(args.burn),
        args.color_completion,
        args.color_n_speciation,
        args.color_birth,
        args.color_death,
        args.color_net,
        args.color_transition,
        args.anc_state_labels,
        args.anc_state_colors,
        args.path_to_custom_functions,
        args.out_images_format
    ]

    # Call the R script using subprocess
    try:
        subprocess.run(["Rscript", os.path.join(args.path_to_custom_functions, "sum_and_plot_traderpros_result.R")] + r_args, check=True)
    except subprocess.CalledProcessError as e:
        print("Error in R script execution:", e)

if __name__ == "__main__":
    main()
