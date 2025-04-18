#!/usr/bin/env python3

import argparse
import os
import subprocess
import re
import shutil
import glob

def create_revscript(args):
    revscript = f"""
    # Parameters
    tree_path = "{args.phylogeny}"
    transition_rates = [{args.transition_rates}]
    root_prior = simplex({args.root_probs})
    sp_completion_rates = [{args.sp_completion_rates}]
    outFolder = "{args.out_dir}/"
    outPrefix = "{args.out_prefix}"
    n_simulations = {args.n_simulations}
    seed_number = {args.seed_number}
    n_unknown_sp = {args.n_unknown_sp}

    # Begin RevScript for simulation
    if (seed_number!=0){{
      seed(seed_number)
    }}
    
    printSeed()
    
    timetree <- readTrees(tree_path)[1]

    for (n in 1:n_simulations){{
        print("Simulation", n, "of", n_simulations, separator=" ")

        Q := fnFreeK(transition_rates, rescale=FALSE)
        
        root_prior_observed_states := simplex(root_prior[1], root_prior[2])
        
        trait_evol ~ dnPhyloCTMCDASiteIID(timetree,
                                          Q,
                                          branchRates=1,
                                          type="Standard",
                                          nSites=1,
                                          rootFrequencies=root_prior_observed_states,
                                          treatAmbiguousAsGap=FALSE)
        
        n_nodes <- timetree.nnodes()
        n_branches <- n_nodes -1
        
        for(i in 1:n_branches) {{
            state_branch_rate[i] := sum(trait_evol.relativeTimeInStates(i,1) * sp_completion_rates)
        }}
        
        for (i in 1:(n_branches)){{
            branch_speciation_events[i] ~ dnPoisson(state_branch_rate[i]*timetree.branchLength(i))
        }}
        
        writeNexus(timetree, filename=outFolder + outPrefix + "_Tree.Sim" + n + ".tre")
        writeNexus(trait_evol, filename=outFolder + outPrefix + "_trait.Sim" + n + ".nexus")
        
        topology[n] <- readBranchLengthTrees(outFolder + outPrefix + "_Tree.Sim" + n + ".tre")[1]
        topology[n].renumberNodes(timetree)
        #topology[n] <- readBranchLengthTrees(tree_path)[1]
        SpeciationBranchTree[n] <- fnTreeAssembly(topology[n], branch_speciation_events)
        writeNexus(SpeciationBranchTree[n], filename=outFolder + outPrefix + "_TrBrLen.Sim" + n + ".tre")
        
        spEventsTree <- readBranchLengthTrees(outFolder + outPrefix + "_TrBrLen.Sim" + n + ".tre")[1]
        distanceMatrix <- fnTreePairwiseDistances(spEventsTree)
        
        num_tips <- spEventsTree.ntips()
        for (tip in 1:num_tips){{
            visited[tip] <- round(-1)
        }}
    
        # Loop through each tip and assign groups based on distance of 0
        for (i in 1:(num_tips)) {{
            # Only assign a new group if the tip has not been assigned yet
            if (visited[i] == -1) {{
                # Assign a new group
                group <- [distanceMatrix.names()[i]]
                visited[i] <- 1
                
                # Find all tips with distance 0 to the current tip and assign the same group
                for (j in 1:num_tips) {{
                    if (j > i && distanceMatrix.matrix()[i][j] == 0) {{
                        group.append(distanceMatrix.names()[j])
                        visited[j] <- 1
                    }}
                }}
    
                if (i==1) {{
                    groups <- [group]
                }}else{{
                    groups.append(group)
                }}
            }}
        }}

        # Create a file with the number of species
        n_species <- groups.size()
        write(n_species,
              filename=outFolder + outPrefix + "." + n_species + ".spp.Sim" + n + ".txt",
              append=FALSE)

        # Create a species matrix file 
        write("taxon",
              "species_code" + "\\n",
              filename=outFolder + outPrefix + ".SpeciesMatrix.Sim" + n + ".txt",
              append=FALSE)
        
        for (i in 1:groups.size()) {{
            for (j in 1:groups[i].size()){{
            write(groups[i][j],
                  i-1,
                  "\\n",
                  filename=outFolder + outPrefix + ".SpeciesMatrix.Sim" + n + ".txt",
                  append=TRUE)
            }}
        }}

        if(n_unknown_sp  > 0){{
          matrix = readDelimitedDataFile(outFolder + outPrefix + ".SpeciesMatrix.Sim" + n + ".txt")
          maxline = matrix.size()-n_unknown_sp
          initial_line ~ dnUniformNatural(2, round(maxline) )
          final_line <-  initial_line + n_unknown_sp -1
          lines <- initial_line:final_line
          write("",filename=outFolder + outPrefix + ".SpeciesMatrixUnkn.Sim" + n + ".txt" )
          for (i in 1:matrix.size()){{
            if(lines.contains(i)){{
              write(matrix[i][1],
                     "?" + "\\n",
                    filename=outFolder + outPrefix + ".SpeciesMatrixUnkn.Sim" + n + ".txt",
                    append=TRUE)
            }}else{{
              write(matrix[i][1],
                     matrix[i][2] + "\\n",
                    filename=outFolder + outPrefix + ".SpeciesMatrixUnkn.Sim" + n + ".txt",
                    append=TRUE)
           }}
          }}
        }}
    }}

    printSeed()
    quit()
    """
    return revscript.strip()


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

def correct_state_numbers(directory):
    """
    Iterates over files in the given directory with the pattern *_trait.Sim
    and replaces 'symbols="0123"' with 'symbols="01"', saving changes to the same file.
    """
    # Pattern to match files
    pattern = os.path.join(directory, '*_trait.Sim*')
    
    # Get all matching files
    files = glob.glob(pattern)
    
    for file_path in files:
        # Read file content
        with open(file_path, 'r') as file:
            content = file.read()
        # Replace the target string
        updated_content = content.replace('symbols="0123"', 'symbols="01"')
        # Write the updated content back to the same file
        with open(file_path, 'w') as file:
            file.write(updated_content)

def organize(directory):
    # Regular expression to capture ".SimX." patterns
    pattern = re.compile(r'\.Sim(\d+)\.')
    
    # Iterate over files in the specified directory
    for filename in os.listdir(directory):
        match = pattern.search(filename)
        if match:
            # Extract the folder name, e.g., "Sim1", "Sim2", etc.
            folder_name = f"Sim{match.group(1)}"
            folder_path = os.path.join(directory, folder_name)
            
            # Create the folder if it does not exist
            os.makedirs(folder_path, exist_ok=True)
            
            # Move the file into the folder
            source_path = os.path.join(directory, filename)
            destination_path = os.path.join(folder_path, filename)
            shutil.move(source_path, destination_path)
            print(f"Moved '{filename}' to '{folder_path}'")

def main():
    parser = argparse.ArgumentParser(description='Simulates data under the Traderpros model.')
    parser.add_argument('-phy', '--phylogeny', type=str, help='Phylogenetic tree in nexus or newick.')
    parser.add_argument('-s', '--seed_number', type=int, default=0, help='Seed number for replication (default: 0).')
    parser.add_argument('-opre', '--out_prefix', type=str, default='Traderpros', help='Output file prefix (default: "Traderpros").')
    parser.add_argument('-odir', '--out_dir', type=str, default='trader_out', help='Output directory (default: "trader_out").')
    parser.add_argument('-tr', '--transition_rates', type=str, required=True, help='A comma-separated list of transition rates between observed states. It should be absolute rates. First value represents transition from 0 to 1 and second value should be the 1 to 0 transition (mandatory).')
    parser.add_argument('-rp', '--root_probs', type=str, required=True, help='A comma-separated list of probabilities for the all 4 root states (Observed 0 Hidden A, Observed 1 Hidden A, Observed 0 Hidden B, Observed 1 Hidden B).  It could be values that sum to 1, or relative weights for a simplex. E.g. -rp 1,0,1,0 will translate to 0.5, 0.0, 0.5, 0.0 (mandatory).')
    parser.add_argument('-sc', '--sp_completion_rates', type=str, required=True, help='A comma-separated list of state specific speciation completion rates for observed states (mandatory).')
    parser.add_argument('-ns', '--n_simulations', type=float, default=1, help='The number of simulations (default: 1).')
    parser.add_argument('-justscript', '--just_script', action="store_true", help='Create script only, do not run (default: False).')
    parser.add_argument("-org", "--organize_folder", action="store_true", help='Organize files into folders by simulation. (default: False)')
    parser.add_argument("-unk", "--n_unknown_sp",  type=int, default=0, help='Number of tips to be set as "unknown" species assignment in the Species Matrix file for testing application of the model for species delimitation (default: 0)')

    args = parser.parse_args()

    # Convert out_dir to absolute path
    args.out_dir = os.path.abspath(args.out_dir)
    args.phylogeny = os.path.abspath(args.phylogeny)

    # Create the revscript with parameters
    revscript = create_revscript(args)
    
    # Check if the directory exists and confirm with the user
    if os.path.exists(args.out_dir):
        print(f"Warning: Directory '{args.out_dir}' already exists. Proceeding may erase previous simulations.")
        if not confirm_proceed():
            exit()

    # Create output directory
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # Write the revscript to the output file
    script_filename = os.path.join(args.out_dir, f"{args.out_prefix}.SCRIPT.Rev")
    with open(script_filename, 'w') as script_file:
        script_file.write(revscript)

    # Run the RevBayes script if not just creating the script
    if not args.just_script:
        subprocess.run(['rb', script_filename], check=True)
    
    # Correct the state numbers in the trait nexus file
    correct_state_numbers(args.out_dir)

    # Organize folder
    if args.organize_folder:
        organize(args.out_dir)
        print("Files organized.")

if __name__ == '__main__':
    main()
