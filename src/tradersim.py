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
    OUT_PREFIX   = "{args.out_prefix}"
    OUT_DIR      = "{args.out_dir}" + "/"
    N_SIMULATIONS = {args.n_simulations}
    SEED_NUMBER = {args.seed_number}

    global_net_pop_formation = {args.global_net_pop_formation}
    global_pop_turnover = {args.global_pop_turnover}
    pop_formation_relative = [{args.pop_formation_relative}]
    extinction_relative = [{args.extinction_relative}]
    pop_formation_hidden = [{args.pop_formation_hidden}]
    extinction_hidden = [{args.extinction_hidden}]
    global_trans_rate = {args.global_trans_rate}
    relative_transition = [{args.relative_transition}]
    transition_hidden = {args.transition_hidden}
    root_prior = [{args.root_probs}]
    rho = {args.rho}
    root = {args.root_age}
    N_TIPS = {args.n_tips}
    global_spCompletion = {args.global_spCompletion}
    spCompletion_relative = [{args.spCompletion_relative}]
    spCompletion_hidden = [{args.spCompletion_hidden}]
    N_UNKNOWN = {args.n_unknown_sp}


    # Begin RevScript for simulation
    if (SEED_NUMBER!=0){{
      seed(SEED_NUMBER)
    }}
    
    printSeed()

    NUM_HIDDEN = pop_formation_hidden.size()
    num_states = pop_formation_relative.size()

    
    # Create global_pop_formation and global_extinction from global_net_pop_formation and global turnover
    denom := abs(1.0 - global_pop_turnover) 
    global_pop_formation := global_net_pop_formation / denom
    global_extinction := (global_pop_turnover * global_net_pop_formation) / denom

    # Create observed rates
    pop_formation_observed := pop_formation_relative * global_pop_formation
    extinction_observed := extinction_relative * global_extinction

    # Create the combined observed and hidden rates
    for (j in 1:NUM_HIDDEN) {{
        for (i in 1:num_states) {{
            index = i+(j*num_states)-num_states
            pop_formation[index] := pop_formation_observed[i] * pop_formation_hidden[j]
            extinction[index] := extinction_observed[i] * extinction_hidden[j]
        }}
    }}
    
    # Q matrix
    ## Observed transition
    transition_rates := relative_transition * global_trans_rate
    
    ## Hidden transition rates
    for (i in 1:(NUM_HIDDEN * (NUM_HIDDEN - 1))) {{
        R[i] := transition_hidden
    }}
    
    # Rate matrix with hidden states
    rate_matrix := fnHiddenStateRateMatrix(transition_rates, R, rescaled=FALSE)

    # Speciation completion rates
    spCompletion_observed := spCompletion_relative * global_spCompletion
    
    for (j in 1:NUM_HIDDEN) {{
        for (i in 1:num_states) {{
            index = i+(j*num_states)-num_states
            spCompletion[index] := spCompletion_observed[i] *spCompletion_hidden[j]
        }}
    }}

    # Write a file with the parameter values              
    write("spCompletion," + spCompletion + "\\n",
          "spCompletion_observed," + spCompletion_observed + "\\n",
          "spCompletion_relative," + spCompletion_relative + "\\n",
          "global_spCompletion," + global_spCompletion + "\\n",
          "spCompletion_hidden," + spCompletion_hidden + "\\n",
          "transition_rates," + transition_rates + "\\n",
          "hidden_transition," + transition_hidden + "\\n",
          "relative_transition," + relative_transition + "\\n",
          "global_trans_rate," + global_trans_rate + "\\n",
          "pop_formation," + pop_formation + "\\n",
          "pop_formation_observed," + pop_formation_observed + "\\n",
          "pop_formation_hidden," + pop_formation_hidden + "\\n",
          "global_pop_formation," + global_pop_formation + "\\n",
          "pop_formation_relative," + pop_formation_relative + "\\n",
          "extinction," + extinction + "\\n",
          "extinction_observed," + extinction_observed + "\\n",
          "extinction_hidden," + extinction_hidden + "\\n",
          "extinction_relative," + extinction_relative + "\\n",
          "global_extinction," + global_extinction + "\\n",
          "root_prior," + root_prior + "\\n",
          "global_pop_turnover," + global_pop_turnover + "\\n",
          "global_net_pop_formation," + global_net_pop_formation + "\\n",
          "N_SIMULATIONS," + N_SIMULATIONS + "\\n",
          "OUT_PREFIX," + OUT_PREFIX + "\\n",
           filename= OUT_DIR + OUT_PREFIX + ".Param.Sim" + ".csv",
           append=FALSE, separator = "")

    timetree ~ dnCDBDP( 
        rootAge           = root,
        speciationRates   = pop_formation,
        extinctionRates   = extinction,
        Q                 = rate_matrix,
        pi                = root_prior,
        rho               = rho,
        simulateCondition = "numTips",
        exactNumLineages  = N_TIPS
        )

    trait_evol ~ dnPhyloCTMCDASiteIID(tree=timetree,
                              Q=rate_matrix,
                              branchRates=1,
                              type="NaturalNumbers",
                              nSites=1,
                              rootFrequencies=root_prior)

    trait_evol.setValue(timetree.getCharData())

    for (n in 1:N_SIMULATIONS){{
        print("Simulation:",n)
        timetree.redraw()
        trait_evol.redraw()
        trait_evol.setValue(timetree.getCharData())

        max_tries <- 100
        tries <- 0
        invariable <- trait_evol.numInvariableBlocks()
        while (invariable != 0 & tries < max_tries ){{
            tries <- tries + 1
            print("Redrawing trait evolution on tree because no trait change was simulated")
            timetree.redraw()
            trait_evol.redraw()
            invariable <- trait_evol.numInvariableBlocks()
        }}

        if (invariable != 0 ){{
          print("No variable trait simulated after ", tries, " tries. ","Try different simulation conditions. Likely trait transitions are too low for the root age and branch lengths", separator="")
          quit()
        }}
 
        writeNexus(timetree, filename= OUT_DIR + OUT_PREFIX + "_Tree.Sim" + n + ".tre")
        writeNexus(trait_evol, filename= OUT_DIR + OUT_PREFIX + "_trait.Sim" + n + ".nexus")
        
        # Write file with state frequency
        state_freq <- timetree.getCharData().getEmpiricalBaseFrequencies()
        write(state_freq, filename= OUT_DIR + OUT_PREFIX + ".StateFreq.Sim" + n + ".txt")

        # Updated version of RevBayes allows for this method, which could be interesting to save trait histories
        #char_hist <- trait_evol.characterHistories() 
        #writeNexus(char_hist, filename= OUT_DIR + OUT_PREFIX + "_TraitHist.Sim" + n + ".tre")
 
        n_nodes <- timetree.nnodes()
        n_branches <- n_nodes -1

        for(i in 1:n_branches) {{
            branch_spCompletion[i] := sum(trait_evol.relativeTimeInStates(i,1) * spCompletion)
        }}

        for (i in 1:(n_branches)){{
            branch_speciation_events[i] ~ dnPoisson(branch_spCompletion[i]*timetree.branchLength(i))
        }}
        
        topology[n] <- readBranchLengthTrees(OUT_DIR + OUT_PREFIX + "_Tree.Sim" + n + ".tre")[1]
        topology[n].renumberNodes(timetree)
    
        SpeciationBranchTree[n] <- fnTreeAssembly(topology[n], branch_speciation_events)
        writeNexus(SpeciationBranchTree[n], filename=OUT_DIR + OUT_PREFIX + "_TrBrLen.Sim" + n + ".tre")
        
        spEventsTree <- readBranchLengthTrees(OUT_DIR + OUT_PREFIX + "_TrBrLen.Sim" + n + ".tre")[1]
        distanceMatrix <- fnTreePairwiseDistances(spEventsTree)
      
        num_tips <- spEventsTree.ntips()
        for (tip in 1:num_tips){{
            visited[tip] <- round(-1)
        }}
    
        # Loop through each tip and assign groups based on distance of 0
        for (i in 1:num_tips) {{
            # Only assign a new group if the tip has not been assigned yet
            if (visited[i] == -1) {{
                # Assign a new group
                group <- [distanceMatrix.names()[i]]
                visited[i] <- 1
                
                # Find all tips with distance 0 to the current tip and assign the same group
                if (i < num_tips) {{
                    for (j in (i+1):num_tips) {{
                        if (distanceMatrix.matrix()[i][j] == 0) {{
                            group.append(distanceMatrix.names()[j])
                            visited[j] <- 1
                        }}
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
        write("Number of species: ", n_species,
              filename= OUT_DIR + OUT_PREFIX + "." + n_species + ".spp.Sim" + n + ".txt",
              append=FALSE)

        # Create a species matrix file               
        write("taxon",
              "species_code" + "\\n",
              filename= OUT_DIR + OUT_PREFIX + ".SpeciesMatrix.Sim" + n + ".txt",
              append=FALSE)
        
        for (i in 1:groups.size()) {{
            for (j in 1:groups[i].size()){{
            write(groups[i][j],
                  i,
                  "\\n",
                  filename= OUT_DIR + OUT_PREFIX + ".SpeciesMatrix.Sim" + n + ".txt",
                  append=TRUE)
            }}
        }}

        if(N_UNKNOWN  > 0){{
          matrix = readDelimitedDataFile(OUT_DIR + OUT_PREFIX + ".SpeciesMatrix.Sim" + n + ".txt")
          maxline = matrix.size()-N_UNKNOWN
          initial_line ~ dnUniformNatural(2, round(maxline) )
          final_line <-  initial_line + N_UNKNOWN -1
          lines <- initial_line:final_line
          write("",filename= OUT_DIR + OUT_PREFIX + ".SpeciesMatrixUnkn.Sim" + n + ".txt" )
          for (i in 1:matrix.size()){{
            if(lines.contains(i)){{
              write(matrix[i][1],
                     "?" + "\\n",
                    filename= OUT_DIR + OUT_PREFIX + ".SpeciesMatrixUnkn.Sim" + n + ".txt",
                    append=TRUE)
            }}else{{
              write(matrix[i][1],
                     matrix[i][2] + "\\n",
                    filename= OUT_DIR + OUT_PREFIX + ".SpeciesMatrixUnkn.Sim" + n + ".txt",
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

def observed_state_numbers(directory, pop_formation_relative, pop_formation_hidden):
    """
    Create a new trait nexus file with only observed states for each *_trait.Sim*.nexus file.
    Maps combined hidden/observed states to the observed state index and writes a new file
    named with the suffix _trait.OBS.Sim.
    """
    def parse_values(value_string):
        return [v.strip() for v in value_string.split(',') if v.strip()]

    observed_states = parse_values(pop_formation_relative)
    hidden_states = parse_values(pop_formation_hidden)
    n_obs = len(observed_states)
    n_hidden = len(hidden_states)
    if n_obs < 1 or n_hidden < 1:
        return

    symbol_string = ''.join(str(i) for i in range(n_obs))
    pattern = os.path.join(directory, '*_trait.Sim*.nexus')
    matrix_start_re = re.compile(r'^\s*matrix\b', re.IGNORECASE)
    symbols_re = re.compile(r'symbols="[^"]*"')

    for file_path in glob.glob(pattern):
        with open(file_path, 'r') as file:
            lines = file.readlines()

        output_lines = []
        in_matrix = False
        for line in lines:
            if not in_matrix:
                if matrix_start_re.match(line):
                    in_matrix = True
                    output_lines.append(line)
                    continue

                replaced_line = symbols_re.sub(f'symbols="{symbol_string}"', line)
                output_lines.append(replaced_line)
                continue

            # Still in matrix block
            if ';' in line:
                in_matrix = False
                output_lines.append(line)
                continue

            stripped = line.strip()
            if not stripped or stripped.startswith('['):
                output_lines.append(line)
                continue

            matrix_match = re.match(r'^(\s*)(\S+)(\s+)(\S.*)$', line)
            if not matrix_match:
                output_lines.append(line)
                continue

            prefix, taxon, sep, states = matrix_match.groups()
            mapped_states = []
            for ch in states:
                if ch.isdigit():
                    state_index = int(ch)
                    if state_index < n_obs * n_hidden:
                        mapped_states.append(str(state_index % n_obs))
                    else:
                        mapped_states.append(ch)
                else:
                    mapped_states.append(ch)

            output_lines.append(f"{prefix}{taxon}{sep}{''.join(mapped_states)}\n")

        output_file = os.path.join(
            directory,
            os.path.basename(file_path).replace('_trait.Sim', '_trait.OBS.Sim', 1)
        )
        with open(output_file, 'w') as out_file:
            out_file.writelines(output_lines)


def replace_sp_by_pop(directory):
    """
    Replace occurrences of 'sp' followed by digits with 'pop' in all files under directory.
    This prevents changing the 'sp' in words such as 'species'.
    """
    regex = re.compile(r'\bsp(?=\d)')
    for root, _, files in os.walk(directory):
        for filename in files:
            file_path = os.path.join(root, filename)
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
                content = file.read()
            updated_content = regex.sub('pop', content)
            if updated_content != content:
                with open(file_path, 'w', encoding='utf-8') as file:
                    file.write(updated_content)

def correct_data_type(directory):
    """
    Iterates over files in the given directory with the pattern *_trait.Sim
    and replaces 'NaturalNumbers' with 'Standard', saving changes to the same file.
    """
    # Pattern to match files
    pattern = os.path.join(directory, '*_trait*')
    
    # Get all matching files
    files = glob.glob(pattern)
    
    for file_path in files:
        # Read file content
        with open(file_path, 'r') as file:
            content = file.read()
        # Replace the target string
        updated_content = content.replace('NaturalNumbers', 'Standard')
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
    parser.add_argument('-s', '--seed_number', type=int, default=0, help='Random seed for replication; use 0 to randomize each run.')
    parser.add_argument('-pre', '--out_prefix', type=str, default='Traderpros', help='Prefix for all output files.')
    parser.add_argument('-o', '--out_dir', type=str, default='trader_out', help='Directory where output files will be written.')
    parser.add_argument('-fg', '--global_net_pop_formation', type=str, required=True, help='Global net population formation rate (birth - death).')
    parser.add_argument('-tg', '--global_pop_turnover', type=str, required=True, help='Global turnover rate (extinction/pop formation ratio).')
    parser.add_argument('-fr', '--pop_formation_relative', type=str, required=True, help='Comma-separated observed-state proportions for population formation; values should sum to 1.')
    parser.add_argument('-er', '--extinction_relative', type=str, required=True, help='Comma-separated observed-state proportions for extinction; values should sum to 1. Numbers of values must match number of pop_formation_relative values.')
    parser.add_argument('-fh', '--pop_formation_hidden', type=str, required=True, help='Comma-separated hidden-state modifiers for population formation; one value per hidden state.')
    parser.add_argument('-eh', '--extinction_hidden', type=str, required=True, help='Comma-separated hidden-state modifiers for extinction; one value per hidden state. Numbers of values must match number of pop_formation_hidden values.')
    parser.add_argument('-gt', '--global_trans_rate', type=float, required=True, help='Global transition rate for the observed trait.')
    parser.add_argument('-tr', '--relative_transition', type=str, required=True, help='Comma-separated relative transition weights for observed-state trait changes. Number of values must be n_observed_states * (n_observed_states - 1).')
    parser.add_argument('-th', '--transition_hidden', type=float, required=True, help='Transition rates between hidden states. Only one absolute value since we assumed equal rates for hidden states. Note: very high hidden transition relative to the observed transition and to the speciation completion might generate errors.')
    parser.add_argument('-rp', '--root_probs', type=str, required=True, help='Comma-separated root prior probabilities for all state combinations (i.e.number of values should be number of hidden * number of observed states).')
    parser.add_argument('-sf', '--rho', type=float, required=True, help='Sampling probability (rho) for observed tips.')
    parser.add_argument('-ra', '--root_age', type=float, default=10, help='(Not in use) Root age for the simulation; retained for RevBayes script flexibility only. Not used in the CLI here because the simulation is conditioned on the number of tips.')
    parser.add_argument('-nt', '--n_tips', type=int, required=True, help='Number of tips in the simulated tree.')
    parser.add_argument('-sg', '--global_spCompletion', type=float, required=True, help='Global speciation completion rate.')
    parser.add_argument('-sr', '--spCompletion_relative', type=str, required=True, help='Comma-separated observed-state proportional weights for speciation completion rates. Numbers of values must match number of pop_formation_relative values.')
    parser.add_argument('-sh', '--spCompletion_hidden', type=str, required=True, help='Comma-separated hidden-state modifiers for speciation completion rates. Numbers of values must match number of pop_formation_hidden values.')
    parser.add_argument('-ns', '--n_simulations', type=float, default=1, help='Number of replicate simulations to run.')
    parser.add_argument('-justscript', '--just_script', action="store_true", help='Write the RevBayes script only and do not execute it.')
    parser.add_argument("-org", "--organize_folder", action="store_true", help='Organize output files into per-simulation folders.')
    parser.add_argument("-unk", "--n_unknown_sp",  type=int, default=0, help='Number of tips to mark as unknown in the species matrix for species delimitation testing.')

    args = parser.parse_args()

    # Convert out_dir to absolute path
    args.out_dir = os.path.abspath(args.out_dir)

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
    
    # Create observed-state-only trait nexus files
    observed_state_numbers(args.out_dir, args.pop_formation_relative, args.pop_formation_hidden)

    # Replace sp code with pop code in output files before organizing
    replace_sp_by_pop(args.out_dir)
    
    # Correct data type
    correct_data_type(args.out_dir)

    # Organize folder
    if args.organize_folder:
        organize(args.out_dir)
        print("Files organized.")

if __name__ == '__main__':
    main()
