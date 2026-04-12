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
    extinction_hidden = [{args.pop_formation_hidden}]
    global_trans_rate = {args.global_trans_rate}
    relative_transition = [{args.relative_transition}]
    transition_hidden = [{args.transition_hidden}]
    root_prior = simplex({args.root_prior})
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
    
    # Create global_pop_formation and global_extinction from global_net_pop_formation and global turnover
    denom := abs(1.0 - global_pop_turnover) 
    global_pop_formation := global_net_pop_formation / denom
    global_extinction := (global_pop_turnover * global_net_pop_formation) / denom

    # Create observed rates
    pop_formation_observed := pop_formation_relative * global_pop_formation
    extinction_observed := extinction_relative * global_extinction

    # Create the combined observed and hidden rates
    for (j in 1:pop_formation_hidden.size()) {{
        for (i in 1:pop_formation_relative.size()) {{
            index = i+(j*pop_formation_relative.size())-pop_formation_relative.size()
            pop_formation[index] := pop_formation_observed[i] * pop_formation_hidden[j]
            extinction[index] := extinction_observed[i] * extinction_hidden[j]
        }}
    }}
    
    # Q matrix
    ## Observed transition
    transition_rates := relative_transition * global_trans_rate
    
    ## Hidden transition rates
    for (i in 1:(pop_formation_hidden.size() * (pop_formation_hidden.size() - 1))) {{
        R[i] := transition_hidden
    }}

    # Rate matrix with hidden states
    rate_matrix := fnHiddenStateRateMatrix(transition_rates, R, rescaled=FALSE)

    # Speciation completion rates
    spCompletion_observed := spCompletion_relative * global_spCompletion
    
    for (j in 1:spCompletion_hidden.size()) {{
        for (i in 1:spCompletion_observed.size()) {{
            index = i+(j*spCompletion_observed.size())-spCompletion_observed.size()
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
          "hidden_transition," + R + "\\n",
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

    for (n in 1:N_SIMULATIONS){{
        
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
                                  rootFrequencies=root_prior,
                                  treatAmbiguousAsGap=TRUE)
        

        writeNexus(timetree, filename= OUT_DIR + OUT_PREFIX + "_Tree.Sim_" + n + ".tre")
        writeNexus(trait_evol, filename= OUT_DIR + OUT_PREFIX + "_trait.Sim_" + n + ".nexus")
        
        # Updated version of RevBayes allows for this method, which could be interesting to save trait histories
        #char_hist <- trait_evol.characterHistories() 
        #writeNexus(char_hist, filename= OUT_DIR + OUT_PREFIX + "_TraitHist.Sim_" + n + ".tre")
 
        n_nodes <- timetree.nnodes()
        n_branches <- n_nodes -1

        for(i in 1:n_branches) {{
            branch_spCompletion[i] := sum(trait_evol.relativeTimeInStates(i,1) * spCompletion)
        }}

        for (i in 1:(n_branches)){{
            branch_speciation_events[i] ~ dnPoisson(branch_spCompletion[i]*timetree.branchLength(i))
        }}
        
        topology[n] <- readBranchLengthTrees(OUT_DIR + OUT_PREFIX + "_Tree.Sim_" + n + ".tre")[1]
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
        for (i in 1:(num_tips-1)) {{
            # Only assign a new group if the tip has not been assigned yet
            if (visited[i] == -1) {{
                # Assign a new group
                group <- [distanceMatrix.names()[i]]
                visited[i] <- 1
                
                # Find all tips with distance 0 to the current tip and assign the same group
                for (j in (i+1):num_tips) {{
                    if (distanceMatrix.matrix()[i][j] == 0) {{
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
