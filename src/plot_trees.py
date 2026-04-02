#!/usr/bin/env python3
"""
Wrapper that runs an R script for Traderpros tree plots.
"""

import argparse
import os
import subprocess
import textwrap
import sys

def make_anc_state_labels(n_observed, n_hidden):
    hidden_letters = [chr(ord("A") + i) for i in range(n_hidden)]
    labels = []
    for h in hidden_letters:
        for o in range(n_observed):
            labels.append(f"Observed {o} Hidden {h}")
    return labels

def make_color_expression(colors):
    # Quote style for R
    return ", ".join(f'"{c}"' for c in colors)

def main():
    parser = argparse.ArgumentParser(description="Plot Traderpros tree outputs with RevGadgets")
    parser.add_argument("-i", "--input", required=True, help="Folder with Traderpros results")
    parser.add_argument("-o", "--output", required=True, help="Output folder name")
    parser.add_argument("-n", "--n_observed_states", type=int, default=2)
    parser.add_argument("-H", "--n_hidden_states", type=int, default=2)
    parser.add_argument("-asc", "--anc_state_colors", nargs="+", default=None)
    parser.add_argument("-rc", "--rate_colors", nargs=2, default=["#2166AC", "#B2182B"])
    parser.add_argument("-l", "--layout", default="rectangular",
                        choices=['rectangular', 'cladogram', 'slanted', 'ellipse',
                                 'roundrect', 'fan', 'circular', 'inward_circular',
                                 'radial', 'equal_angle', 'daylight', 'ape'])
    parser.add_argument("-b", "--burn_in", type=float, default=0.1)
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

    # R-compatible variables
    input_dir = os.path.abspath(args.input).replace("\\", "/")
    output_dir = os.path.abspath(args.output).replace("\\", "/")
    rc_colors = make_color_expression(args.rate_colors)

    # If no anc colors specified, will be NULL in R and code defaults
    if args.anc_state_colors is None:
        asc_expr = "NULL"
    else:
        asc_expr = "c(" + make_color_expression(args.anc_state_colors) + ")"

    # Prepare R script text
    r_script = f"""
    library(ggplot2)
    library(RevGadgets)
    library(coda)
    library(RColorBrewer)
    library(viridis)

    input <- "{input_dir}"
    output <- "{output_dir}"
    n_observed_states <- {args.n_observed_states}
    n_hidden_states <- {args.n_hidden_states}
    layout_type <- "{args.layout}"
    rate_colors <- c({rc_colors})
    anc_state_colors <- {asc_expr}
    burn <- {args.burn_in}

    # files
    spBrRateTreeFile <- list.files(path = input, pattern = ".Protracted.MAP.tre", full.names = TRUE)[1]
    ancStatesTreeFile_cond <- list.files(path = input, pattern = ".traits.MAP.cond.tree", full.names = TRUE)[1]
    ancStatesTreeFile_marg <- list.files(path = input, pattern = ".traits.MAP.marg.tree", full.names = TRUE)[1]
    branchRatesFile <- list.files(path = input, pattern = ".BirthDeathBrRates.log", full.names = TRUE)[1]

    # dynamic labels
    hidden_letters <- toupper(letters[1:n_hidden_states])
    states <- expand.grid(obs = 0:(n_observed_states - 1), hid = hidden_letters)
    anc_state_labels <- paste0(states$obs, states$hid)

    if (is.null(anc_state_colors)) {{
      if (length(anc_state_labels) == 4) {{
        anc_state_colors <- c("#0072B2", "#E69F00", "#009E73", "#CC3311")
      }} else {{
        anc_state_colors <- viridis::turbo(length(anc_state_labels))
      }}
    }}

    if (!dir.exists(output)) dir.create(output, recursive = TRUE)

    # protracted rate tree
    print("Plotting tree with speciation completion rates")
    spBrRateTree <- readTrees(spBrRateTreeFile)
    SpCompletionRate <- plotTree(tree = spBrRateTree,
                                 node_age_bars = FALSE,
                                 node_pp = FALSE,
                                 tip_labels = TRUE,
                                 tip_labels_size = 3,
                                 color_branch_by = "branch_spCompletion",
                                 branch_color = rate_colors,
                                 line_width = 2.5) +
      ggplot2::theme(legend.position.inside = c(.1, .9))
    ggsave(SpCompletionRate, file = file.path(output, paste0("Traderpros.SpCompletionRate.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
    ggsave(SpCompletionRate, file = file.path(output, paste0("Traderpros.SpCompletionRate.MAP.tre.tiff")), width = 18, height = 24, device = "tiff")

    print("Plotting tree with speciation completion probabilities")
    SpCompletionProbs <- plotTree(tree = spBrRateTree,
                                 node_age_bars = FALSE,
                                 node_pp = FALSE,
                                 tip_labels = TRUE,
                                 tip_labels_size = 3,
                                 color_branch_by = "p_speciation",
                                 branch_color = rate_colors,
                                 line_width = 2.5) +
      ggplot2::theme(legend.position.inside = c(.1, .9))
    ggsave(SpCompletionProbs, file = file.path(output, paste0("Traderpros.SpeciationProb.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
    ggsave(SpCompletionProbs, file = file.path(output, paste0("Traderpros.SpeciationProb.MAP.tre.tiff")), width = 18, height = 24, device = "tiff")

    
    print("Plotting tree with number of speciation events")
    spEventsTreeFile <- list.files(path = input, pattern = ".SpEvents.MAP.tre", full.names = TRUE)[1]
    spEventsTree <- readTrees(spEventsTreeFile)
    sp_eventsTree <- plotTree(tree = spEventsTree,
                              node_age_bars = FALSE,
                              node_pp = FALSE,
                              tip_labels = TRUE,
                              tip_labels_size = 3,
                              color_branch_by = "branch_speciation_events",
                              branch_color = rate_colors,
                              line_width = 2.5) +
      ggplot2::theme(legend.position.inside = c(.1, .9))
    ggsave(sp_eventsTree, file = file.path(output, "Traderpros.SpEvents.tree.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(sp_eventsTree, file = file.path(output, "Traderpros.SpEvents.tree.tiff"), width = 18, height = 24, device = "tiff")

    # Create named state_labels vector for processAncStates
    state_names <- as.character(0:(length(anc_state_labels)-1))
    named_state_labels <- stats::setNames(anc_state_labels, state_names)

    # Ensure anc_state_colors length matches number of states
    if (length(anc_state_colors) < length(anc_state_labels)) {{
      anc_state_colors <- viridis::turbo(length(anc_state_labels))
    }}

    print("Processing ancestral states")
    ancStatesTree_cond <- processAncStates(ancStatesTreeFile_cond,
                                           state_labels = named_state_labels)

    print("Plotting tree with ancestral states")
    ancStatesPlot <- plotAncStatesMAP(t = ancStatesTree_cond,
                                      cladogenetic = FALSE,
                                      node_size = c(1,5),
                                      node_color = anc_state_colors,
                                      node_color_as = "state",
                                      node_size_as = "state_posterior",
                                      tip_states_size = 2.5,
                                      tip_labels_size = 5,
                                      timeline = TRUE,
                                      geo = FALSE) +
      ggplot2::theme(legend.position = "left", legend.title = element_blank())

    ggsave(ancStatesPlot, file = file.path(output, "Traderpros.ASR.cond.MAP.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(ancStatesPlot, file = file.path(output, "Traderpros.ASR.cond.MAP.tiff"), width = 18, height = 24, device = "tiff")
    
    # Plot tree with BD Branch rates for birth and death
    print("Reading and saving Birth and Death tree with branch rates")
    branch_data <- readTrace(branchRatesFile)
    branchTree <- processBranchData(spBrRateTree, branch_data, burnin = burn, parnames = c("avg_lambda", "avg_mu", "num_shifts"), summary = "median", net_div = TRUE)
  
    birth <- plotTree(tree = branchTree,
                      node_age_bars = FALSE,
                      node_pp = FALSE,
                      tip_labels = TRUE,
                      tip_labels_size = 5,
                      color_branch_by = "avg_lambda",
                      branch_color = rate_colors,
                      line_width = 2.5) +
      ggplot2::theme(legend.position.inside = c(.1, .9))
  
    #plot(birth)
  
      ggsave(birth, file = file.path(output, "Traderpros.PopFormationTree.pdf"), width = 18, height = 24, device = "pdf")
      ggsave(birth, file = file.path(output, "Traderpros.PopFormationTree.tiff"), width = 18, height = 24, device = "png")
  
  
    ## Plot tree with death rates
    extinction <- plotTree(tree = branchTree, 
                           node_age_bars = FALSE,
                           node_pp = FALSE, 
                           tip_labels = TRUE,
                           tip_labels_size = 5,
                           color_branch_by = "avg_mu",
                           branch_color=rate_colors,
                           line_width = 2.5) + 
      ggplot2::theme(legend.position.inside = c(.1, .9))
  
      #plot(extinction)
  
      ggsave(extinction, file = file.path(output, "Traderpros.ExtinctionTree.pdf"), width = 18, height = 24, device = "pdf")
      ggsave(extinction, file = file.path(output, "Traderpros.ExtinctionTree.tiff"), width = 18, height = 24, device = "png")
  
    # Plot tree with net diversification rates 
    # (if you are using a population tree this is the net population fragmentation)
    net <- plotTree(tree = branchTree, 
                    node_age_bars = FALSE,
                    node_pp = FALSE, 
                    tip_labels = TRUE,
                    tip_labels_size = 5,
                    color_branch_by = "net_div",
                    branch_color= rate_colors,
                    line_width = 2.5) + 
      ggplot2::theme(legend.position.inside = c(.1, .9))
  
      #plot(net)
  
      ggsave(net, file = file.path(output, "Traderpros.netPopFormation.pdf"), width = 18, height = 24, device = "pdf")
      ggsave(net, file = file.path(output, "Traderpros.netPopFormation.tiff"), width = 18, height = 24, device = "png")

      print("Done")
    """
    # save temporary R script
    script_path = os.path.join(args.output, "plot_trees_generated.R")
    with open(script_path, "w") as f:
        f.write(textwrap.dedent(r_script))

    # execute Rscript
    try:
        subprocess.run(["Rscript", script_path], check=True)
    except subprocess.CalledProcessError as e:
        print("R script execution failed:", e)
        sys.exit(1)

if __name__ == "__main__":
    main()