#' Summarize and Plot Results from TraderPros Output
#'
#' @description This function processes TraderPros output data by summarizing, plotting, 
#' and saving several key aspects of the analysis results, including state-dependent 
#' probabilities, trace summaries, and ancestral state reconstructions.
#'
#' @param out_dir Output directory to save results. Default is `./outputstrader`.
#' @param in_dir Directory containing TraderPros outputs to be processed.
#' @param prefix_used_in_traderpros Prefix used for the TraderPros run output files.
#' @param burn Burn-in percentage for processing the trace file. Default is 0.1.
#' @param color_completion Vector of colors (length 2) for speciation completion rate. Default is `c("#A6CEE3", "#1F78B4")`.
#' @param color_n_speciation Vector of colors (length 2) for the number of speciation events. Default is `c("#BDBDBD", "#000000")`.
#' @param color_birth Vector of colors (length 2) for the birth rate. Default is `c("#DADAEB", "#3F007D")`.
#' @param color_death Vector of colors (length 2) for death rate. Default is `c("#FB9A99", "#E31A1C")`.
#' @param color_net Vector of colors (length 2) for net fragmentation rate. Default is `c("#FFFF99", "#B15928")`.
#' @param color_transition Vector of colors (length 2) for observed state transition rates. Default is `c("#FDBF6F","#FF7F00")`.
#' @param anc_state_labels Named vector (length 4) of labels for the ancestral states. Default labels provided.
#' @param anc_state_colors Vector (length 4) of colors for the ancestral states. Default is `c("#D73027","#A50026","#4575B4", "#313695")`.
#' @param path_to_custom_functions Path to directory containing custom functions (`custom_plot_trace.R`, `summary_stats_from_trace.R`, `plot_rj.R`). Default is `./`.
#' @param out_images_format Output image format, either `"pdf"`, `"png"`, or `"both"`.
#'
#' @return A list containing ggplot objects for each type of plot created.
#' @export
#'
#' @examples
#' sum_and_plot_traderpros_result(
#'   out_dir = "./outputstrader",
#'   in_dir = "./traderpros_data",
#'   prefix_used_in_traderpros = "run1",
#'   burn = 0.1,
#'   out_images_format = "both"
#' )
#'

sum_and_plot_traderpros_result <- function(out_dir = "./outputstrader",
                                           in_dir,
                                           prefix_used_in_traderpros,
                                           burn = 0.1,
                                           color_completion = c("#2166AC", "#B2182B"),
                                           color_n_speciation = c("#2166AC", "#B2182B"),
                                           color_birth = c("#2166AC", "#B2182B"),
                                           color_death = c("#2166AC", "#B2182B"),
                                           color_net = c("#2166AC", "#B2182B"),
                                           color_transition = c("#E08214", "#2D004B"),
                                           anc_state_labels = c("State 0 Hidden A",
                                                                "State 1 Hidden A",
                                                                "State 0 Hidden B",
                                                                "State 1 Hidden B"),
                                           anc_state_colors = c("#542788", "#2D004B", "#FDB863", "#E08214"),
                                           path_to_custom_functions = "./",
                                           out_images_format = "both") {

  # Load required libraries
  library(ggplot2)
  library(RevGadgets)
  library(coda)
  library(RColorBrewer)
  library(viridis)

  # Load custom functions
  source(file.path(path_to_custom_functions, "custom_plot_trace.R"))
  source(file.path(path_to_custom_functions, "summary_stats_from_trace.R"))
  source(file.path(path_to_custom_functions, "plot_rj.R"))

  # Ensure output directory exists
  if (!dir.exists(out_dir)) dir.create(out_dir)

  # File paths for TraderPros outputs
  spBrRateTreeFile <- file.path(in_dir, paste0(prefix_used_in_traderpros, ".Protracted.MAP.tre"))
  spEventsTreeFile <- file.path(in_dir, paste0(prefix_used_in_traderpros, ".SpEvents.MAP.tre"))
  ancStatesTreeFile_cond <- file.path(in_dir, paste0(prefix_used_in_traderpros, ".traits.MAP.cond.tree"))
  traceFile <- file.path(in_dir, paste0(prefix_used_in_traderpros, ".model.log"))
  branchRatesFile <- file.path(in_dir, paste0(prefix_used_in_traderpros, ".BirthDeathBrRates.log"))

  # Plot tree with speciation rates
  print("Reading and saving Protracted tree with annotations")
  spBrRateTree <- readTrees(spBrRateTreeFile)
  SpCompletionRate <- plotTree(tree = spBrRateTree,
                               node_age_bars = FALSE,
                               node_pp = FALSE,
                               tip_labels = TRUE,
                               tip_labels_size = 3,
                               color_branch_by = "state_branch_rate",
                               branch_color = color_completion,
                               line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(SpCompletionRate)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(SpCompletionRate, file = file.path(out_dir, paste0(prefix_used_in_traderpros, ".SpBrachRate.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(SpCompletionRate, file = file.path(out_dir, paste0(prefix_used_in_traderpros, ".SpBrachRate.MAP.tre.png")), width = 18, height = 24, device = "png")
  }

  # Plot tree with probability of speciation
  SpCompletionProbs <- plotTree(tree = spBrRateTree,
                               node_age_bars = FALSE,
                               node_pp = FALSE,
                               tip_labels = TRUE,
                               tip_labels_size = 3,
                               color_branch_by = "p_speciation",
                               branch_color = color_completion,
                               line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(SpCompletionProbs)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(SpCompletionProbs, file = file.path(out_dir, paste0(prefix_used_in_traderpros, ".SpeciationProb.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(SpCompletionProbs, file = file.path(out_dir, paste0(prefix_used_in_traderpros, ".Speciation.MAP.tre.png")), width = 18, height = 24, device = "png")
  }

  # Plot tree with number of speciation events
  print("Reading and saving tree with number of speciation events")
  spEventsTree <- readTrees(spEventsTreeFile)
  sp_eventsTree <- plotTree(tree = spEventsTree,
                            node_age_bars = FALSE,
                            node_pp = FALSE,
                            tip_labels = TRUE,
                            tip_labels_size = 3,
                            color_branch_by = "branch_speciation_events",
                            branch_color = color_n_speciation,
                            line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(sp_eventsTree)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(sp_eventsTree, file = file.path(out_dir, "traderpros.SpEvents.tree.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(sp_eventsTree, file = file.path(out_dir, "traderpros.SpEvents.tree.png"), width = 18, height = 24, device = "png")
  }

  # Plot ancestral states tree with hidden states
  ancStatesTree_cond <- processAncStates(ancStatesTreeFile_cond,
                                         state_labels = c("0" = anc_state_labels[1],
										                  "1" = anc_state_labels[2],
														  "2" = anc_state_labels[3],
														  "3" = anc_state_labels[4]))
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

  plot(ancStatesPlot)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(ancStatesPlot, file = file.path(out_dir, "AncEstTree.cond.MAP.wHiddenn.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(ancStatesPlot, file = file.path(out_dir, "AncEstTree.cond.MAP.wHiddenn.png"), width = 18, height = 24, device = "png")
  }

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
                    branch_color = color_birth,
                    line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(birth)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(birth, file = file.path(out_dir, "BirthRateTree.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(birth, file = file.path(out_dir, "BirthRateTree.png"), width = 18, height = 24, device = "png")
  }

  ## Plot tree with death rates
  extinction <- plotTree(tree = branchTree, 
                         node_age_bars = FALSE,
                         node_pp = FALSE, 
                         tip_labels = TRUE,
                         tip_labels_size = 5,
                         color_branch_by = "avg_mu",
                         branch_color=color_death,
                         line_width = 2.5) + 
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(extinction)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(extinction, file = file.path(out_dir, "extinctionTree.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(extinction, file = file.path(out_dir, "extinctionTree.png"), width = 18, height = 24, device = "png")
  }

  # Plot tree with net diversification rates 
  # (if you are using a poplation tree this is the net population fragmentation)

  net <- plotTree(tree = branchTree, 
                  node_age_bars = FALSE,
                  node_pp = FALSE, 
                  tip_labels = TRUE,
                  tip_labels_size = 5,
                  color_branch_by = "net_div",
                  branch_color= color_net,
                  line_width = 2.5) + 
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(net)

  if (out_images_format == "both" || out_images_format == "pdf") {
    ggsave(net, file = file.path(out_dir, "netDivTree.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(net, file = file.path(out_dir, "netDivTree.png"), width = 18, height = 24, device = "png")
  }

  # Process traces and plot custom trace plots
  print("Reading traces and saving posterior density figures")
  traceModel <- readTrace(traceFile, burnin = burn)
  match_strings <- c("^birth_observed\\[\\d+\\]$",
                     "^extinction_observed\\[\\d+\\]$",
                     "^stateSpecificRates\\[\\d+\\]$",
                     "^transition_rates\\[\\d+\\]$",
                     "^birth_hidden\\[\\d+\\]$",
                     "^extinction_hidden\\[\\d+\\]$")

  names_to_save <- c("birth_obs", "extinction_obs", "Sp_completion",
                     "State_Transition", "Birth_Hiden", "Extinction_Hidden")
c(anc_state_colors[2],anc_state_colors[4])
  
  colors_to_use <- cbind(c(anc_state_colors[2],anc_state_colors[4]), c(anc_state_colors[2],anc_state_colors[4]),
                         c(anc_state_colors[2],anc_state_colors[4]), 
                         color_transition, 
                         c(anc_state_colors[1],anc_state_colors[3]),
                         c(anc_state_colors[2],anc_state_colors[4]))
  #colors_to_use <- cbind(color_birth, color_death, color_completion, 
  #                      color_transition, color_birth,color_death   )

  for (i in 1:length(match_strings)) {
    custom_plot_trace(trace_object = traceModel,
                      variable_match = match_strings[i],
                      colors = colors_to_use[, i],
                      labels_name = c(paste(names_to_save[i], "_1", sep = ""), paste(names_to_save[i], "_2", sep = "")),
                      fig_out_name = file.path(out_dir, paste0(prefix_used_in_traderpros, names_to_save[i])),
                      out_format = c("pdf", "png"))
  }

  # Summary of rates
  print("Calculating rates summary stats")
  summary_df <- summary_stats_from_trace(trace_object = traceModel,
                           variable_match = match_strings,
                           out_name = file.path(out_dir, paste0(prefix_used_in_traderpros, ".RatesSummStats.csv")))

  # rjMCMC probabilities
  print("Saving reversible jumping probabilities")
  rj_variables = c("is_bd_state_dependent",
                  "is_birth_state_dependent",
                  "is_death_state_dependent",
                  "is_reversible",
                  "is_spCompletion_state_dependent",
                  "is_death_hidden",
                  "is_birth_hidden")

  for (i in 1:length(rj_variables)) {
    plot_rj(trace_object = traceModel,
            variable=rj_variables[i],
            colors = c("black", "grey"),
            labels_name = c("No", "Yes"),
            fig_out_name = file.path(out_dir,
                                     paste0(prefix_used_in_traderpros,
                                     "_" , rj_variables[i])
                                     ),
            out_format = out_images_format,
            burn_in = 0.1)
  }

  # Return the list of created plots and dataframe 
  return(list(
    traceModel = traceModel,
    SpCompletionRate = SpCompletionRate,
    SpCompletionProbs = SpCompletionProbs,
    sp_eventsTree = sp_eventsTree,
    ancStatesPlot = ancStatesPlot,
    net = net,
    extinction = extinction,
    birth = birth,
    summary_df = summary_df
  ))
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
args_list <- list(
  out_dir = args[1],
  in_dir = args[2],
  prefix_used_in_traderpros = args[3],
  burn = as.numeric(args[4]),
  color_completion = strsplit(args[5], ",")[[1]],
  color_n_speciation = strsplit(args[6], ",")[[1]],
  color_birth = strsplit(args[7], ",")[[1]],
  color_death = strsplit(args[8], ",")[[1]],
  color_net = strsplit(args[9], ",")[[1]],
  color_transition = strsplit(args[10], ",")[[1]],
  anc_state_labels = strsplit(args[11], ",")[[1]],
  anc_state_colors = strsplit(args[12], ",")[[1]],
  path_to_custom_functions = args[13],
  out_images_format = args[14]
)

do.call(sum_and_plot_traderpros_result, args_list)