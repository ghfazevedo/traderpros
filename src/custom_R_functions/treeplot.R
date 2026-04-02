  # Load required libraries
  library(ggplot2)
  library(RevGadgets)
  library(coda)
  library(RColorBrewer)
  library(viridis)


  # File paths for TraderPros outputs
  spBrRateTreeFile <- list.files(path = input, pattern = ".Protracted.MAP.tre", full.names = TRUE)[1]
  ancStatesTreeFile_cond <- list.files(path = input, pattern = ".traits.MAP.cond.tree", full.names = TRUE)[1]
  ancStatesTreeFile_marg <- list.files(path = input, pattern = ".traits.MAP.marg.tree", full.names = TRUE)[1]
  branchRatesFile <- list.files(path = input, pattern =".BirthDeathBrRates.log", full.names = TRUE)[1]

  # Plot tree with speciation rates
  print("Reading and saving Protracted tree with annotations")
  spBrRateTree <- readTrees(spBrRateTreeFile)
  SpCompletionRate <- plotTree(tree = spBrRateTree,
                               node_age_bars = FALSE,
                               node_pp = FALSE,
                               tip_labels = TRUE,
                               tip_labels_size = 3,
                               color_branch_by = "state_branch_rate",
                               branch_color = rate_colors,
                               line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(SpCompletionRate)

    ggsave(SpCompletionRate, file = file.path(output, paste0(prefix_used_in_traderpros, ".SpBrachRate.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
    ggsave(SpCompletionRate, file = file.path(output, paste0(prefix_used_in_traderpros, ".SpBrachRate.MAP.tre.tiff")), width = 18, height = 24, device = "png")


  # Plot tree with probability of speciation
  SpCompletionProbs <- plotTree(tree = spBrRateTree,
                               node_age_bars = FALSE,
                               node_pp = FALSE,
                               tip_labels = TRUE,
                               tip_labels_size = 3,
                               color_branch_by = "p_speciation",
                               branch_color = rate_colors,
                               line_width = 2.5) +
    ggplot2::theme(legend.position.inside = c(.1, .9))

  plot(SpCompletionProbs)

    ggsave(SpCompletionProbs, file = file.path(output, paste0(prefix_used_in_traderpros, ".SpeciationProb.MAP.tre.pdf")), width = 18, height = 24, device = "pdf")
    ggsave(SpCompletionProbs, file = file.path(output, paste0(prefix_used_in_traderpros, ".Speciation.MAP.tre.png")), width = 18, height = 24, device = "tiff")


  # Plot tree with number of speciation events
  print("Reading and saving tree with number of speciation events")
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

  plot(sp_eventsTree)

    ggsave(sp_eventsTree, file = file.path(output, "traderpros.SpEvents.tree.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(sp_eventsTree, file = file.path(output, "traderpros.SpEvents.tree.png"), width = 18, height = 24, device = "tiff")


  # Plot ancestral states tree with hidden states
  hidden_letters <- toupper(letters[1:n_hidden_states])
  states <- expand.grid(
    obs = 0:(n_observed_states - 1),
    hid = hidden_letters
  )
  anc_state_labels <- paste0(states$obs, states$hid)

  if (!is.null("anc_state_colors")){
    if (length(anc_state_labels)==4){
      anc_state_colors <- c( "#56B4E9", "#EE7733", "#0072B2", "#CC3311")
    }else {
      anc_state_colors <- turbo(length(anc_state_labels))
    }
  }

  ancStatesTree_cond <- processAncStates(ancStatesTreeFile_cond,
                                         state_labels = anc_state_labels)

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
    ggsave(ancStatesPlot, file = file.path(output, "AncEstTree.cond.MAP.wHiddenn.pdf"), width = 18, height = 24, device = "pdf")
  }
  if (out_images_format == "both" || out_images_format == "png") {
    ggsave(ancStatesPlot, file = file.path(output, "AncEstTree.cond.MAP.wHiddenn.png"), width = 18, height = 24, device = "png")
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

    ggsave(birth, file = file.path(output, "BirthRateTree.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(birth, file = file.path(output, "BirthRateTree.png"), width = 18, height = 24, device = "png")


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


    ggsave(extinction, file = file.path(output, "extinctionTree.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(extinction, file = file.path(output, "extinctionTree.png"), width = 18, height = 24, device = "png")

  # Plot tree with net diversification rates 
  # (if you are using a population tree this is the net population fragmentation)

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

    ggsave(net, file = file.path(output, "netDivTree.pdf"), width = 18, height = 24, device = "pdf")
    ggsave(net, file = file.path(output, "netDivTree.png"), width = 18, height = 24, device = "png")

