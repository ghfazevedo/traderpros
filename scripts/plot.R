library(ggplot2)
library(RevGadgets)
library(RColorBrewer)
library(viridis)

# output and input directories
outFigsDir = "../outfigs/"
inputFilesDir = "../outputs/"

# files
treeFile = paste(inputFilesDir,"traderpros.SpCompletionRates.MAP.tre", sep = "")
eventsTreeFile = paste(inputFilesDir, "traderpros.SpeciationEvents.MAP.tre", sep = "")
traceFile = paste(inputFilesDir,  "traderpros.model.log", sep = "")

# burnin percentege
burn = 0.10

# Plot tree with rates
tree = readTrees(treeFile)

pal_name <- "RdYlBu"
ratecolors <- brewer.pal(n = 11, name = pal_name)

SpCompletionRate <- plotTree(tree = tree, 
                   node_age_bars = FALSE,
                   node_pp = FALSE, 
                   tip_labels = TRUE,
                   tip_labels_size = 3,
                   color_branch_by = "state_branch_rate",
                   branch_color=c(ratecolors[11], ratecolors[2]),
                   line_width = 0.8) + 
  ggplot2::theme(legend.position=c(.1, .9));SpCompletionRate

ggsave(SpCompletionRate,file=paste(outFigsDir, "traderpros.SpCompRates.tree.pdf", sep=""), width = 18, height = 24) 
ggsave(SpCompletionRate,file=paste(outFigsDir, "traderpros.SpCompRates.tree.png", sep=""), width = 18, height = 24) 


# Plot number of speciation events
eventsTree = readTrees(eventsTreeFile)

pal_name <- "RdYlBu"
colors <- brewer.pal(n = 11, name = pal_name)

sp_eventsTree <- plotTree(tree=eventsTree,
         node_age_bars = FALSE,
         node_pp = FALSE,
         tip_labels = TRUE,
         tip_labels_size = 3,
         color_branch_by = "branch_speciation_events",
         branch_color=c(colors[10], colors[2]),
         line_width = 0.8) +
  ggplot2::theme(legend.position=c(.1, .9));sp_eventsTree         

ggsave(sp_eventsTree,file=paste(outFigsDir, "traderpros.SpEvents.tree.pdf", sep=""), width = 18, height = 24) 
ggsave(sp_eventsTree,file=paste(outFigsDir, "traderpros.SpEvents.tree.png", sep=""), width = 18, height = 24) 


# Plot Trace
traceModel <- readTrace(traceFile)


plotSpCompl <- plotTrace(trace = traceModel, 
                       vars = c("stateSpecificRates[1]","stateSpecificRates[2]"));plotSpCompl 

ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.pdf", sep="")) 
ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.png", sep=""))


# Plot reversible jump hypotheses
plotSPComplRJ <- plotTrace(trace = traceModel, 
                             vars = c("is_spCompletion_state_dependent"));plotSPComplRJ
ggsave(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.posterior.pdf", sep="")) 
ggsave(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.posterior.png", sep="")) 


plotSPComplPost <- plotTrace(trace = traceModel, 
                           vars = c("SpeciationComplRate"));plotSPComplPost
ggsave(file=paste(outFigsDir, "traderpros.SPCompletion.posterior.pdf", sep="")) 


# Plot troglomorphism Reversible jump probabilities
plotTrogloRJ <- plotTrace(trace = traceModel, 
                           vars = c("is_troglo_reversible"));plotTrogloRJ
ggsave(file=paste(outFigsDir, "traderpros.RJTroglRever.posterior.pdf", sep="")) 


plotRelTrasn <- plotTrace(trace = traceModel, 
                          vars = c("relative_transition[1]","relative_transition[2]"));plotRelTrasn
ggsave(file=paste(outFigsDir, "traderpros.RelatTrans.posterior.pdf", sep="")) 

plotTrasn <- plotTrace(trace = traceModel, 
                          vars = c("global_trans_rate"));plotTrasn
ggsave(file=paste(outFigsDir, "traderpros.TransRates.posterior.pdf", sep="")) 



# Plot speciation events for target branches
spEventProbBr115 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[115]"));spEventProbBr115

spEventProbBr116 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[116]"));spEventProbBr116

spEventProbBr117 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[117]"));spEventProbBr117

spEventProbBr118 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[118]"));spEventProbBr118

spEventProbBr119 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[119]"));spEventProbBr119

spEventProbBr127 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[127]"));spEventProbBr127

spEventProbBr128 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[128]"));spEventProbBr128

spEventProbBr129 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[129]"));spEventProbBr129

spEventProbBr141 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[141]"));spEventProbBr141

spEventProbBr140 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[140]"));spEventProbBr140
