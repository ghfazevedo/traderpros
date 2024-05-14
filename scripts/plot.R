library(ggplot2)
library(RevGadgets)
library(coda)
library(RColorBrewer)
library(viridis)
#library(lessR)

# output and input directories
outFigsDir = "../outfigs/"
inputFilesDir = "../outputs/"

ifelse(!dir.exists(file.path(outFigsDir)), dir.create(file.path(outFigsDir)), "Output directory exists.Rename it or make sure files within will have diferent names")


# files
treeFile = paste(inputFilesDir,"traderpros.SpCompletionRates.MAP.tre", sep = "")
eventsTreeFile = paste(inputFilesDir, "traderpros.SpeciationEvents.MAP.tre", sep = "")
traceFile = paste(inputFilesDir,  "traderpros.model.log", sep = "")
ratesTraceFile = paste(inputFilesDir,  "traderpros.Rates.log", sep = "")

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
         branch_color=c(colors[9], colors[1]),
         line_width = 0.8) +
  ggplot2::theme(legend.position=c(.1, .9));sp_eventsTree         

ggsave(sp_eventsTree,file=paste(outFigsDir, "traderpros.SpEvents.tree.pdf", sep=""), width = 18, height = 24) 
ggsave(sp_eventsTree,file=paste(outFigsDir, "traderpros.SpEvents.tree.png", sep=""), width = 18, height = 24) 


# Plot Trace
traceModel <- readTrace(traceFile)
traceRate <- readTrace(ratesTraceFile)

trace_quant_MCMC <- as.mcmc(traceRate[[1]])
ESS<-as.data.frame(effectiveSize(trace_quant_MCMC))
traceplot(trace_quant_MCMC)

RatesSummary<-as.data.frame(summarizeTrace(trace = traceRate, 
                                           vars =  c("stateSpecificRates[1]",
                                                     "stateSpecificRates[2]",
                                                     "SpeciationComplRate")))
RatesSummary


plotSpCompl <- plotTrace(trace = traceModel, 
                       vars = c("stateSpecificRates[1]","stateSpecificRates[2]","SpeciationComplRate")) 
plotSpCompl


ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.pdf", sep="")) 
ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.png", sep=""))


# Plot reversible jump hypotheses
plotSPComplRJ <- plotTrace(trace = traceModel, 
                             vars = c("is_spCompletion_state_dependent"));plotSPComplRJ
ggsave(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.posterior.pdf", sep="")) 
ggsave(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.posterior.png", sep="")) 


StateDependentProbability <- as.data.frame(plotSPComplRJ[[1]]$data)
StateDependentProbability$State <- replace(StateDependentProbability$State, StateDependentProbability$State == "0", "Independent")
StateDependentProbability$State <- replace(StateDependentProbability$State, StateDependentProbability$State == "1", "Dependent")


png(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.probability.png"))
pie(StateDependentProbability$Probability, labels = paste(StateDependentProbability$State, " ", round(StateDependentProbability$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.RJSPSompStateDep.probability.pdf"))
pie(StateDependentProbability$Probability, labels = paste(StateDependentProbability$State, " ", round(StateDependentProbability$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()



# Plot speciation events for target branches
spEventPosteriorBr115 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[115]"));spEventPosteriorBr115

spEventPosteriorBr116 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[116]"));spEventPosteriorBr116

spEventPosteriorBr117 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[117]"));spEventPosteriorBr117

spEventPosteriorBr118 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[118]"));spEventPosteriorBr118

spEventPosteriorBr119 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[119]"));spEventPosteriorBr119

spEventPosteriorBr127 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[127]"));spEventPosteriorBr127

spEventPosteriorBr128 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[128]"));spEventPosteriorBr128

spEventPosteriorBr129 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[129]"));spEventPosteriorBr129

spEventPosteriorBr141 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[141]"));spEventPosteriorBr141

spEventPosteriorBr140 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[140]"));spEventPosteriorBr140

spEventPosteriorBr142 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[142]"));spEventPosteriorBr140



# Plot the probability of no speciation vs speciation (rjMCMC results)

# Branch 115
spEventProbBr115 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br115"));spEventProbBr115

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr115.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr115.posterior.png", sep=""))

spEventProbBr115_table <- as.data.frame(spEventProbBr115[[1]]$data)
spEventProbBr115_table$State <- replace(spEventProbBr115_table$State, spEventProbBr115_table$State == "0", "No Speciation")
spEventProbBr115_table$State <- replace(spEventProbBr115_table$State, spEventProbBr115_table$State == "1", "Speciation")

pie(spEventProbBr115_table$Probability, main = "Probability of speciation Branch 115", labels = paste(spEventProbBr115_table$State, " ", round(spEventProbBr115_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr115.png"))
pie(spEventProbBr115_table$Probability, labels = paste(spEventProbBr115_table$State, " ", round(spEventProbBr115_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr115.pdf"))
pie(spEventProbBr115_table$Probability, labels = paste(spEventProbBr115_table$State, " ", round(spEventProbBr115_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

# Branch 116
spEventProbBr116 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br116"));spEventProbBr116

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr116.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr116.posterior.png", sep=""))

spEventProbBr116_table <- as.data.frame(spEventProbBr116[[1]]$data)
spEventProbBr116_table$State <- replace(spEventProbBr116_table$State, spEventProbBr116_table$State == "0", "No Speciation")
spEventProbBr116_table$State <- replace(spEventProbBr116_table$State, spEventProbBr116_table$State == "1", "Speciation")

pie(spEventProbBr116_table$Probability, main = "Probability of speciation Branch 116", labels = paste(spEventProbBr116_table$State, " ", round(spEventProbBr116_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr116.png"))
pie(spEventProbBr116_table$Probability, labels = paste(spEventProbBr116_table$State, " ", round(spEventProbBr116_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr116.pdf"))
pie(spEventProbBr116_table$Probability, labels = paste(spEventProbBr116_table$State, " ", round(spEventProbBr116_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

# Branch 117
spEventProbBr117 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br117"));spEventProbBr117

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr117.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr117.posterior.png", sep=""))

spEventProbBr117_table <- as.data.frame(spEventProbBr117[[1]]$data)
spEventProbBr117_table$State <- replace(spEventProbBr117_table$State, spEventProbBr117_table$State == "0", "No Speciation")
spEventProbBr117_table$State <- replace(spEventProbBr117_table$State, spEventProbBr117_table$State == "1", "Speciation")

pie(spEventProbBr117_table$Probability, main = "Probability of speciation Branch 117", labels = paste(spEventProbBr117_table$State, " ", round(spEventProbBr117_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr117.png"))
pie(spEventProbBr117_table$Probability, labels = paste(spEventProbBr117_table$State, " ", round(spEventProbBr117_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr117.pdf"))
pie(spEventProbBr117_table$Probability, labels = paste(spEventProbBr117_table$State, " ", round(spEventProbBr117_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 118
spEventProbBr118 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br118"));spEventProbBr118

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr118.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr118.posterior.png", sep=""))

spEventProbBr118_table <- as.data.frame(spEventProbBr118[[1]]$data)
spEventProbBr118_table$State <- replace(spEventProbBr118_table$State, spEventProbBr118_table$State == "0", "No Speciation")
spEventProbBr118_table$State <- replace(spEventProbBr118_table$State, spEventProbBr118_table$State == "1", "Speciation")

pie(spEventProbBr118_table$Probability, main = "Probability of speciation Branch 118", labels = paste(spEventProbBr118_table$State, " ", round(spEventProbBr118_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr118.png"))
pie(spEventProbBr118_table$Probability, labels = paste(spEventProbBr118_table$State, " ", round(spEventProbBr118_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr118.pdf"))
pie(spEventProbBr118_table$Probability, labels = paste(spEventProbBr118_table$State, " ", round(spEventProbBr118_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 119
spEventProbBr119 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br119"));spEventProbBr119

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr119.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr119.posterior.png", sep=""))

spEventProbBr119_table <- as.data.frame(spEventProbBr119[[1]]$data)
spEventProbBr119_table$State <- replace(spEventProbBr119_table$State, spEventProbBr119_table$State == "0", "No Speciation")
spEventProbBr119_table$State <- replace(spEventProbBr119_table$State, spEventProbBr119_table$State == "1", "Speciation")

pie(spEventProbBr119_table$Probability, main = "Probability of speciation Branch 119", labels = paste(spEventProbBr119_table$State, " ", round(spEventProbBr119_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr119.png"))
pie(spEventProbBr119_table$Probability, labels = paste(spEventProbBr119_table$State, " ", round(spEventProbBr119_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr119.pdf"))
pie(spEventProbBr119_table$Probability, labels = paste(spEventProbBr119_table$State, " ", round(spEventProbBr119_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 127
spEventProbBr127 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br127"));spEventProbBr127

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr127.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr127.posterior.png", sep=""))

spEventProbBr127_table <- as.data.frame(spEventProbBr127[[1]]$data)
spEventProbBr127_table$State <- replace(spEventProbBr127_table$State, spEventProbBr127_table$State == "0", "No Speciation")
spEventProbBr127_table$State <- replace(spEventProbBr127_table$State, spEventProbBr127_table$State == "1", "Speciation")

pie(spEventProbBr127_table$Probability, main = "Probability of speciation Branch 127", labels = paste(spEventProbBr127_table$State, " ", round(spEventProbBr127_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr127.png"))
pie(spEventProbBr127_table$Probability, labels = paste(spEventProbBr127_table$State, " ", round(spEventProbBr127_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr127.pdf"))
pie(spEventProbBr127_table$Probability, labels = paste(spEventProbBr127_table$State, " ", round(spEventProbBr127_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 128
spEventProbBr128 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br128"));spEventProbBr128

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr128.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr128.posterior.png", sep=""))

spEventProbBr128_table <- as.data.frame(spEventProbBr128[[1]]$data)
spEventProbBr128_table$State <- replace(spEventProbBr128_table$State, spEventProbBr128_table$State == "0", "No Speciation")
spEventProbBr128_table$State <- replace(spEventProbBr128_table$State, spEventProbBr128_table$State == "1", "Speciation")

pie(spEventProbBr128_table$Probability, main = "Probability of speciation Branch 128", labels = paste(spEventProbBr128_table$State, " ", round(spEventProbBr128_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr128.png"))
pie(spEventProbBr128_table$Probability, labels = paste(spEventProbBr128_table$State, " ", round(spEventProbBr128_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr128.pdf"))
pie(spEventProbBr128_table$Probability, labels = paste(spEventProbBr128_table$State, " ", round(spEventProbBr128_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 129
spEventProbBr129 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br129"));spEventProbBr129

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr129.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr129.posterior.png", sep=""))

spEventProbBr129_table <- as.data.frame(spEventProbBr129[[1]]$data)
spEventProbBr129_table$State <- replace(spEventProbBr129_table$State, spEventProbBr129_table$State == "0", "No Speciation")
spEventProbBr129_table$State <- replace(spEventProbBr129_table$State, spEventProbBr129_table$State == "1", "Speciation")

pie(spEventProbBr129_table$Probability, main = "Probability of speciation Branch 129", labels = paste(spEventProbBr129_table$State, " ", round(spEventProbBr129_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr129.png"))
pie(spEventProbBr129_table$Probability, labels = paste(spEventProbBr129_table$State, " ", round(spEventProbBr129_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr129.pdf"))
pie(spEventProbBr129_table$Probability, labels = paste(spEventProbBr129_table$State, " ", round(spEventProbBr129_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

# Branch 140
spEventProbBr140 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br140"));spEventProbBr140

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr140.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr140.posterior.png", sep=""))

spEventProbBr140_table <- as.data.frame(spEventProbBr140[[1]]$data)
spEventProbBr140_table$State <- replace(spEventProbBr140_table$State, spEventProbBr140_table$State == "0", "No Speciation")
spEventProbBr140_table$State <- replace(spEventProbBr140_table$State, spEventProbBr140_table$State == "1", "Speciation")

pie(spEventProbBr140_table$Probability, main = "Probability of speciation Branch 140", labels = paste(spEventProbBr140_table$State, " ", round(spEventProbBr140_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr140.png"))
pie(spEventProbBr140_table$Probability, labels = paste(spEventProbBr140_table$State, " ", round(spEventProbBr140_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr140.pdf"))
pie(spEventProbBr140_table$Probability, labels = paste(spEventProbBr140_table$State, " ", round(spEventProbBr140_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 141
spEventProbBr141 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br141"));spEventProbBr141

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr141.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr141.posterior.png", sep=""))

spEventProbBr141_table <- as.data.frame(spEventProbBr141[[1]]$data)
spEventProbBr141_table$State <- replace(spEventProbBr141_table$State, spEventProbBr141_table$State == "0", "No Speciation")
spEventProbBr141_table$State <- replace(spEventProbBr141_table$State, spEventProbBr141_table$State == "1", "Speciation")

pie(spEventProbBr141_table$Probability, main = "Probability of speciation Branch 141", labels = paste(spEventProbBr141_table$State, " ", round(spEventProbBr141_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr141.png"))
pie(spEventProbBr141_table$Probability, labels = paste(spEventProbBr141_table$State, " ", round(spEventProbBr141_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr141.pdf"))
pie(spEventProbBr141_table$Probability, labels = paste(spEventProbBr141_table$State, " ", round(spEventProbBr141_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


# Branch 142
spEventProbBr142 <- plotTrace(trace = traceModel, 
                              vars = c("speciation_br142"));spEventProbBr142

#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr142.posterior.pdf", sep="")) 
#ggsave(file=paste(outFigsDir, "traderpros.spEventProbBr142.posterior.png", sep=""))

spEventProbBr142_table <- as.data.frame(spEventProbBr142[[1]]$data)
spEventProbBr142_table$State <- replace(spEventProbBr142_table$State, spEventProbBr142_table$State == "0", "No Speciation")
spEventProbBr142_table$State <- replace(spEventProbBr142_table$State, spEventProbBr142_table$State == "1", "Speciation")

pie(spEventProbBr142_table$Probability, main = "Probability of speciation Branch 142", labels = paste(spEventProbBr142_table$State, " ", round(spEventProbBr142_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))

png(file=paste(outFigsDir, "traderpros.spEventProbBr142.png"))
pie(spEventProbBr142_table$Probability, labels = paste(spEventProbBr142_table$State, " ", round(spEventProbBr142_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()

pdf(file=paste(outFigsDir, "traderpros.spEventProbBr142.pdf"))
pie(spEventProbBr142_table$Probability, labels = paste(spEventProbBr142_table$State, " ", round(spEventProbBr142_table$Probability, 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))
dev.off()


######################################
# Species delimitation probabilities #
######################################

# Probability C. neovespera (branch 115) be the same species as C. bullis (branch 140)
P_synonym =  spEventProbBr115_table$Probability[1] * spEventProbBr140_table$Probability[1]

# Probability C. neovespera (branch 115) be different species species as C. bullis (branch 140)
P_valid =  spEventProbBr115_table$Probability[1] * spEventProbBr140_table$Probability[2] +
           spEventProbBr115_table$Probability[2] * spEventProbBr140_table$Probability[1] +
           spEventProbBr115_table$Probability[2] * spEventProbBr140_table$Probability[2]

lab = c("synonym", "not synonym")
pie(c(P_synonym,P_valid), main = "Delimitation probability", labels = paste( lab, " ", round(c(P_synonym,P_valid), 2), "%", sep=""),  border="white", col=c(colors[1],colors[11] ))


# Probability of the three possible speciation histories that would make C. neovespera a valid species
P_history_1 =  spEventProbBr115_table$Probability[1] * spEventProbBr140_table$Probability[2] 
P_history_2 =  spEventProbBr115_table$Probability[2] * spEventProbBr140_table$Probability[1] 
P_history_3 =  spEventProbBr115_table$Probability[2] * spEventProbBr140_table$Probability[2]

paste("The probability of the speciation history 1, 2 and 3 are", round(P_history_1, 3), round(P_history_2, 3), round(P_history_3, 3), "respectively", sep=" ")



# Probability new populations being a new species (branch 119, 118, 127) different from C. baronia (branch 116, 117, 128)
P_new =  spEventProbBr119_table$Probability[1] * spEventProbBr118_table$Probability[1] * spEventProbBr127_table$Probability[2] *
         spEventProbBr116_table$Probability[1] * spEventProbBr117_table$Probability[1] * spEventProbBr128_table$Probability[1] +
         spEventProbBr119_table$Probability[1] * spEventProbBr118_table$Probability[1] * spEventProbBr127_table$Probability[1] *
         spEventProbBr116_table$Probability[1] * spEventProbBr117_table$Probability[1] * spEventProbBr128_table$Probability[2] + 
         spEventProbBr119_table$Probability[1] * spEventProbBr118_table$Probability[1] * spEventProbBr127_table$Probability[2] *
         spEventProbBr116_table$Probability[1] * spEventProbBr117_table$Probability[1] * spEventProbBr128_table$Probability[2]

P_new


P_notNew =  spEventProbBr119_table$Probability[1] * spEventProbBr118_table$Probability[1] * spEventProbBr127_table$Probability[1] *
            spEventProbBr116_table$Probability[1] * spEventProbBr117_table$Probability[1] * spEventProbBr128_table$Probability[1]


paste("Probability of being new: ", P_new)
paste("Probability of being the same: ", P_notNew)







