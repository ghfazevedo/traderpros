# <ins>Tra</ins>it <ins>De</ins>pendent <ins>R</ins>ates <ins>Pro</ins>tracted <ins>S</ins>peciation (TraDeRProS)

## Introduction

## TraDeRProS in RevBayes

### Getting started  

Clone this directory, start RevBayes from your [scripts](./scripts) directory. 

```R
git clone LINK

cd traderpros/scripts

rb 
```

Set seed if you want to replicate the analysis.
```R
seed(1)
```  

Set the path to inputs and outputs  

```R
workDirPath = "../"
dataDirPath = workDirPath + "data/"
outDirPath = workDirPath + "outputs/"
outputPrefix = "traderpros"
popMapPath = dataDirPath + "popmap.tsv"
treePath = dataDirPath + "BPPConstraint.MCC.CAH.nex"
traitPath = dataDirPath + "CicTroglomorphism.nex"
speciesDataPath=dataDirPath+"SpeciesMatrix.txt"
```

Create helpers variables vor the moves and monitors.
```R
moves = VectorMoves()
monitors = VectorMonitors()
```  

### The Data
We first read the tree, the population to species map, the trait data and get some variables to use later. The population to species map is a tab delimited file, with headers, which the first column is the population (the name of the tips in the tree file) and the second column is the species from which the population belongs to. If the species assignment is uncertain, the species name should be "unknown". The algorithm will average across possible assignments, including the possibility of them being new species. You can use the most probable history of speciation summarized at the end of the analyses to delimit species.  


```R
# Read the tree file
tree <- readTrees(treePath)[1]

# Save tree with branch indexes (if you want to check the probability of speciation event on a branch later)
writeNexus(tree, filename=outDirPath+outputPrefix+"_indexed.tree.nex")

# Get the number of nodes and branches
n_nodes <- tree.nnodes()
n_branches <- n_nodes -1

# Read population to species map
popMap = readDataDelimitedFile(popMapPath, header=TRUE)
```
>Note: I tried to use the ReadTaxonData() and a loop with the taxon() function to assign species to tip similar to how is done for the dnMultiSpeciesCoalescent(), but it did not work well.  

```R
# Read trait data
trait <- readDiscreteCharacterData(traitPath)

# Get number of states
n_states <- trait.getStateDescriptions().size()
```
  
### The transition matrix
We first define the prior on the global trait transition rate. Although it is possible to independently infer the transition rates here with a less informative distribution, we will use a distribution based on the results of a previous analyses using HiSSE [here](). It has been shown before that transition rates might be wrongly inferred if it is influencing the tree shape in a state dependent birth and death model (REF). It might also be possible to use a birth and death tree process jointly with our protracted speciation to estimate rates in a trait dependent protracted speciation birth and death. However that could be time consuming and our data might not have enough power to infer a hidden state protractes speciation completion rate (that needs to be explored). Therefore we prefer to first jointly inferred the transition, extinction (death) and population formation (birth) in a HiSSE model, and use the posterior distribution to inform our transition rate here.

```R
global_trans_rate <-0.0331
relative_transition <- v(0.8128, 0.1872)
```

>If you prefer not to estimate the transition rates in a previous SSE analysis, you can co-estimate the transition rates together with the protracted speciation completion rate using the code below.
>>```
>>>traRate_expoDistrib_parameter <- 10
>>>global_trans_rate ~  dnExp(traRate_expoDistrib_parameter)
>>>moves.append( mvScale( global_trans_rate, weight=1 ) )
>>```
>  
>Also, if you are estimating the rates you should draw the relative transition rates from a distribution, such as a flat Dirichlet distribution. We can also create a reversible jump to test for reversibility in your trait
>>```R
>>n_rates = n_states * (n_states-1)
>>relative_transition ~ dnReversibleJumpMixture( simplex(1,0), dnDirichlet(rep(1, n_states)), p=0.5 )
>>
>>moves.append( mvRJSwitch(relative_transition , weight=1.0) )
>>moves.append( mvDirichletSimplex( relative_transition, weight=1 ) )
>>
>>is_troglo_reversible := ifelse( relative_transition == simplex(1,0), 0.0, 1.0)
>>```


We can now create our transition rate matrix.
```R
Q := fnFreeK(relative_transition, rescale=TRUE)
```

We also  create a prior for the root state frequency. We ca assume equal probabilities, sample from a prior distribution or fix it. We based on our previous HiSSE result to set it to a fixed value.
```R
rf <- simplex(0.6009,0.3991)
```
>Use the code below to estimate the root frequencies
>```
>rf_prior <- rep(1, n_states)
>rf ~ dnDirichlet( rf_prior )
>moves.append( mvDirichletSimplex( rf, weight=1 ) )
>```

We will use data augmentation to sample the trait history evolution along the tree, similar to what is done [here](https://revbayes.github.io/tutorials/cont_traits/state_dependent_bm.html).  

Now we can create our augmented data matrix DAG node of our model that sample the trait evolution history. We clamp the node to our observed trait data.
```R
trait_evol ~ dnPhyloCTMCDASiteIID(tree,
                                  Q,
                                  branchRates=global_trans_rate,
                                  type="Standard",
                                  nSites=1,
                                  rootFrequencies=rf)

trait_evol.clamp(trait)
```

We need to create moves for the character history along the branches.
```R
moves.append( mvCharacterHistory(ctmc=trait_evol, 
                                 qmap_site=Q,
                                 graph="node",
                                 proposal="rejection",
                                 weight=round(n_nodes/2)) )
moves.append( mvCharacterHistory(ctmc=trait_evol,
                                 qmap_site=Q,
                                 graph="branch",
                                 proposal="rejection",
                                 weight=round(n_branches/2)) )
```

### Speciation Completion Rates
We sample a global rate from an exponential distribution.

```R
spcompletion_expoDistParam <- 1
SpeciationComplRate ~ dnExponential(spcompletion_expoDistParam)
```  
We need to create moves to the speciation completion rates.
```R
moves.append(mvScale(SpeciationComplRate, weight=1))
```

The global rate is the sum of each state specific rate. To determine the vector of state specific rates, we first draw proportional weights from a Dirichlet distribution with flat prior. Then, we multiply each rate by the global rate to have each state specific rates. We use a reversible jump to test for state dependent rates (that is, if the two rates are the same or not).

```R
stateSpecificRateWeights ~ dnReversibleJumpMixture( simplex(rep(1, n_states)), dnDirichlet(rep(1, n_states)), p=0.5 )

moves.append( mvRJSwitch(stateSpecificRateWeights, weight=1.0) )
moves.append(mvSimplex(stateSpecificRateWeights, weight=1))

is_spCompletion_state_dependent := ifelse( stateSpecificRateWeights == simplex(rep(1, n_states)), 0.0, 1.0)

# If you do not want to use a reversible jump, just use the commented code below, instead.
#stateSpecificRateWeights ~ dnDirichlet(rep(1, n_states))
#moves.append(mvSimplex(stateSpecificRateWeights, weight=1))

stateSpecificRates := stateSpecificRateWeights * SpeciationComplRate
```

The branch specific rates are determined by the relative time spent on each trait in that branch. The method `trait_evol.relativeTimeInStates(i,1)` calculates the proportion of time the characater 1 spent on each state in the branch *i*.
```R
for(i in 1:n_branches) {
    state_branch_rate[i] := sum(trait_evol.relativeTimeInStates(i,1) * stateSpecificRates)
}
```

## The Trait Dependent Rates Protracted Speciation Tree
We first use our population map to replace the tip names for species names.
```R
TreeSpeciesTip <- tree
for (i in 1:popMap.size()){
    TreeSpeciesTip.setTaxonName(popMap[i][1], popMap[i][2] )
    }
# We save the tree with species names to check.
write(TreeSpeciesTip, filename=outDirPath+outputPrefix+".SpeciesName.tree" )
```

Our model is similar to the one used in [DELINEATE](https://jeetsukumaran.github.io/delineate/) ([Sukumaran et al. 2021][1]) in the way that assumes that the number of speciation events in each branch is Poisson distributed. The difference is that the rates in each branch will vary according to the state in that branch instead of a single global rate. That is, ```N_speciation_events_branch[i] ~ Poisson(state_branch_rate[i]*branch_length[i])``` where *i* is the branch in the tree, and *state_branch_rate[i]* is the rate accounting for the time spent in each state in the branch *i* (created in the step above). We know that speciation completion did not happened on branches that connect populations, so we can clamp the number of speciation events to 0. If the branches connects populations of two different species, at least one speciation event must have happened along the path connecting these two heterospecific populations. Branches that connect populations with unknown status or that connects more than two species are not used, since we cannot constrain the number of speciation events in each branch [FIGURES]().  The number of speciation events on those branches can be co-estimated with the rates for species delimitation.  

> Note: The code should work with paraphyletic species, although it is necessary to explore the effects of paraphyletic species constraints more clearly. Paraphyletic species can be inferred, though, likely with no problem.

First we create an stochastic variable to hold the number of speciation events on each branch
```R
for (i in 1:(n_branches)){
    branch_speciation_events[i] ~ dnPoisson(state_branch_rate[i]*TreeSpeciesTip.branchLength(i))
}
```

We now iterate over the tree nodes to clamp the number of speciation events to zero if the node has child branches with descendants belonging to the same species and create moves o the number of speciation events on the other branches.

```R
iteration=1
for (i in (tree.ntips()+1):n_nodes){
    j=iteration++
    # Visit every internal node i of the tree and get child branches indices and branch lengths
    c1[j] = TreeSpeciesTip.child(i,1)
    c2[j]  = TreeSpeciesTip.child(i,2)

    # Get the descendants of each child to know if they connect same or different species
    c1descendants[j] = TreeSpeciesTip.getDescendantTaxa(c1[j])
    c2descendants[j] = TreeSpeciesTip.getDescendantTaxa(c2[j])

    # Get a list of descendant species to know if node connects population of same or different species
    c1splist[j] = [c1descendants[j][1].getName()]
    if (c1descendants[j].size() > 1){
        for (n in 2:c1descendants[j].size()){
        c1splist[j].append([c1descendants[j][n].getName()])
        }
    }
    c1splist[j].unique()

    c2splist[j] = [c2descendants[j][1].getName()]
    if (c2descendants[j].size() > 1){
        for (n in 2:c2descendants[j].size()){
        c2splist[j].append([c2descendants[j][n].getName()])
        }
    }
    c2splist[j].unique()

    # If node connects populations of same species, clamp the observed value to 0 on both daughter branches to indicate that no speciation can be observed. It is only applicable if the species has known assignments. If node connects population of different species we assign moves to it so the number of speciation events can be sampled during MCMC.
    if ((c1splist[j].size() == 1 & c2splist[j].size() == 1 & c1splist[j].contains("unknown")==FALSE & c2splist[j].contains("unknown")==FALSE) & ( c1splist[j][1] == c2splist[j][1])) {
        branch_speciation_events[c1[j]].clamp(0)
        branch_speciation_events[c2[j]].clamp(0)
    } else {
        moves.append(mvRandomNaturalWalk(branch_speciation_events[c1[j]], weight=1))
        moves.append(mvRandomNaturalWalk(branch_speciation_events[c2[j]], weight=1))
        # we set a initial value just to facilitate initialization of the MCMC
        branch_speciation_events[c1[j]].setValue(1)
        branch_speciation_events[c2[j]].setValue(1)
    }
}
```

We can add a reversible jump to track the probability of having speciation events on branches of interest. You can check the branch index in the [indexed tree created in the beginning of the tutorial](outputs/traderpros_indexed.tree.nex). This is optional, since you can see the species delimitation result in the MAP tree with speciation events at the end. Also, if you have many unknown populations, the number of branches can be large and the interpretation might be tricky. However this approach might be interesting for some cases if you want to investigate the support in favor of a specific hypothesis of species limits around a few branches. Here we will test the hypothesis of *Cicurina neovespera* being conspecific with *C. bullis*, that is, no speciation event on branches number 115 (*C. neovespera* branch) and 140 (MRCA of all *C. bullis* populations).


```R
branch_speciation_events[115] ~ dnReversibleJumpMixture( 0, dnPoisson(state_branch_rate[115]*TreeSpeciesTip.branchLength(115)), p=0.5 )
moves.append( mvRJSwitch(branch_speciation_events[115], weight=1.0) )
speciation_br115 := ifelse(branch_speciation_events[115] == 0, 0, 1)
branch_speciation_events[140] ~ dnReversibleJumpMixture( 0, dnPoisson(state_branch_rate[140]*TreeSpeciesTip.branchLength(140)), p=0.5 )
moves.append( mvRJSwitch(branch_speciation_events[140], weight=1.0) )
speciation_br140 := ifelse(branch_speciation_events[140] == 0, 0, 1)
```
> Note: The moves when the node is in Poisson distribution were already set on the previous part of the code.  

Now that we clamped the observed branches with no speciation, we need to constrain our search of speciation histories to speciation events on branches that makes our species limits congruent with the observed assignments of populations to species.  
To do that, we first create a new tree which the branch lengths corresponds to the number of speciation events sampled above. In this tree, branches connecting populations of same species will have a length of 0  
```R
# This tree is used mainly to create a binary speciation tree below that is used for constraining the tree search. It can also be used to check mean (or median) number of speciation events on the branches. 
topologyPopNames <- readBranchLengthTrees(treePath)[1]
topologyPopNames.renumberNodes(tree)
SpeciationBranchTreePopNames := fnTreeAssembly(topologyPopNames, branch_speciation_events)

# We modify the tree above to represent a tree with binary speciation events; that is,  branch length=0 if no speciation and 1 if one or more speciation events.
for (i in 1:branch_speciation_events.size()){
    if (branch_speciation_events[i] > 0 ){
        branch_speciation_binary[i] <- abs(1)
    } else {
        branch_speciation_binary[i] <- abs(0)
    }
}
SpeciationBranchTreeBinary := fnTreeAssembly(SpeciationBranchTreePopNames, branch_speciation_binary)
```

Now we constrain our search of speciation events trees to only sample speciation events that are compatible with the observed species, i.e., trees which branches that connects different species having a length of 0 would have also 0 probability.
For this we use a species data matrix in which the trait states indicate the species assignments, and populations with unknown identity are coded with "?" (see the [species data matrix file](data/SpeciesMatrix.txt)). This is exactly the same information as the population to species map, but numbers are used instead of species names/codes.
We then create a phyloCTMC object which "evolves" the species assignments in the tree with speciation events as brach lengths.  

```R
# We need to inform the maximum number of species in the tree. This includes the "unknown" species that could possible be new too. For example, if you have 20 known species and 6 populations which you do not known the assignments and suspects that each of these populations could be a new species, you have a maximum of 26 species in your tree.
max_number_of_species=26

# read the matrix with species data:
speciesData = readDelimitedCharacterData(file=speciesDataPath, stateLabels=max_number_of_species )

# Create a transition matrix with equal rates of transition between all species labels
Qs <- fnFreeK(rep(1, max_number_of_species^2), rescale=FALSE)

speciation ~ dnPhyloCTMC(tree=SpeciationBranchTreePopNames, Q=Qs, type="NaturalNumbers")
```
We clamp the object to the species label matrix data. 
```R
speciation.clamp(speciesData)
```

We can now finalize our model 
```R
mymodel = model(branch_speciation_events)
```

### The MCMC
Define the number of generations, burn in and print frequency
```R
## Number of generations
n_gen = 500000
## Print to screen every x generations
n_print_screen = 100
## Print parameters to file every x generations
n_print_log = 100
## Print tree to file every x generations
n_print_tree=100
## Number of independent runs
n_runs = 1
## Create checkpoint file every x generations
n_check_interval = 1000
## Burn in
burnin_percentage = 0.10 
```

Monitor variables of interest
```R
# Monitor all parameters
monitors.append( mnModel(file=outDirPath+outputPrefix+".model.log", 
                          printgen=n_print_log) )

# monitor rates
monitors.append( mnFile(SpeciationComplRate, stateSpecificRates,
                        filename=outDirPath+outputPrefix+".Rates.log",
                        printgen=n_print_log) )

# monitor reversible jump hypotheses to test for equal rates and for speciation events on branches of interest.
monitors.append( mnFile(is_spCompletion_state_dependent,
                        speciation_br115,
                        speciation_br140
                        filename=outDirPath+outputPrefix+".RJ.log",
                        printgen=n_print_log) )


# Monitor tree with speciation events and branch lengths
monitors.append( mnFile(SpeciationBranchTreePopNames, 
                        filename=outDirPath+outputPrefix+".EventsTree.log", 
                        printgen=n_print_tree) )

monitors.append( mnFile(SpeciationBranchTreeBinary, 
                        filename=outDirPath+outputPrefix+".EventsTreeBinary.log", 
                        printgen=n_print_tree) )

# Monitor trait evolution (I am saving this, but it needs some coding for summarizing the character history results that I do not know)
monitors.append( mnFile( tree, filename=outDirPath+outputPrefix+".tre", 
                        printgen=n_print_tree) )

monitors.append( mnCharacterHistorySummary( filename=outDirPath+outputPrefix+".history.txt",
                                            ctmc=trait_evol,
                                            tree=tree,
                                            printgen=n_print_tree ) )

monitors.append(mnCharHistoryNewick(ctmc=trait_evol,
                                    tree=tree,
                                    filename=outDirPath+outputPrefix+".history.nwck", 
                                    printgen=n_print_tree ) )



# Monitor tree with speciation completion branch rates
monitors.append( mnExtNewick(
                 filename=outDirPath+outputPrefix+".SpCompletionRates.trees",
                 isNodeParameter=TRUE,
                 printgen=n_print_tree,
                 separator=TAB,
                 tree=tree,
                 state_branch_rate
                 ))

# Monitor the original tree with branch lengths in time together with the number of speciation events associated with each branch.
monitors.append( mnExtNewick(
                 filename=outDirPath+outputPrefix+".SpeciationEvents.trees",
                 isNodeParameter=TRUE,
                 printgen=n_print_tree,
                 separator=TAB,
                 tree=tree, 
                 branch_speciation_events
                 ) )

# Monitor only the speciation events on each branch (this can be used to check the posterior distribution of number of speciation events and guide species delimitation)
monitors.append( mnFile(branch_speciation_events,
                        filename=outDirPath+outputPrefix+".SpeciationEvents.log",
                        printgen=n_print_log) )

```


Now we can create the MCMC object and run

```R
mymcmc = mcmc(mymodel, monitors, moves, nruns=n_runs, combine="mixed")

mymcmc.run(n_gen, checkpointFile=outDirPath+outputPrefix+".checkpoint", 
           checkpointInterval= n_check_interval )
```

## Process results
We can now summarize the results

```R
# Summarize original tree with number of speciation events associated with each branch. This is used for species delimitation.
treeTrace = readTreeTrace(outDirPath+outputPrefix+".SpeciationEvents.trees",
            burnin=burnin_percentage )
mapTree(treeTrace,
        file=outDirPath+outputPrefix+".SpeciationEvents.MAP.tre",
        mean=FALSE
        )
# The output needs some modification. We use a bash command to do it. If you are using windows, you may need to modify it.
system("sed -i 's/=1}/}/g' ../outputs/*.SpeciationEvents.MAP.tre")


# Summarize tree with speciation completion rates
treeTraceRates = readTreeTrace(outDirPath+outputPrefix+".SpCompletionRates.trees",
            burnin=burnin_percentage )
mapTree(treeTraceRates,
        file=outDirPath+outputPrefix+".SpCompletionRates.MAP.tre")
```

> You can also summarize the tree with branch lengths in number of speciation events. Note that the branch length would be the median (or mean if you choose) number of speciation events, and it might not be equal to the MAP number of events summarized above. You can also use this for guiding your conclusions about species limits. 
>```R
>treeTrace1 = readTreeTrace(outDirPath+outputPrefix+".EventsTree.log",
>            burnin=burnin_percentage )
>mapTree(treeTrace1,
>        file=outDirPath+outputPrefix+".EventsTree.MAP.tre",
>        ccAges=TRUE,
>        ccp=TRUE,
>        conditionalAges=TRUE,
>        mean=FALSE,
>        sampledAncestors=FALSE,
>        positiveBranchLengths=FALSE)
>
>treeTraceBinary = readTreeTrace(outDirPath+outputPrefix+".EventsTreeBinary.log",
>            burnin=burnin_percentage )
>mapTree(treeTraceBinary,
>        file=outDirPath+outputPrefix+".EventsTreeBinary.MAP.tre",
>        ccAges=TRUE,
>        ccp=TRUE,
>        conditionalAges=TRUE,
>        mean=FALSE,
>        sampledAncestors=FALSE,
>        positiveBranchLengths=FALSE)
>```

## Plot results with RevGadgets

Open R


```R
library(ggplot2)
library(RevGadgets)
library(coda)
library(RColorBrewer)
library(viridis)

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
```

The the branch specific rates can be seen in the figure below.
![Tree with speciation completion rates](./outfigs/traderpros.SpCompRates.tree.png) 

You can compare with the trait evolution estimated [here]() with the HiSSE model.
![Tree with trait evolution](./outfigs/AncEstTree.cond.MAP.png) 
  

You can see the number of estimated speciation events on each branch. This results is also the maximum a posteriori species delimitation, if you had populations with unknown species assignments. 

```R
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
```

Branches connected by 0 speciation events belongs to the same species. This is the most probable species limits.

![Tree with estimated number of speciation events](./outfigs/traderpros.SpEvents.tree.png) 

  
You can also plot the posterior distribution of parameters.
```R
# Read Trace and get summary stats
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


# Plot state specific speciation completion rate posterior distribution
plotSpCompl <- plotTrace(trace = traceModel, 
                       vars = c("stateSpecificRates[1]","stateSpecificRates[2]","SpeciationComplRate")) 
plotSpCompl


ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.pdf", sep="")) 
ggsave(file=paste(outFigsDir, "traderpros.StateSpecificRates.posterior.png", sep=""))
```

This is the posterior distribution for each state specific rates and the global rate
![state specific speciation completion rate posterior distribution ](./outfigs/traderpros.StateSpecificRates.posterior.png) 

```R
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
```
This is the probability of the model being state dependent 
![state specific speciation completion rate posterior distribution ](./outfigs/%20traderpros.RJSPSompStateDep.probability.png) 
  

Plot probability of speciation events on target branches
```R
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

```

This is the reversible jump result testing the probability of speciation happening on branch 115.
![Speciation on Branch 115](./outfigs/%20traderpros.spEventProbBr115.png)


You can do the same for the other branches
```R
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
```

Other possble plots below.
```R
spEventPosteriorBr115 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[115]"));spEventPosteriorBr115

spEventPosteriorBr140 <- plotTrace(trace = traceModel, 
                              vars = c("branch_speciation_events[140]"));spEventPosteriorBr140
```


# References
Sukumaran J, Holder MT, Knowles LL (2021) Incorporating the speciation process into species delimitation. [PLOS Computational Biology 17(5): e1008924](ttps://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008924) 
