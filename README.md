# <ins>Tra</ins>it <ins>De</ins>pendent <ins>R</ins>ates <ins>Pro</ins>tracted <ins>S</ins>peciation (Traderpros)

## Introduction
Traderpros is a package of python programs to estimate parameters, delimit species and simulate data under the Trait Dependent Rates Protracted Speciation model in a Bayesian framework as proposed in [Azevedo et al. (2024)](). There are also programs to summarize results and plot graphs. Traderpros programs allow users to simply run the model and check results using unix command line interface.  


### The Traderpros Model
The protracted speciation model was first proposed by [REF]() and have been further developed ([REF]()) and used for understand diversification and for species delimitation ([REF]()). Our approach consist of modeling the population tree branching pattern as a Birth and Death model (or any derivatives) and the number of speciation completion events on the brach as extended Poisson process similar to the one used in [DELINEATE](https://jeetsukumaran.github.io/delineate/) ([Sukumaran et al. 2021][1]). The difference is that the speciation completion rate in each branch will vary according to the state in that branch instead of a single global rate. That is, ```N_speciation_events_branch[i] ~ Poisson(state_branch_rate[i]*branch_length[i])``` where *i* is the branch in the tree, and *state_branch_rate[i]* is the rate accounting for the time spent in each trait state in the branch *i*. To sample the trait history evolution along the tree and obtain the the time spent in each states on a branch, we use a data augmentation approach in similar way that is done for continuous trait evolution by [May and Moore (2020)](https://doi.org/10.1093/sysbio/syz069)  (see also [RevBayes tutorial on data augmentation](https://revbayes.github.io/tutorials/cont_traits/state_dependent_bm.html)). For the Birth and Death process, we decided to use a Binary State Dependent Speciation and Extinction with hidden rates (HiSSE) because it has been shown to improve parameter estimation, including the trait transition rates which are important for the trait history data augmentation. However, any kind of Birth and Death process can be use [SSE models examples](). Although the python wrapper program we developed only implements HiSSE, you can modify the RevBayes script to adapt it to your data needs. It is important notice also that our current implementation does not include Hidden states for the protracted part of the model. We stress that that could be a further development and that tests regarding the power and identifiability of the model should be performed. Our preliminary exploration suggested that the use of four states with our tree were not able to inform the model and the posterior was identical to the prior, and bigger tress with may be necessary for accurate infer parameters when a hidden (or more than 2 observed states) are used.  

To estimate the speciation completion rate, we used constraints informed by the population to species map (species matrix file). We know that speciation completion did not happened on branches that connect populations, so we can clamp the number of speciation events to 0. If the branches connects populations of two different species, at least one speciation event must have happened along the path connecting these two heterospecific populations. Branches that connect populations with unknown status or that connects more than two species are not used as constraints and the posterior distribution of the number of speciation completion events on those branches are estimate ([FIGURES]()). Therefore, a population tree with many populations per species must be used. Some population with unknown identity status can be used and the results of the probability of speciation on the branches connection these populations to others can be used for species delimitation (see below).  More detailed explanation of the algorithm can be seen at the [RevBayes script Tutorial](./TutorialRev).

>Note that our package is limited to the model priors and parameters provided in [Azevedo et al. (2024)](). We strongly recommend you to check the [RevBayes Tutorial](./TutorialRev/) for better understanding of the steps and customization of models to better use it with your data. The use of *traderpros* program with the argument *--just_script* produce the RevBayes script without running the analysis it self, and it can be useful for customization.

### The programs
[traderpros](./src/traderpros.py) and [traderpros_sim](./src/traderpros_sim.py) generates Rev scripts and calls [RevBayes](https://revbayes.github.io/) to run them.  [traderpros](./src/traderpros.py) is used for estimating parameters and fos species delimitation. [traderpros_sim](./src/traderpros_sim.py) is used for simulations.
  
[sum_and_plot](./src/sum_and_plot.py) is a wrapper of [customized R functions](./src/custom_R_functions/) that uses [RevGadgets](), [coda]() and [ggplot2]() to summarize and plot posterior distribution of parameters and probability pie charts of the different models tested with Reversible Jump MCMC (e.g. if a irreversible model of trait is more likely; if there is association between traits and rates, etc...).  
  
[conspecific_probs](./src/conspecific_probs.py) and [conspecific_binary](./src/conspecific_binary.py) are functions that use [DendroPy]() library for visualization of species delimitation results.  [conspecific_probs](./src/conspecific_probs.py) generates a heatmap with the probability of two tips (populations or individuals) be conspecific.  [conspecific_binary](./src/conspecific_binary.py) shows if two tips were estimated in the maximum *a posteriori* result as being conspecific or not.
  
All these dependencies must be installed and accessible from command line (they should be added to your system $PATH variable). See [Installation instructions](#installation). If you use Traderpros, please cite those dependencies.
  
Citations:
  


  

## Installation  
If you do not have R and RevBayes installed and added to your $PATH, please see instructions at [R project](https://cran.r-project.org/bin/linux/ubuntu/fullREADME.html) and [RevBayes](https://revbayes.github.io/download) ([or compile from source](https://revbayes.github.io/compile-linux)) pages. This programs were tested with R version 4.4.2 and RevBayes version 1.2.5. 

Open R and Install R packages dependencies if not done yet.

>> If RevGadgets installation fails, check [instructions on the website](https://revbayes.github.io/tutorials/intro/revgadgets.html). It might be related to the "magick" R package that depends on external software [ImageMagick](https://imagemagick.org/script/download.php). 


```R
install.packages("ggplot2")
install.packages("coda")
install.packages("devtools")
devtools::install_github("YuLab-SMU/ggtree")
devtools::install_github("cmt2/RevGadgets")

# RevGadgets requires "magick" and "ggimage"
# If installation fails it might be related to those.
#install.packages("magick")
#devtools::install_github("GuangchuangYu/ggimage")
```

We suggest you to install the Traderpros in a conda environment. If you do not have conda installed in your computer, check [miniconda](https://docs.anaconda.com/miniconda/install/#quick-command-line-install) for conda installation instructions. You can use the code below to create a environment with all the python dependencies we need.

>Note: that we are using the package versions which we tested traderpros. You can try newer versions, but if it does not work, try with the specified versions)

```bash
conda create -n traderenv python=3.13.0 
conda activate traderenv
pip install seaborn==0.13.2 numpy==2.1.3 matplotlib==3.9.2 dendropy==5.0.1
```

Clone this github project and install Traderpros.

```bash
git clone https://github.com/ghfazevedo/traderpros
cd traderpros
pip install .
```

The  installation should check for the dependencies and make sure you can access R and RevBayes from the command line. You can test the program with the commands below. It should run a quick analysis.

> Note: Sometimes the analyses will freeze in the initial steps of the MCMC chain. It usually happens when you have another terminal with RevBayes opened or if you opened RevBayes in that terminal before. If that happens, close the terminals and try again. 

```bash
traderpros  -opre Test \
            -odir Test \
            -t CicurinaData/BPPConstraint.MCC.CAH.nex \
            -trpt CicurinaData/CicTroglomorphism.nex \
            -spmt CicurinaData/CicSpeciesMatrix.txt \
            -tsph no \
            -npop 457 \
            -nhdd 2 \
            -mxsp 24 \
            -ngen 100 \
            -prsc 1 \
            -prlg 1 \
            -prtr 1 \
            -chkp 10

sum_and_plot \
     --out_dir TestFIGS \
     --in_dir  Test \
     --prefix_used_in_traderpros  Test \
     --burn  0.1 \
     --out_images_format  "both" \
     --path_to_custom_functions  src/custom_R_functions

conspecific_probs \
     -pt Test/Test.Protracted.MAP.tre \
     -od TestFIGS \
     -pre Test

conspecific_binary \
     -ant Test/Test.SpEvents.MAP.tre \
     -od TestFIGS \
     -pre Test
```
  
## Reproducing analyses from Azevedo et al.
To reproduce the analyses present in [Azevedo et al. ()](), run the following codes.

```bash
traderpros  -opre CicRun01 \
            -odir CicRun01 \
            -t CicurinaData/BPPConstraint.MCC.CAH.nex \
            -trpt CicurinaData/CicTroglomorphism.nex \
            -spmt CicurinaData/CicSpeciesMatrix.txt \
            -tsph no \
            -npop 457 \
            -nhdd 2 \
            -mxsp 24 \
            -ngen 100000 \
            -prsc 10 \
            -prlg 100 \
            -prtr 100 \
            -chkp 10 \
            --char_move_intensity 0.75 \
            -seed_number 


traderpros  -opre CicRun02 \
            -odir CicRun02 \
            -t data/BPPConstraint.MCC.CAH.nex \
            -trpt data/CicTroglomorphism.nex \
            -spmt data/CicSpeciesMatrix.txt \
            -tsph no \
            -npop 457 \
            -nhdd 2 \
            -mxsp 24 \
            -ngen 100000 \
            -prsc 10 \
            -prlg 100 \
            -prtr 100 \
            -chkp 1000 \
			--char_move_intensity 0.75 \
            -seed_number 
```

## What you need to run your data
For running Traderpros you will need:  
1. A dated tree (in nexus or newick format) with many populations of each species as tips, [like this one for example](./data/BPPConstraint.MCC.CAH.nex)
2. A trait matrix in nexus format associating each tip with the state found in that population, [like this one](./data/CicTroglomorphism.nex). Note that for now, missing data are not accepted and may cause issues.
3. A tab delimited file mapping each population to a species code, which are represented by natural numbers starting from 0. Every population belonging to the same species should have the same species code. See [this file for example](./data/CicSpeciesMatrix.txt)
4. You also should have a approximated idea of how many population probably exists in the clade (argument *-npop*). We acknowledge that this may be difficult, but we encourage to try different assumptions and see how it affects your conclusions. This should have an stronger effect on the estimative of Birth and Death parameters. If you think you have no good estimate and you think this can influence too much your results, you can try to modify the RevScript to run the analyses without the Birth and Death part of the model. 
5. You should have and expected the maximum number of sampled species in your clade (argument *-mxsp*), if you have populations with unknown species assignments. If all populations are known to belong to a species, this is the number of sampled species.

### Command reference
Type name of the program with the flag *-h*  to see all options and defaults.

#### traderpros
```
usage: traderpros [-h] [-s SEED_NUMBER] [-opre OUT_PREFIX] [-odir OUT_DIR] -t TREE_PATH -trpt TRAIT_PATH -spmt SP_MATRIX [-tsph {yes,no}] [-ptph TENSOR_PATH]
                  -npop NUM_TOTAL_POPULATIONS [-nhdd NUM_HIDDEN] [-mhdd MEAN_HIDDEN_HYPERPRIOR] [-nproc NUM_PROCESSORS] [-trpr TRANSITION_PRIOR_PARAM]
                  [-btpr BIRTH_PRIOR_PARAM] [-dtpr DEATH_PRIOR_PARAM] [-sppr SPECIATION_PRIOR_PARAM] -mxsp MAX_NUM_SPECIES [-mvit CHAR_MOVE_INTENSITY] [-ngen N_GEN]
                  [-prsc PRINT_SCREEN] [-prlg PRINT_LOG] [-prtr PRINT_TREE] [-nrun N_RUNS] [-chkp CHECK_POINT_INTERVAL] [-burn BURN_IN] [-justscript JUST_SCRIPT]

Estimate parameters of a trait dependent protracted speciation model.

options:
  -h, --help            show this help message and exit
  -s, --seed_number SEED_NUMBER
                        Seed number for replication (default: 0).
  -opre, --out_prefix OUT_PREFIX
                        Output file prefix (default: "Traderpros").
  -odir, --out_dir OUT_DIR
                        Output directory (default: "trader_out").
  -t, --tree_path TREE_PATH
                        Tree file path (mandatory).
  -trpt, --trait_path TRAIT_PATH
                        Trait data file path (mandatory).
  -spmt, --sp_matrix SP_MATRIX
                        Species matrix file path (mandatory).
  -tsph, --use_tensor {yes,no}
                        Specify "yes" to use TensorPhylo or "no" to not use it (default: "no"). TensorPhylo is a library to speed up likelihood computation for SSE models. It need
                        to be downloaded and installed separatly. See details on TensorPhylo at https://bitbucket.org/mrmay/tensorphylo/. Note: As for 20 Nov 2024 - Sometimes it
                        throws an error and interrupts analysis.
  -ptph, --tensor_path TENSOR_PATH
                        TensorPhylo plugin path (required if use_tensor=True).
  -npop, --num_total_populations NUM_TOTAL_POPULATIONS
                        Total number of populations (mandatory).
  -nhdd, --num_hidden NUM_HIDDEN
                        Number of hidden states (default: 2).
  -mhdd, --mean_hidden_hyperprior MEAN_HIDDEN_HYPERPRIOR
                        Mean hidden hyperprior (default: 0.587405).
  -nproc, --num_processors NUM_PROCESSORS
                        Number of processors (default: 4). only used if using TensorPhylo plugin.
  -trpr, --transition_prior_param TRANSITION_PRIOR_PARAM
                        Transition prior parameter (default: 10).
  -btpr, --birth_prior_param BIRTH_PRIOR_PARAM
                        Birth prior parameter (default: 1).
  -dtpr, --death_prior_param DEATH_PRIOR_PARAM
                        Death prior parameter (default: 1).
  -sppr, --speciation_prior_param SPECIATION_PRIOR_PARAM
                        Speciation prior parameter (default: 1).
  -mxsp, --max_num_species MAX_NUM_SPECIES
                        Max number of species (mandatory).
  -mvit, --char_move_intensity CHAR_MOVE_INTENSITY
                        Character move intensity (default: 0.5).
  -ngen, --n_gen N_GEN  Number of generations (default: 10000).
  -prsc, --print_screen PRINT_SCREEN
                        Print screen interval (default: 100).
  -prlg, --print_log PRINT_LOG
                        Print log interval (default: 100).
  -prtr, --print_tree PRINT_TREE
                        Print tree interval (default: 100).
  -nrun, --n_runs N_RUNS
                        Number of parallel runs (default: 1).
  -chkp, --check_point_interval CHECK_POINT_INTERVAL
                        Check point interval (default: 1000).
  -burn, --burn_in BURN_IN
                        Burn-in percentage (default: 0.10).
  -justscript, --just_script JUST_SCRIPT
                        Create script only, do not run (default: False).
```
#### sum_and_plot
```
usage: sum_and_plot [-h] [-od OUT_DIR] -in IN_DIR -pre PREFIX_USED_IN_TRADERPROS [-b BURN] [--color_completion COLOR_COMPLETION]
                    [--color_n_speciation COLOR_N_SPECIATION] [--color_birth COLOR_BIRTH] [--color_death COLOR_DEATH] [--color_net COLOR_NET]
                    [--color_transition COLOR_TRANSITION] [--anc_state_labels ANC_STATE_LABELS] [--anc_state_colors ANC_STATE_COLORS]
                    [--path_to_custom_functions PATH_TO_CUSTOM_FUNCTIONS] [--out_images_format {pdf,png,both}]

Python wrapper for sum_and_plot_traderpros_result in R.

options:
  -h, --help            show this help message and exit
  -od, --out_dir OUT_DIR
                        Output directory
  -in, --in_dir IN_DIR  Input directory
  -pre, --prefix_used_in_traderpros PREFIX_USED_IN_TRADERPROS
                        Prefix used in TraderPros output
  -b, --burn BURN       Burn-in percentage
  --color_completion COLOR_COMPLETION
                        Comma-separated colors for completion rate
  --color_n_speciation COLOR_N_SPECIATION
                        Comma-separated colors for number of speciation events
  --color_birth COLOR_BIRTH
                        Comma-separated colors for birth rates
  --color_death COLOR_DEATH
                        Comma-separated colors for death rates
  --color_net COLOR_NET
                        Comma-separated colors for net diversification
  --color_transition COLOR_TRANSITION
                        Comma-separated colors for transitions
  --anc_state_labels ANC_STATE_LABELS
                        Comma-separated ancestral state labels
  --anc_state_colors ANC_STATE_COLORS
                        Comma-separated colors for ancestral states
  --path_to_custom_functions PATH_TO_CUSTOM_FUNCTIONS
                        Path to custom functions
  --out_images_format {pdf,png,both}
                        Output image format
```
#### conspecific_binary
```
usage: conspecific_binary [-h] -spt SPECIATION_EVENTS_TREE [-od OUT_DIR] [-pre PREFIX] [-tax KEEP_ONLY_TAXA] [-fgs FIGSIZE] [-fts FONTSIZE]

Create plot of conspecificity.

options:
  -h, --help            show this help message and exit
  -spt, --speciation_events_tree SPECIATION_EVENTS_TREE
                        Path to tree with speciation events (Nexus format)
  -od, --out_dir OUT_DIR
                        Output directory
  -pre, --prefix PREFIX
                        Prefix for output files
  -tax, --keep_only_taxa KEEP_ONLY_TAXA
                        Comma-separated list of focal taxa to plot. At least two should be provided
  -fgs, --figsize FIGSIZE
                        Figure size as width,height (default: '12,10')
  -fts, --fontsize FONTSIZE
                        Font size for labels and ticks (default: 10)
```
#### conspecific_probs
```
usage: conspecific_probs [-h] -pt PROTRACTED_TREE [-od OUT_DIR] [-pre PREFIX] [-tax KEEP_ONLY_TAXA] [-fgs FIGSIZE] [-fts FONTSIZE]

Calculate conspecific probabilities and create heatmap

options:
  -h, --help            show this help message and exit
  -pt, --protracted_tree PROTRACTED_TREE
                        Path to protracted tree (Nexus format)
  -od, --out_dir OUT_DIR
                        Output directory
  -pre, --prefix PREFIX
                        Prefix for output files
  -tax, --keep_only_taxa KEEP_ONLY_TAXA
                        Comma-separated list of focal taxa to plot. At least two should be provided
  -fgs, --figsize FIGSIZE
                        Figure size as width,height (default: '12,10')
  -fts, --fontsize FONTSIZE
                        Font size for labels and ticks (default: 10)
```

#### traderpros_sim
```
usage: traderpros_sim [-h] [-s SEED_NUMBER] [-opre OUT_PREFIX] [-odir OUT_DIR] -bo BIRTH_OBSERVED -eo EXTINCTION_OBSERVED -bh BIRTH_HIDDEN -eh EXTINCTION_HIDDEN
                      -tr TRANSITION_RATES -ht HIDDEN_TRANS -rp ROOT_PROBS -ra ROOT_AGE [-ps POPSAMPLING_PROB] -nt N_TIPS -sc SP_COMPLETION_RATES [-ns N_SIMULATIONS]
                      [-justscript] [-org] [-unk N_UNKNOWN_SP]

Simulates data under the Traderpros model.

options:
  -h, --help            show this help message and exit
  -s, --seed_number SEED_NUMBER
                        Seed number for replication (default: 0).
  -opre, --out_prefix OUT_PREFIX
                        Output file prefix (default: "Traderpros").
  -odir, --out_dir OUT_DIR
                        Output directory (default: "trader_out").
  -bo, --birth_observed BIRTH_OBSERVED
                        A comma-separated list of state specific birth rates for observed states (mandatory).
  -eo, --extinction_observed EXTINCTION_OBSERVED
                        A comma-separated list of state specific death rates for observed states (mandatory).
  -bh, --birth_hidden BIRTH_HIDDEN
                        A comma-separated list of state specific birth rates for hidden states. The mean of rates must be 1. If you wish a model with no hidden rates,
                        chose 1,1 (mandatory).
  -eh, --extinction_hidden EXTINCTION_HIDDEN
                        A comma-separated list of state specific death rates for hidden states. The mean of rates must be 1. If you wish a model with no hidden rates,
                        chose 1,1 (mandatory).
  -tr, --transition_rates TRANSITION_RATES
                        A comma-separated list of transition rates between observed states. It should be absolute rates. First value represents transition from 0 to 1
                        and second value should be the 1 to 0 transition (mandatory).
  -ht, --hidden_trans HIDDEN_TRANS
                        A comma-separated list of transition rates between hidden states. It should be absolute rates. First value represents transition from A to B
                        and second value should be the B to 1 transition (mandatory).
  -rp, --root_probs ROOT_PROBS
                        A comma-separated list of probabilities for the all 4 root states (Observed 0 Hidden A, Observed 1 Hidden A, Observed 0 Hidden B, Observed 1
                        Hidden B). It could be values that sum to 1, or relative weights for a simplex. E.g. -rp 1,0,1,0 will translate to 0.5, 0.0, 0.5, 0.0
                        (mandatory).
  -ra, --root_age ROOT_AGE
                        The age of the root of the tree (mandatory).
  -ps, --popsampling_prob POPSAMPLING_PROB
                        Probability of sampling a population (default: 1.0).
  -nt, --n_tips N_TIPS  The number of tips in the final tree (mandatory).
  -sc, --sp_completion_rates SP_COMPLETION_RATES
                        A comma-separated list of state specific speciation completion rates for observed states (mandatory).
  -ns, --n_simulations N_SIMULATIONS
                        The number of simulations (default: 1).
  -justscript, --just_script
                        Create script only, do not run (default: False).
  -org, --organize_folder
                        Organize files into folders by simulation. (default: False)
  -unk, --n_unknown_sp N_UNKNOWN_SP
                        Number of tips to be set as "unknown" species assignment in the Species Matrix file for testing application of the model for species
                        delimitation (default: 0)
```

##### Simulating data as in Azevedo et al.
```
traderpros_sim \
  -opre sim4  \
  -odir dataSim4 \
  -bo 2.239,4.314 \
  -eo 0.341,0.324 \
  -bh 0.331,1.669 \
  -eh 0.743,1.257 \
  -tr 0.247,0.00465 \
  -ht 0.71,0.71 \
  -rp 0.444,0.294,0.0924,0.17 \
  -ra 10 \
  -ps 0.27 \
  -nt 124 \
  -sc 0.7,0.5   \
  -ns 5 \
  -org \
  -unk 3
```


# References
Sukumaran J, Holder MT, Knowles LL (2021) Incorporating the speciation process into species delimitation. [PLOS Computational Biology 17(5): e1008924](ttps://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008924) 
