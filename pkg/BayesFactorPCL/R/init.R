
BFAnovaInit<-function()
{
  tclRequire("BWidget")
} 

BFAnovaFonts=function(){
list(
fontHeading = tkfont.create(family="times",size=16,weight="bold")
)
}

BFAnovaHelp=function(){
list(
ColSelHelp = "Each analysis requires at least four columns:\n1. The response column is the column containing the response of the participant on each trial (0='same',1='change').\n2. The change column is the column containing the correct answer (0=same,1=change).\n3. The set size column contains the how many items the participant was asked to remember in the change detection task.\n4. Once you've added the three columns above, choose the columns of interest, and add them as either categorical or continuous. A categorical variable will be modelled with random effects; a continous variable will be modelled with a slope.",

EffSelHelp="The effect selection screen is where you build the second level of the hierarchical model. Using the columns of interest specified in the previous screen, WOMMBAT lists all possible effects using those columns. If you'd like to include an effect in the model, simply decide which parameter(s) (K, Z, or G) will include the effect and click 'Affects X' to add it.\n\nAt the bottom are some aditional settings:\n1. Use Zone-out parameter - Will the Zone-out parameter be estimated? For some designs, estimation of this parameter is not possible. Unclicking the box will set the Z parameter to a high value, with a very informative prior.\n2. Model on K - Different experimental designs require different models on capacity. For uncued designs, a high-threshold model is appropriate (Pashler, 1988); for cued designs, a double-high threshold is appropriate (Cowan et al., 1995).", 

PriorSetHelp="Most of the priors can be specified on the prior specification screen.\n1. The same inverse gamma prior is placed on all variances which are not modelled in a covariance matrix.\n2. A Wishart prior with an identity scale matrix is placed on all covariance matrices. The degrees of freedom may be specified here, and the same value is used for all covariance matrices.\n3. A normal prior is placed on the grand mean capacity, and the grand mean Z and G parameters. Keep in mind when specifying priors on K, Z, and G that K is on a units of capacity scale; Z and G are on the logit scale (0 -> probability .5; 2 -> probability .88). The default priors are fairly broad. They should be made more informative if necessary.\n4. A normal hierarchical prior is placed on groups of slope parameters. The prior on the variance of the slope parameters is the inverse gamma mentioned in (1). The prior on the mean is flat (noninformative). These settings cannot be changed, but should not be problematic.",

MCMCHelp="The MCMC settings screen allows you to tweak the Markov Chain Monte Carlo settings for optimal parameter estimation.\n1. General options : The number of iterations specifies the length of the MCMC chain used to obtain marginal posterior distributions. The longer, the better, but longer chains take more time. The burnin is the number of initial iterations thrown away when computing posterior means. If progress>0, then a text progress bar will let you know when percentage is done.\n2. Test run : Because it is often useful to do a small initial run to check the quality of the MCMC chain, it is wise to select 'Test run' before doing a long run. If 'Test run' is selected, after the run a summary of the chains will be shown. Following this, you will be returned to the MCMC settings screen so that you can tune your chain parameters.\n3. optim settings : Numerical optimization is used to select reasonable starting values for the MCMC chains. If the process is taking too long, it may be turned off by setting the max. optim() iterations to 0. If the starting values returned are unreasonable, the algorithm may need more iterations to converge. In this case, increase the number of iterations.\n4. Hybrid Monte Carlo settings : Hybrid Monte Carlo is used to obtain estimates of all parameters, excluding covariances and slope means. The quality of the samples depends in these settings. If chains are too highly autocorrelated, or the acceptance rate is outside of (0.6,0.9), change the HMC parameters and try a few test runs.",


InstrHelp="Click help at any time for help regarding the options in a screen.",

OutputHelp="If you would like to analyze the output in a program other than R, or save the entire analysis (with settings) for later, you may output some files to CSV. If output is enabled, a small information file is written. Additional files can be selected to be written, as shown on the output screen. If output is disabled, no files are written.",

InitInstructions="Welcome to the WOrking Memory Modelling using Bayesian Analysis Techniques (WOMMBAT) system. If you have not specified a filename to analyze, you can open a file after you click 'Next.' Currently, WOMMBAT only reads CSV files, and they must have the variable names as the first row. If in doubt, open the file with a text editor and check. Each row of the CSV file should correspond to a single trial in your experiment. Click next to select a file.\n\nThis message may be skipped by using the skipInstr=TRUE argument next time you start the GUI."
)			
}
