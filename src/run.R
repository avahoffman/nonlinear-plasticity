###########################################################################################
# EXECUTION SCRIPT FOR THE PROJECT PIPELINE
###########################################################################################

# Set working directory for the repository (should be the git repo):
wd <-
  "/Users/avahoffman/Dropbox/Research/Andropogon_geno_drought_study_2014/nonlinear-plasticity"

setwd(wd)
sessionInfo()

# General functions and configuration
source("src/utils.R")

# Specific functions
source("src/theoretical_fig.R")
source("src/data_clean.R")
source("src/corr_analysis.R")
source("src/pca.R")
source("src/measure_models.R")

###########################################################################################

# Make the theoretical figure
make_theor_fig(outfile = "figures/theoretical_fig.pdf")

# Clean data
clean_biomass_data()
clean_phys_data()
clean_all_plant_data()
clean_recovery_data()

# Run correlations to see if perhaps we should select only a few key traits
# TODO: Rerun these before publication / make sure indices are right in the corr_analysis script.
run_corr_tests()

# Run principal components analysis to determine traits of interest
produce_prcomps()

# Modeling ----
do_measure_mcmc_sampling()
