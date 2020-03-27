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
source("src/timeline.R")
source("src/data_clean.R")
source("src/corr_analysis.R")
source("src/pca.R")
source("src/measure_models.R")
source("src/theta_model.R")
source("src/measure_plot.R")
source("src/breakpoint_analysis.R")
source("src/phenotype_plot.R")


###########################################################################################

# Make the theoretical figure
make_theor_fig(outfile = "figures/theoretical_fig.pdf")

# Make experimental timeline
make_timeline()

# Clean data
clean_biomass_data()
clean_phys_data()
clean_all_plant_data()
clean_recovery_data()

# Run correlations to see if perhaps we should select only a few key traits
run_corr_tests()

# Run principal components analysis to determine traits of interest
produce_prcomps()

# Modeling ----
do_measure_mcmc_sampling() # Measures
do_flwr_rh_mcmc_sampling() # Thetas

# Breakpoint analysis ----
run_breakpoint_analysis()

# Plot main effects ----
make_fig_effects(outfile = "figures/fig_effects_and_plasticity_functions.pdf")

# Plot phenotypes by treatment ====
cycle_phenotype_plots()