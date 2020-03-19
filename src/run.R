###########################################################################################
# EXECUTION SCRIPT FOR THE PROJECT PIPELINE
###########################################################################################

# Set working directory for the repository (should be the git repo):
wd <-
  "/Users/avahoffman/Dropbox/Research/Andropogon_geno_drought_study_2014/nonlinear-plasticity"

setwd(wd)
sessionInfo()

# General functions and configuration
# source("utils/data_utils.R")

# Specific functions
source("src/theoretical_fig.R")

###########################################################################################

# Make the theoretical figure
make_theor_fig(outfile = "figures/theoretical_fig.pdf")