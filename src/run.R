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
source("src/plot_theoretical.R")
source("src/plot_timeline.R")
source("src/data_clean.R")
source("src/corr_analysis.R")
source("src/pca.R")
source("src/measure_models.R")
source("src/theta_model.R")
source("src/measure_plot.R")
source("src/breakpoint_analysis.R")
source("src/phenotype_plot.R")
source("src/lda.R")


###########################################################################################

# Make the theoretical figure
make_theor_fig(outfile = "figures/theoretical_fig.pdf")

# Make experimental timeline to go in the supplementary material
make_timeline()

# Clean data
clean_biomass_data()
clean_phys_data()
clean_all_plant_data()
clean_recovery_data()

# Run correlations to check for covariance among phenotypes
# run_corr_tests()

# Run principal components analysis to determine traits of interest
# produce_prcomps()

# Modeling ----
do_measure_mcmc_sampling() # Measures
do_flwr_rh_mcmc_sampling() # Thetas

# Breakpoint analysis ----
run_breakpoint_analysis()

# Plot main effects ----
make_effect_plot(genotype_comparison = T)
ggsave(file = "figures/genotype_effects_20VWC.pdf", height = 9, width=10)
make_effect_plot()
ggsave(file = "figures/treatment_effects.pdf", height = 7, width=7)

# Plot LDAs ----
gather_lda_plots(outfile = "figures/genotype_LDAs.pdf")

# Breakpoint heatmap ----
make_breakpoint_plot()
ggsave(file = "figures/breakpoints.pdf", height = 7, width = 7)

# Recovery effects
make_effect_plot(recovery = T)
ggsave(file = "figures/recovery_trt.pdf", height = 4, width=7)

# Plot phenotypes by treatment ====
cycle_phenotype_plots()
growth_summary()
instantaneous_summary()
cumulative_summary()
recovery_summary()