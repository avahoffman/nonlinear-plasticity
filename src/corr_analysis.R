###########################################################################################
## BEFORE EXAMINING INDIVIDUAL PHENOTYPIC MEASURES, WE WANT TO DETERMINE WHETHER TRAITS
## ARE CORRELATED - IF THEY ARE WE PROBABLY WANT TO DO A SUBSET
###########################################################################################


run_corr_tests <-
  function() {
    # This function plots and writes all correlations for subsets of traits
    # to determine whether we should do just a few traits or not
    
    df <- read.csv("data/biomass_plants.csv", header = T)
    df <- na.omit(df)
    pdf(file = "figures/correlations/cor_biomass_plants.pdf",
        height = 8,
        width = 8)
    plot(df[, 8:ncol(df)])
    dev.off()
    write.csv(cor(df[, 8:ncol(df)]), file = "output/correlations/cor_biomass_plants.csv")
    setwd(wd)
    
    df <- read.csv("data/phys_plants.csv", header = T)
    df <- na.omit(df)
    pdf(file = "figures/correlations/cor_phys_plants.pdf",
        height = 30,
        width = 30)
    plot(df[, 10:ncol(df)])
    dev.off()
    write.csv(cor(df[, 10:ncol(df)]), file = "output/correlations/cor_phys_plants.csv")
    setwd(wd)
    
    df <- read.csv("data/all_plants.csv", header = T)
    df <- na.omit(df)
    pdf(file = "figures/correlations/cor_growth_plants.pdf",
        height = 15,
        width = 15)
    plot(df[, 6:ncol(df)])
    dev.off()
    write.csv(cor(df[, 6:ncol(df)]), file = "output/correlations/cor_all_plants.csv")
    setwd(wd)
    
    df <- read.csv("data/recovery_plants.csv", header = T)
    df <- na.omit(df)
    pdf(file = "figures/correlations/cor_recovery_plants.pdf",
        height = 8,
        width = 8)
    plot(df[, 10:ncol(df)])
    dev.off()
    write.csv(cor(df[, 10:ncol(df)]), file = "output/correlations/cor_recovery_plants.csv")
    setwd(wd)
  }