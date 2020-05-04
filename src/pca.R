###########################################################################################
## PRINCIPAL COMPONENTS WERE USED TO DETERMINE TRAITS OF INTEREST
###########################################################################################

getprcomps <-
  function(df,
           limits) {
    # This function runs PCA on subsets of data
    
    df <- na.omit(df)
    df.responsevars <- df[, limits]
    means <- apply(as.matrix(df.responsevars), 2, mean)
    normalized.df <-
      matrix(nrow = nrow(df.responsevars),
             ncol = ncol(df.responsevars))
    for (i in 1:ncol(df.responsevars)) {
      ## RESPONSES MUST BE SCALED
      coldata <- df.responsevars[, i] / means[i]
      normalized.df[, i] <- coldata
    }
    colnames(normalized.df) <- names(df.responsevars)
    pr.cmp2 <- prcomp(normalized.df)
    pr.cmp <- princomp(normalized.df)
    print(summary(pr.cmp2))
    return(pr.cmp2)
  }


produce_prcomps <-
  function() {
    # Wrapper function to run PCA for all data subsets
    
    write.csv(getprcomps(
      df = read.csv(file = "data/biomass_plants.csv", header = T),
      limits = c(8:23) 
    )$rotation[, 1:3], # Select response vars to use
    file = "output/pca/biomass_plants_pca_rotation.csv")
    
    write.csv(getprcomps(
      df = read.csv(file = "data/phys_plants_clean.csv", header = T),
      limits = c(54:63)
    )$rotation[, 1:3], # Select response vars to use
    file = "output/pca/phys_plants_pca_rotation.csv")
    
    write.csv(getprcomps(
      df = read.csv(file = "data/all_plants_clean.csv", header = T),
      limits = c(27:30)
    )$rotation[, 1:3], # Select response vars to use
    file = "output/pca/all_plants_pca_rotation.csv")
    
    write.csv(getprcomps(
      df = read.csv(file = "data/recovery_plants_clean.csv", header = T),
      limits = c(12, 13, 15, 16, 18:21, 32:35) # Select response vars to use
    )$rotation[, 1:3], #excluding flowering params used in later models
    file = "output/pca/recovery_plants_pca_rotation.csv")
  }
