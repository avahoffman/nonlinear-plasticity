###########################################################################################
## PRINCIPAL COMPONENTS WERE USED TO DETERMINE TRAITS OF INTEREST
###########################################################################################

getprcomps <-
  function(df,
           limits = c(8:ncol(df))) {
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
    
    getprcomps(
      df = read.csv(file = "data/biomass_plants.csv", header = T),
      limits = c(8:23)
    )
    getprcomps(df = read.csv(file = "data/phys_plants.csv", header = T),
               limits = c(10:59))
    getprcomps(df = read.csv(file = "data/all_plants.csv", header = T),
               limits = c(6:25))
    getprcomps(
      df = read.csv(file = "data/recovery_plants.csv", header = T),
      limits = c(10:12, 14, 15, 17:26)
    ) #excluding flowering params used in later models
  }
