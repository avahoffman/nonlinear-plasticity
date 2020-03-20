###########################################################################################
# UTILITY FUNCTIONS
###########################################################################################


run_mcmc <-
  function(comp,
           infile,
           response,
           iter = 10000) {
    # This function runs the MCMC sampler and makes a fit summary
    
    # Read in data
    df <-
      read.csv(infile, header = T)
    
    # Filter out NA response variable
    response.1 <-
      df[, c(response)]
    df.filt <-
      df[!(response.1 = is.na(response.1)),]
    
    # Make design matrix
    dm <-
      model.matrix(~ as.factor(trt) * as.factor(geno), data = df.filt) # Model Matrix
    
    # New dependent variable with NA removed
    response.var <-
      df.filt[, c(response)]
    
    # Declare model components
    model.components <-
      list(
        'N' = nrow(df.filt), # 
        'y' = response.var,
        'X' = dm,
        'J' = ncol(dm)
      )
    
    ##SAMPLE
    fit <-
      sampling(
        comp,
        data = model.components,
        iter = iter,
        warmup = iter / 2,
        thin = 1,
        chains = 3
      )
    summ_fit <-
      summary(fit)
    
    return(as.data.frame(summ_fit$summary))
  }


gather_posterior_data <-
  function(summ_fit,
           response) {
    ests <-
      # Collect mean, 25-75, CI, and Rhat
      (summ_fit[grep("b", rownames(summ_fit)), c(1, 5, 6, 4, 8, 10)])
    descr <- gsub("as.factor", "", colnames(dm))
    descr <- gsub("geno)3", "geno)5", descr)
    ests <- cbind(rep(response, nrow(ests)), descr, ests)
    
    # Calculate Pr
    Pr_calc <- data.frame()
    Pr_vals <-  as.data.frame(ests[, 6:7])
    for (r in 1:nrow(Pr_vals)) {
      max_ <- max(Pr_vals[r,])
      min_ <- min(Pr_vals[r,])
      # Calculate the proportion of overlap on zero
      if (max_ < 0 & min_ < 0) {
        pr. <- 1
      } else if (max_ > 0 & min_ > 0) {
        pr. <- 1
      } else if (abs(max_) > abs(min_)) {
        pr. <- abs(max_) / (abs(min_) + abs(max_))
      } else if (abs(max_) < abs(min_)) {
        pr. <- abs(min_) / (abs(min_) + abs(max_))
      } else {
        pr. <- "NA"
      }
      Pr_calc[r, 1] <- round(pr., 2)
    }
    ests <- cbind(ests, Pr_calc)
    colnames(ests)[9] <- "Pr"
    colnames(ests)[1] <- "measure"
    
    return(ests)
  }


make_normality_plots <-
  function(summ_fit,
           response) {
    sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
    rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
    pdf(
      file = paste(
        "figures/normality_test_plots/",
        response,
        "_normtest.pdf",
        sep = ""
      ),
      height = 7,
      width = 7
    )
    test_norm(rdf$mean)
    dev.off()
  }