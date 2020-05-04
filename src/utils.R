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
        'N' = nrow(df.filt),
        #
        'y' = response.var,
        'J' = ncol(dm),
        'X' = dm
      )
    
    # SAMPLE
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
      as.data.frame(summary(fit)$summary)
    
    return(summ_fit)
  }


gather_posterior_data <-
  function(summ_fit,
           response,
           recovery = F) {
    ests <- as.data.frame(rbind(
      # Collect mean, 25-75, CI, and Rhat
      (summ_fit[grep("b", rownames(summ_fit)), c(1, 3, 5, 7, 4, 8, 10)]),
      (summ_fit[grep("G", rownames(summ_fit)), c(1, 3, 5, 7, 4, 8, 10)]),
      (summ_fit[grep("T", rownames(summ_fit)), c(1, 3, 5, 7, 4, 8, 10)]),
      (summ_fit[grep("I", rownames(summ_fit)), c(1, 3, 5, 7, 4, 8, 10)]),
      (summ_fit[grep("Y", rownames(summ_fit)), c(1, 3, 5, 7, 4, 8, 10)])
    ))
    ests <-
      ests[!(grepl("beta", rownames(ests))), ] # Drop model param "beta"
    descr <- gsub("as.factor", "", colnames(dm))
    descr <- c(gsub("geno)3", "geno)5", descr),
               rep("no_descr", 20))
    
    param <-
      c(
        rep("beta", 15),
        "G11[R]-G2[R]",
        "G11[R]-G5",
        "G2[R]-G5",
        "trt_effect",
        "int_effect",
        rep("posterior value", 15)
      )
    
    geno <-
      c(rep("none", 20), rep("G11", 5), rep("G2", 5), rep("G5", 5))
    
    trt <- c(rep("none", 20), rep(c("10", "15", "20", "25", "Sat'd"), 3))
    
    ests <- cbind(param, ests, geno, trt)
    
    if (!(recovery)) {
      ests <- cbind(rep(response, nrow(ests)), descr, ests)
    } else {
      ests <-
        cbind(rep(paste(response, "_recovery", sep = ""), nrow(ests)), descr, ests)
    }
    
    # Calculate Pr
    Pr_calc <- data.frame()
    Pr_vals <-  as.data.frame(ests[, 8:9])
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
    colnames(ests)[13] <- "Pr"
    colnames(ests)[1] <- "measure"
    
    return(ests)
  }


make_normality_plots <-
  function(summ_fit,
           response) {
    # This function makes the posterior checks
    
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


run_flwr_rh_mcmc <-
  function(comp,
           response,
           iter = 10000) {
    # This function does sampling for flowering or resprouting data
    
    # Read in data
    df <-
      read.csv("data/recovery_plants_clean.csv", header = T)
    response <- response
    
    # Filter out NA response variable
    response.1 <-
      df[, c(response)]
    df.filt <-
      df[!(response.1 = is.na(response.1)),]
    
    # New dependent variable with NA removed
    response.var <-
      df.filt[, c(response)]
    
    # Declare model components
    model.components <- list(
      'N' = nrow(df.filt),
      'w' = response.var,
      'geno' = df.filt$geno,
      'G' = 3
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


gather_flwr_rh_posterior_data <-
  function(summ_fit,
           response) {
    ests <- as.data.frame(
      # Collect mean, 25-75, CI, and Rhat
      (summ_fit[grep("theta", rownames(summ_fit)),
                c(1, 3, 5, 7, 4, 8, 10)]))
    descr <- c("geno 11", "geno 2", "geno 5")
    ests <- cbind(rep(response, nrow(ests)), descr, ests)
    
    colnames(ests)[1] <- "measure"
    
    return(ests)
  }


g_legend <- function(a.gplot) {
  tmp <-
    ggplot_gtable(ggplot_build(a.gplot))
  leg <-
    which(sapply(tmp$grobs,
                 function(x)
                   x$name) == "guide-box")
  legend <-
    tmp$grobs[[leg]]
  
  return(legend)
}
