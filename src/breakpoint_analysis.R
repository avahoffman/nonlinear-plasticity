###########################################################################################
## DETERMINE BREAKPOINTS IN GENOTYPE TREATMENT FUNCTION
###########################################################################################
# Very nice information about segmented here: https://rpubs.com/MarkusLoew/12164
library(segmented)


get_breakpoints <-
  function(df, response) {
    # This function detects breakpoints in the treatment regression using the segmented
    # package
    # df: data frame from which response variable is collected
    # response: response variable (string)
    
    # Fit glm
    fit <-
      glm(as.formula(paste(response, " ~ trt", sep = "")), data = df)
    
    # Find breakpoints
    seg_fit <-
      suppressWarnings(segmented(
        fit,
        seg.Z = ~ trt,
        # Starting point is 22 %VWC according to Hoover et al. 2014 this is a reasonable assumption
        psi = 20,
        # Only one breakpoint, for simplicity
        npsi = 1,
        # Bootstrap params
        control = seg.control(n.boot = 50, seed = 2)
      ))
    
    plot.segmented(seg_fit)
    lines.segmented(seg_fit)
  
    # 95% CI lower
    lower <- confint.segmented(seg_fit)[2]
    
    # 95% CI upper
    upper <- confint.segmented(seg_fit)[3]
    
    # If no breakpoint detected, return zeroes
    if (is.null(seg_fit$psi)) {
      return(c(0, 0, 0, 0, 0))
    } else {
      breakpoint <- c(seg_fit$psi, # Estimate and St.Error
                      lower, # 95% CI lower
                      upper # 95% CI upper
      )
      return(breakpoint)
    }
  }

cycle_genotypes <-
  function(infile, response, subset) {
    # This function ..
    #
    # infile: subpath of input file
    # response: string or vector of strings indicating traits on which to perform
    # breakpoint analysis
    
    # Iterate through responses
    for (resp in response) {
      # Iterate through genotypes
      for (i in 1:3) {
        df <-
          read.csv(infile, header = T)
        # Remove Saturated treatment since we don't know the exact %VWC
        df_one_geno <- df[(df$geno == i & df$trt != 30),]
        if (i > 1) {
          breakpoints <-
            rbind(breakpoints, get_breakpoints(df_one_geno, resp))
        } else {
          breakpoints <- get_breakpoints(df_one_geno, resp)
        }
      }
      
      breakpoints <- as.data.frame(breakpoints)
      # Needed in case a 3x3 matrix of zeros with no colnames is made!
      # 
      colnames(breakpoints) <-
        c("Initial", "Est.", "St.Err", "lower", "upper")
      
      # Add genotype indicator
      breakpoints$geno <- c(11, 2, 5)
      
      # Add columns for the response variable and plant growth subset
      breakpoints$measure <- c(rep(resp, 3))
      breakpoints$subset <- c(rep(subset, 3))
      
      # Combine everything
      if (resp != response[1]) {
        breakpoints_df <- rbind(breakpoints_df, breakpoints)
      } else {
        breakpoints_df <- breakpoints
      }
    }
    
    # If error is higher than 3, it's probably not a real breakpoint
    # breakpoints_df[(breakpoints_df$St.Err > 3), ]$Est. <- 0
    
    # Rename row names
    rownames(breakpoints_df) <- seq(1, nrow(breakpoints_df))
    
    return(breakpoints_df)
  }


run_breakpoint_analysis <-
  function() {
    # This function runs breakpoint analysis on subsets of responses and
    # writes them to file for future use in plotting, etc.
    
    df <-
      rbind(
        
        # Recovery subset
        cycle_genotypes(
          infile = "data/recovery_plants_clean.csv",
          subset = "Recovery",
          response =  c(
            "urgr",
            "uH",
            "max_rgr",
            "max_H",
            "B_total",
            "B_below",
            "B_above",
            "B_A"
          )
        ),
        
        # Cumulative subset
        cycle_genotypes(
          infile = "data/biomass_plants_clean.csv",
          subset = "Cumulative",
          response =  c(
            "ssa",
            "srl",
            "SLA",
            "Rv",
            "Rt",
            "Rsa_la",
            "Rsa",
            "Rl",
            "Rd",
            "LA",
            "DMCv",
            "DMCr",
            "Bv",
            "Br",
            "B_A"
          )
        ),
        
        # Instantaneous subset
        cycle_genotypes(
          infile = "data/phys_plants_clean.csv",
          subset = "Instantaneous",
          response =  c(
            "uWUEi",
            "ugs",
            "ufv",
            "uAnet",
            "max_WUEi",
            "max_gs",
            "max_fv",
            "max_Anet"
          )
        ),
        
        # Growth subset
        cycle_genotypes(
          infile = "data/all_plants_clean.csv",
          subset = "Growth",
          response =  c("urgr",
                        "uH",
                        "max_rgr",
                        "max_H")
        )
      )
    
    write.csv(df, file = "output/breakpoint_analysis.csv")
  }
