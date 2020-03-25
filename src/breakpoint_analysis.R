###########################################################################################
## DETERMINE BREAKPOINTS IN GENOTYPE TREATMENT FUNCTION
###########################################################################################
# Very nice information about segmented here: https://rpubs.com/MarkusLoew/12164

get_breakpoints <- 
  function(df, response){
      fit <- lm(as.formula(paste(response," ~ trt", sep = "")), data = df)
      seg_fit <- suppressWarnings(segmented::segmented(fit, 
                           seg.Z = ~ trt))
      if (is.null(seg_fit$psi)){
        return(c(0, 0, 0))
      } else {
        breakpoint <- seg_fit$psi
        return(breakpoint)
      }
  }

cycle_genotypes <- 
  function(infile, response){
    # This function ..
    # 
    # infile: subpath of input file
    # response: string or vector of strings indicating traits on which to perform 
    # breakpoint analysis
    for (resp in response){
          for (i in 1:3){
            df <-
              read.csv(infile, header = T)
            df_one_geno <- df[(df$geno == i),]
            if (i > 1){
              breakpoints <- rbind(breakpoints, get_breakpoints(df_one_geno, resp))
            } else {
              breakpoints <- get_breakpoints(df_one_geno, resp)
            }
          }
      breakpoints <- as.data.frame(breakpoints)
      breakpoints$geno <- c(11, 2, 5)
      breakpoints$measure <- c(rep(resp, 3))
      if (resp != response[1]){
        breakpoints_df <- rbind(breakpoints_df, breakpoints)
      } else {
          breakpoints_df <- breakpoints
        }
    }
    rownames(breakpoints_df) <- seq(1, nrow(breakpoints_df))
    return(breakpoints_df)
  }

cycle_genotypes(infile = "data/recovery_plants_clean.csv", response =  c("Bv", "Bf"))
