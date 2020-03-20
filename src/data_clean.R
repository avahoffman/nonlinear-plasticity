###########################################################################################
## REMOVE OUTLIERS AND MAKE NECESSARY CALCULATIONS
###########################################################################################
library(outliers)
library(dplyr)


clean_biomass_data <-
  function() {
    infile <- "biomass_plants"
    
    df <-
      read.csv(file = paste("data/", infile, ".csv", sep = ""),
               header = T)
    
    for (measure in 9:23) {
      d <- df[, measure]
      d[outliers::scores(na.omit(d), prob = 0.99) == 1] <-
        mean(d[outliers::scores(na.omit(d), prob = 0.99) == 0])
      df[, measure] <- d
    }
    
    write.csv(df, file = paste("data/", infile, "_clean.csv", sep = ""))
  }    
    
    
clean_phys_data <-
  function() {
    infile <- "phys_plants"
    
    df <-
      read.csv(file = paste("data/", infile, ".csv", sep = ""),
               header = T)
    
    for (measure in 11:45) {
      d <- df[, measure]
      d[outliers::scores(na.omit(d), prob = 0.98) == 1] <-
        mean(d[outliers::scores(na.omit(d), prob = 0.98) == 0])
      df[, measure] <- d
    }
    
    # Anet and gs should never be less than zero
    df$Anet5_0907[df$Anet5_0907 < 0] <- NA
    df$Anet6_0913[df$Anet6_0913 < 0] <- NA
    df$Anet7_0920[df$Anet7_0920 < 0] <- NA
    df$Anet8_0925[df$Anet8_0925 < 0] <- NA
    df$Anet9_1004[df$Anet9_1004 < 0] <- NA
    df$Anet10_1011[df$Anet10_1011 < 0] <- NA
    df$gs5_0907[df$gs5_0907 < 0] <- NA
    df$gs6_0913[df$gs6_0913 < 0] <- NA
    df$gs7_0920[df$gs7_0920 < 0] <- NA
    df$gs8_0925[df$gs8_0925 < 0] <- NA
    df$gs9_1004[df$gs9_1004 < 0] <- NA
    df$gs10_1011[df$gs10_1011 < 0] <- NA
    df$WUEi5_0907[df$WUEi5_0907 < 0] <- NA
    df$WUEi6_0913[df$WUEi6_0913 < 0] <- NA
    df$WUEi7_0920[df$WUEi7_0920 < 0] <- NA
    df$WUEi8_0925[df$WUEi8_0925 < 0] <- NA
    df$WUEi9_1004[df$WUEi9_1004 < 0] <- NA
    df$WUEi10_1011[df$WUEi10_1011 < 0] <- NA
    
    # Create new variables
    
    for (i in 1:nrow(df)) {
      # Measure means
      df$uAnet[i] <-
        mean(
          c(
            df$Anet5_0907[i],
            df$Anet6_0913[i],
            df$Anet7_0920[i],
            df$Anet8_0925[i],
            df$Anet9_1004[i],
            df$Anet10_1011[i]
          ),
          na.rm = T
        )
      df$ugs[i] <-
        mean(c(
          df$gs5_0907[i],
          df$gs6_0913[i],
          df$gs7_0920[i],
          df$gs8_0925[i],
          df$gs9_1004[i],
          df$gs10_1011[i]
        ),
        na.rm = T)
      df$ufv[i] <-
        mean(c(
          df$fv5_0907[i],
          df$fv6_0913[i],
          df$fv7_0920[i],
          df$fv8_0925[i],
          df$fv9_1004[i],
          df$fv10_1011[i]
        ),
        na.rm = T)
      df$uE[i] <-
        mean(c(
          df$E5_0907[i],
          df$E6_0913[i],
          df$E7_0920[i],
          df$E8_0925[i],
          df$E9_1004[i],
          df$E10_1011[i]
        ),
        na.rm = T)
      df$uWUEi[i] <-
        mean(c(
          df$WUEi5_0907[i],
          df$WUEi6_0913[i],
          df$WUEi7_0920[i],
          df$WUEi8_0925[i],
          df$WUEi9_1004[i],
          df$WUEi10_1011[i]
        ),
        na.rm = T)
      
      # measure max
      df$max_Anet[i] <-
        max(
          c(
            df$Anet5_0907[i],
            df$Anet6_0913[i],
            df$Anet7_0920[i],
            df$Anet8_0925[i],
            df$Anet9_1004[i],
            df$Anet10_1011[i]
          ),
          na.rm = T
        )
      df$max_gs[i] <-
        max(c(
          df$gs5_0907[i],
          df$gs6_0913[i],
          df$gs7_0920[i],
          df$gs8_0925[i],
          df$gs9_1004[i],
          df$gs10_1011[i]
        ),
        na.rm = T)
      df$max_fv[i] <-
        max(c(
          df$fv5_0907[i],
          df$fv6_0913[i],
          df$fv7_0920[i],
          df$fv8_0925[i],
          df$fv9_1004[i],
          df$fv10_1011[i]
        ),
        na.rm = T)
      df$max_E[i] <-
        max(c(
          df$E5_0907[i],
          df$E6_0913[i],
          df$E7_0920[i],
          df$E8_0925[i],
          df$E9_1004[i],
          df$E10_1011[i]
        ),
        na.rm = T)
      df$max_WUEi[i] <-
        max(c(
          df$WUEi5_0907[i],
          df$WUEi6_0913[i],
          df$WUEi7_0920[i],
          df$WUEi8_0925[i],
          df$WUEi9_1004[i],
          df$WUEi10_1011[i]
        ),
        na.rm = T)
    }
    
    write.csv(df, file = paste("data/", infile, "_clean.csv", sep = ""))
    
  }


clean_all_plant_data <-
  function() {
    infile <- "all_plants"
    
    df <-
      read.csv(file = paste("data/", infile, ".csv", sep = ""),
               header = T)
    
    for (measure in 7:25) {
      d <- df[, measure]
      d[outliers::scores(na.omit(d), prob = 0.999) == 1] <-
        mean(d[outliers::scores(na.omit(d), prob = 0.999) == 0])
      df[, measure] <- d
    }
    
    for (i in 1:nrow(df)) {
      # Measure means
      df$uH[i] <-
        mean(
          c(
            df$H2_0821[i],
            df$H4_0904[i],
            df$H6_0918[i],
            df$H7_0926[i],
            df$H8_1002[i],
            df$H9_1009[i],
            df$H10_1015[i]
          ),
          na.rm = T
        )
      df$urgr[i] <-
        mean(
          c(
            df$rgr_0904[i],
            df$rgr_0918[i],
            df$rgr_0926[i],
            df$rgr_1002[i],
            df$rgr_1009[i],
            df$rgr_1015[i]
          ),
          na.rm = T
        )
      
      # Measure maximums
      df$max_H[i] <-
        max(
          c(
            df$H2_0821[i],
            df$H4_0904[i],
            df$H6_0918[i],
            df$H7_0926[i],
            df$H8_1002[i],
            df$H9_1009[i],
            df$H10_1015[i]
          ),
          na.rm = T
        )
      df$max_rgr[i] <-
        max(
          c(
            df$rgr_0904[i],
            df$rgr_0918[i],
            df$rgr_0926[i],
            df$rgr_1002[i],
            df$rgr_1009[i],
            df$rgr_1015[i]
          ),
          na.rm = T
        )
    }
    
    write.csv(df, file = paste("data/", infile, "_clean.csv", sep = ""))
  }   


clean_recovery_data <-
  function() {
    infile <- "recovery_plants"
    
    df <-
      read.csv(file = paste("data/", infile, ".csv", sep = ""),
               header = T)
    
    for (measure in 11:30) {
      d <- df[, measure]
      d[outliers::scores(na.omit(d), prob = 0.999) == 1] <-
        mean(d[outliers::scores(na.omit(d), prob = 0.999) == 0])
      df[, measure] <- d
    }
    
    for (i in 1:nrow(df)) {
      # Measure means
      df$uH[i] <-
        mean(
          c(
            df$H11_1023[i],
            df$H12_1104[i],
            df$H13_1111[i],
            df$H14_1117[i],
            df$H15_1124[i]
          ),
          na.rm = T
        )
      df$urgr[i] <-
        mean(
          c(
            df$rgr11_1023[i],
            df$rgr12_1104[i],
            df$rgr13_1111[i],
            df$rgr14_1117[i],
            df$rgr15_1124[i]
          ),
          na.rm = T
        )
      
      # Measure max
      df$max_H[i] <-
        max(
          c(
            df$H11_1023[i],
            df$H12_1104[i],
            df$H13_1111[i],
            df$H14_1117[i],
            df$H15_1124[i]
          ),
          na.rm = T
        )
      df$max_rgr[i] <-
        max(
          c(
            df$rgr11_1023[i],
            df$rgr12_1104[i],
            df$rgr13_1111[i],
            df$rgr14_1117[i],
            df$rgr15_1124[i]
          ),
          na.rm = T
        )
    }
    
    write.csv(df, file = paste("data/", infile, "_clean.csv", sep = ""))
  }


clean_biomass_data()
clean_phys_data()
clean_all_plant_data()
clean_recovery_data()
