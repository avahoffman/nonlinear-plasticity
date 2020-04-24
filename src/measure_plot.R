###########################################################################################
## CODE FOR PLOTTING
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)

# Colors ----
strong_effect_color <- "#d0d1e6"
weak_effect_color <- "#fdd49e"
no_effect_color <- "#fc8d59"
colfunc <- colorRampPalette(c('#fff7fb', "#a6bddb"))
breakpoint_pal <- c(colfunc(14), '#fdd49e')
zero_line_col <- "#bdbdbd"


clean_posterior_data_for_plotting <-
  function(param_) {
    # This function gathers data from the posterior distributions from the modelling step.
    #
    # param_: string indicating the type of parameter to be filtered out of the final
    # posterior distribution data, eg. "trt_effect"
    
    # Data containing full / parse-able names for phenotypic measures
    orders <-
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>%
      dplyr::mutate(measure = as.character(measure))
    
    df <-
      # Read posterior data
      read.csv("output/posterior_output.csv", header = T) %>%
      
      # Filter by parameter
      dplyr::filter(param == param_) %>%
      
      # Convert several columns to numeric values so that boxplots can be created
      dplyr::mutate(`X2.50.` = as.numeric(as.character(`X2.50.`)))  %>%
      dplyr::mutate(`X97.50.` = as.numeric(as.character(`X97.50.`)))  %>%
      dplyr::mutate(mean = as.numeric(as.character(mean)))  %>%
      dplyr::mutate(sd = as.numeric(as.character(sd)))  %>%
      dplyr::mutate(`X75.` = as.numeric(as.character(`X75.`)))  %>%
      dplyr::mutate(`X25.` = as.numeric(as.character(`X25.`)))  %>%
      dplyr::mutate(Pr = as.numeric(as.character(Pr))) %>%
      
      # Create a boolean for strong or no support for difference from zero
      dplyr::mutate(Pr_yn = (Pr > 0.95)) %>%
      
      # Update here includes Pr of strong, moderate, or no support for difference from zero
      dplyr::mutate(Pr_yn = ifelse(Pr > 0.95, "strong", ifelse(Pr > 0.9, "moderate", "none"))) %>%
      
      # Attach to data containing parse-able names
      dplyr::full_join(orders, by = "measure")
    
    # Remove "recovery" epithet from several measures (not desired in plotting, instead use a facet)
    df$measure <- gsub("_recovery", "", df$measure)
    
    return(df)
  }


make_effect_plot <-
  function(genotype_comparison = F,
           recovery = F) {
    # This function makes the effect plots for treatment (general) OR each pairwise comparison
    # for the different genotypes
    
    effect_names <- c(
      `trt_effect` = "Treatment~effect",
      `G11[R]-G2[R]` = "G11[R]-G2[R]",
      `G11[R]-G5` = "G11[R]-G5",
      `G2[R]-G5` = "G2[R]-G5",
      #`int_effect` = "Interaction~effect", # Excluded from the labeller for now
      `Growth` = "Growth",
      `Instantaneous` = "Instantaneous",
      `Cumulative` = "Cumulative",
      `Recovery` = "Recovery"
    )
    
    # Gather each subset of data depending on the param / effect desired for plotting
    # Divide up so that treatment effects can be plotted separately from geno effects
    if (!(genotype_comparison)) {
      # Get treatment data only
      df <- clean_posterior_data_for_plotting("trt_effect")
      weak_effect_color <-
        no_effect_color # No weak effects, so replace that color with none
    } else {
      # Get geno data only
      geno_data1 <-
        clean_posterior_data_for_plotting("G11[R]-G2[R]")
      geno_data2 <- clean_posterior_data_for_plotting("G11[R]-G5")
      geno_data3 <- clean_posterior_data_for_plotting("G2[R]-G5")
      df <- rbind(geno_data1, geno_data2, geno_data3)
    }
    
    # Make a reordered factor to order facets (top)
    df$param_f = factor(df$param,
                        levels = c('trt_effect', 'G11[R]-G2[R]', 'G11[R]-G5', 'G2[R]-G5'))
    
    if (!(recovery)) {
      # DROP recovery
      df <-
        df %>% dplyr::filter(facet_left != "Recovery")
    } else {
      # OR, only look at recovery
      df <-
        df %>% dplyr::filter(facet_left == "Recovery")
      
      # Drop total since I didn't really use this in the experimental phase
      df <- df[(df$measure != "B_total"),]
    }
    
    # Make a reordered factor to order facets (left)
    df$facet_left_f = factor(df$facet_left,
                             levels = c('Growth', 'Instantaneous', 'Cumulative', 'Recovery'))
    
    # Order full name so it falls alphabetically downward on the y axis
    df_sort <-
      within(df, short_parse <-
               ordered(short_parse, levels = rev(sort(
                 unique(short_parse)
               ))))
    
    # Order level of evidence in the legend
    df_sort <-
      within(df_sort, Pr_yn <-
               ordered(Pr_yn, levels = c("strong", "moderate", "none")))
    
    # Make the plot
    gg <-
      ggplot(data = df_sort, aes(y = short_parse)) +
      
      # Add a zero line and x axis label
      geom_vline(xintercept = 0, color = zero_line_col) +
      xlab("Standard deviations") +
      
      # Add boxplots
      geom_boxplot(
        aes(
          xmin = `X2.50.` / sd,
          xmax = `X97.50.` / sd,
          xmiddle = mean / sd,
          xupper = `X75.` / sd,
          xlower = `X25.` / sd,
          fill = Pr_yn
        ),
        stat = "identity"
      ) +
      
      # Add general theme
      theme_cowplot() +
      
      # Parse y labels
      scale_y_discrete(breaks = levels(df_sort$short_parse),
                       labels = parse(text = levels(df_sort$short_parse))) +
      
      # Custom boxplot fill
      scale_fill_manual(
        name = "Effect",
        values = c(strong_effect_color, weak_effect_color, no_effect_color)
      ) +
      
      #Customize legend
      theme(
        legend.position = "right",
        axis.title.y = element_blank(),
        strip.placement = "outside"
      )
    
    # Facet according to effect type (param_f) and phenotypic measure grouping (facet_left_f)
    if (!(recovery)) {
      # Keep left facet if it's not recovery, since it's useful
      gg <-
        gg +
        
        # Add facets
        facet_grid(
          facet_left_f ~ param_f,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = as_labeller(effect_names, label_parsed)
        )
    } else {
      # The left facet just says recovery for these, it's not really helpful
      gg <-
        gg +
        
        # Add facets
        facet_grid(
          . ~ param_f,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = as_labeller(effect_names, label_parsed)
        )
    }
    
    return(gg)
  }

