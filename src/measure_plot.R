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
    if (!(genotype_comparison)) {
      df <- clean_posterior_data_for_plotting("trt_effect")
      weak_effect_color <-
        no_effect_color # No weak effects, so replace that color with none
    } else {
      geno_data1 <- clean_posterior_data_for_plotting("G11[R]-G2[R]")
      geno_data2 <- clean_posterior_data_for_plotting("G11[R]-G5")
      geno_data3 <- clean_posterior_data_for_plotting("G2[R]-G5")
      df <- rbind(geno_data1, geno_data2, geno_data3)
    }
    
    # Make a reordered factor to order facets (top)
    df$param_f = factor(df$param,
                        levels = c('trt_effect', 'G11[R]-G2[R]', 'G11[R]-G5', 'G2[R]-G5'))
    
    # DROP recovery for now
    if (!(recovery)) {
      df <-
        df %>% dplyr::filter(facet_left != "Recovery")
    } else {
      df <-
        df %>% dplyr::filter(facet_left == "Recovery")
      # Drop total since I didn't really use this in the experimental phase
      df <- df[(df$measure != "B_total"), ]
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
      gg <-
        gg +
        facet_grid(
          facet_left_f ~ param_f,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = as_labeller(effect_names, label_parsed)
        )
    } else {
      gg <-
        gg +
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


clean_breakpoint_data_for_plotting <-
  function() {
    # This function gathers breakpoint results (output) for plotting
    
    # Data containing full / parse-able names for phenotypic measures
    orders <-
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>%
      dplyr::mutate(measure = as.character(measure))
    
    # Remove "recovery" epithet from several measures (not desired in plotting, instead use a facet)
    orders$measure <- gsub("_recovery", "", orders$measure)
    
    # Attach to data containing parse-able names
    df <-
      read.csv("output/breakpoint_analysis.csv", header = T)  %>%
      dplyr::full_join(orders, by = c("measure", "subset" = "facet_left"))
    
    # Round estimates to the nearest percent
    df$Est. <- round(df$Est., 0)
    
    # Fill zeros with NAs since we don't care to plot them in the heatmap
    df$Est.[(df$Est. == 0)] <- "Linear"
    
    return(df)
  }


make_breakpoint_plot <-
  function(recovery = F) {
    # This function makes the breakpoint plots for the different genotypes
    
    # Read in data
    df <- clean_breakpoint_data_for_plotting()
    
    # DROP recovery for now
    if (!(recovery)) {
      df <-
        df %>% dplyr::filter(subset != "Recovery")
    } else {
      df <-
        df %>% dplyr::filter(subset == "Recovery")
      breakpoint_pal <- c(colfunc(5), '#fdd49e')
      # Drop total since I didn't really use this in the experimental phase
      df <- df[(df$measure != "B_total"), ]
    }
    
    # Make a reordered factor to order facets
    df$facet_left_f = factor(df$subset,
                             levels = c('Growth', 'Instantaneous', 'Cumulative', 'Recovery'))
    
    # Make a reordered factor to order facets
    df$facet_geno = factor(df$geno,
                           levels = c('2', '11', '5'))
    df$facet_geno <- gsub("2", "G2[R]", df$facet_geno)
    df$facet_geno <- gsub("11", "G11[R]", df$facet_geno)
    df$facet_geno <- gsub("5", "G5", df$facet_geno)
    
    # Order full name so it falls alphabetically downward on the y axis
    df_sort <-
      within(df, short_parse <-
               ordered(short_parse, levels = rev(sort(
                 unique(short_parse)
               ))))
    
    # Make the plot
    gg <-
      # Fill the heatmap tile with the magnitude of the estimate
      ggplot(data = df_sort, aes(
        y = short_parse,
        x = as.factor(facet_geno),
        fill = `Est.`
      )) +
      
      # General theme
      theme_cowplot() +
      
      # "Heatmap" style
      geom_tile() +
      
      # X Axis label
      xlab("Breakpoint estimate") +
      
      # Parse y labels
      scale_y_discrete(breaks = levels(df_sort$short_parse),
                       labels = parse(text = levels(df_sort$short_parse))) +
      
      # Custom heatmap colors (the Linear one is done manually according to colors above..)
      scale_fill_manual(values = c(breakpoint_pal)) +
      
      # Where it's not linear, add the estimate and percentage symbol
      geom_text(data = df_sort[(df_sort$Est. != "Linear"), ],
                aes(label = paste(Est., "%", sep = ""))) +
      
      # Where it's linear, just say "Linear"
      geom_text(data = df_sort[(df_sort$Est. == "Linear"), ], aes(label = Est.)) +
      
      # Remove labels, gridlines, etc
      theme(
        legend.position = "none",
        axis.title.y = element_text(color = "transparent"),
        strip.placement = "outside",
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent")
      )
    
    # Facet according to genotype and phenotypic measure grouping (facet_left_f)
    if (!(recovery)) {
      gg <-
        gg +
        facet_grid(
          facet_left_f ~ facet_geno,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = label_parsed
        )
    } else {
      gg <-
        gg +
        facet_grid(
          . ~ facet_geno,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = label_parsed
        )
    }
    
    return(gg)
  }


make_fig_effects <-
  function(outfile = NA) {
    # This function gathers the two subplots, makes one big figure,
    # and saves it if outfile is specified
    
    # Plot two subplots in appropriate widths
    gd <-
      plot_grid(
        make_effect_plot(genotype_comparison = T),
        make_breakpoint_plot(),
        nrow = 1,
        rel_widths = c(0.7, 0.3),
        labels = c("(a)", "(b)")
      )
    
    # Write file or return gridded plot
    if (is.na(outfile)) {
      return(gd)
    } else {
      ggsave(
        plot = gd,
        filename = outfile,
        height = 9,
        width = 18
      )
    }
  }
