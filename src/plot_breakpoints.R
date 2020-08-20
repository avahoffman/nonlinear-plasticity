###########################################################################################
## CODE FOR PLOTTING BREAKPOINT ANALYSIS
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)

# Colors ----
strong_effect_color <- "#d0d1e6"
weak_effect_color <- "#fdd49e"
no_effect_color <- "#fc8d59"
colfunc <- colorRampPalette(c("#fff7fb", '#74a9cf'))
breakpoint_pal <- c(colfunc(14), '#fdd49e')
zero_line_col <- "#bdbdbd"


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
    
    df$interval <- df$upper - df$lower
    
    # Fill zeros with NAs since we don't care to plot them in the heatmap
    df$Est.[(df$`Est.` == 0 | 
               # df$lower < 5 | 
               # df$upper > 35 | 
               # df$`St.Err` > 4 |
               df$pvalue > 0.05 )] <- "Linear"
    
    return(df)
  }


make_breakpoint_plot <-
  function(recovery = F) {
    # This function makes the breakpoint plots for the different genotypes
    
    # Read in data
    df <- clean_breakpoint_data_for_plotting()
    breakpoint_pal <- c(colfunc(length(unique(df$Est.)) - 1), '#fdd49e')
    
    if (!(recovery)) {
      # DROP recovery for now
      df <-
        df %>% dplyr::filter(subset != "Recovery")
    } else {
      # OR only look at recovery
      df <-
        df %>% dplyr::filter(subset == "Recovery")
      # Drop total since I didn't really use this in the experimental phase
      df <- df[(df$measure != "B_total"),]
    }
    
    # Make a reordered factor to order facets
    df$facet_left_f = factor(df$subset,
                             levels = c('Growth', 'Instantaneous', 'Cumulative', 'Recovery'))
    
    # Make a reordered factor to order facets
    df$facet_geno = factor(df$geno,
                           levels = c('2', '11', '5'))
    df$facet_geno <- gsub("2", "G2[S]", df$facet_geno)
    df$facet_geno <- gsub("11", "G11[S]", df$facet_geno)
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
      geom_text(data = df_sort[(df_sort$Est. != "Linear"),],
                aes(label = paste(Est., "%", sep = ""))) +
      
      # Where it's linear, just say "Linear"
      geom_text(data = df_sort[(df_sort$Est. == "Linear"),], aes(label = Est.)) +
      
      # Remove labels, gridlines, etc
      theme(
        legend.position = "none",
        axis.title.y = element_text(color = "transparent"),
        strip.placement = "outside",
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent")
      ) + theme_sigmaplot(ticklen = -0.15)
    
    # Facet according to genotype and phenotypic measure grouping (facet_left_f)
    if (!(recovery)) {
      # Keep left facet if it's not recovery, since it's useful
      gg <-
        gg +
        
        # Add facets
        facet_grid(
          facet_left_f ~ facet_geno,
          scales = "free",
          space = "free",
          switch = "y",
          labeller = label_parsed
        )
    } else {
      # Facets not super useful on lefthand side since it's just recovery
      gg <-
        gg +
        
        # Add facets
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
