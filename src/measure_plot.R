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
breakpoint_pal <- c(colfunc(6), '#fdd49e')
zero_line_col <- "#bdbdbd"


clean_posterior_data_for_plotting <-
  function(param_) {
    # This function ..
    
    orders <-
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>%
      dplyr::mutate(measure = as.character(measure))
    
    df <-
      read.csv("output/posterior_output.csv", header = T) %>%
      dplyr::filter(param == param_) %>%
      dplyr::mutate(`X2.5.` = as.numeric(as.character(`X2.5.`)))  %>%
      dplyr::mutate(`X97.5.` = as.numeric(as.character(`X97.5.`)))  %>%
      dplyr::mutate(mean = as.numeric(as.character(mean)))  %>%
      dplyr::mutate(sd = as.numeric(as.character(sd)))  %>%
      dplyr::mutate(`X75.` = as.numeric(as.character(`X75.`)))  %>%
      dplyr::mutate(`X25.` = as.numeric(as.character(`X25.`)))  %>%
      dplyr::mutate(Pr = as.numeric(as.character(Pr))) %>%
      dplyr::mutate(Pr_yn = (Pr > 0.95)) %>%
      dplyr::mutate(Pr_yn = ifelse(Pr > 0.95, "strong", ifelse(Pr > 0.9, "moderate", "none"))) %>%
      dplyr::full_join(orders, by = "measure")
    
    df$measure <- gsub("_recovery", "", df$measure)
    
    return(df)
  }


make_effect_plot <-
  function() {
    # This function..
    
    effect_names <- c(
      `trt_effect` = "Treatment effect",
      `geno_effect` = "Genotype effect",
      `int_effect` = "Interaction effect",
      `Growth` = "Growth",
      `Instantaneous` = "Instantaneous",
      `Cumulative` = "Cumulative",
      `Recovery` = "Recovery"
    )
    
    trt_data <- clean_posterior_data_for_plotting("trt_effect")
    geno_data <- clean_posterior_data_for_plotting("geno_effect")
    int_data <- clean_posterior_data_for_plotting("int_effect")
    df <- rbind(trt_data, geno_data, int_data)
    
    # Make a reordered factor to order facets
    df$param_f = factor(df$param,
                        levels = c('trt_effect', 'geno_effect', 'int_effect'))
    
    # Make a reordered factor to order facets
    df$facet_left_f = factor(df$facet_left,
                             levels = c('Growth', 'Instantaneous', 'Cumulative', 'Recovery'))
    
    # Order full name so it falls alphabetically downward on the y axis
    df_sort <-
      within(df, short_parse <-
               ordered(short_parse, levels = rev(sort(
                 unique(short_parse)
               ))))
    
    df_sort <-
      within(df_sort, Pr_yn <-
               ordered(Pr_yn, levels = c("strong","moderate","none")))
    
    
    gg <- ggplot(data = df_sort, aes(y = short_parse)) +
      geom_vline(xintercept = 0, color = zero_line_col) +
      xlab("Standard deviations") +
      geom_boxplot(
        aes(
          xmin = `X2.5.` / sd,
          xmax = `X97.5.` / sd,
          xmiddle = mean / sd,
          xupper = `X75.` / sd,
          xlower = `X25.` / sd,
          fill = Pr_yn
        ),
        stat = "identity"
      ) +
      theme_cowplot() +
      facet_grid(
        facet_left_f ~ param_f,
        scales = "free",
        space = "free",
        switch = "y",
        labeller = as_labeller(effect_names)
      ) +
      scale_y_discrete(breaks = levels(df_sort$short_parse),
                       labels = parse(text = levels(df_sort$short_parse))) +
      scale_fill_manual(name = "Effect",
                        values = c(strong_effect_color, weak_effect_color, no_effect_color)) +
      theme(
        legend.position = "right",
        axis.title.y = element_blank(),
        strip.placement = "outside"
      )
    
    return(gg)
  }


clean_breakpoint_data_for_plotting <-
  function() {
    orders <-
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>%
      dplyr::mutate(measure = as.character(measure))
    
    orders$measure <- gsub("_recovery", "", orders$measure)
    
    df <-
      read.csv("output/breakpoint_analysis.csv", header = T)  %>%
      dplyr::full_join(orders, by = c("measure", "subset" = "facet_left"))
    
    df$Est. <- round(df$Est., 0)
    
    # Fill zeros with NAs since we don't care to plot them in the heatmap
    df$Est.[(df$Est. == 0)] <- "Linear"
    
    return(df)
  }


make_breakpoint_plot <-
  function() {
    df <- clean_breakpoint_data_for_plotting()
    
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
    
    gg <-
      ggplot(data = df_sort, aes(
        y = short_parse,
        x = as.factor(facet_geno),
        fill = `Est.`
      )) +
      theme_cowplot() +
      geom_tile() +
      facet_grid(
        facet_left_f ~ facet_geno,
        scales = "free",
        space = "free",
        labeller = label_parsed
      ) +
      xlab("Breakpoint estimate") +
      scale_y_discrete(breaks = levels(df_sort$short_parse),
                       labels = parse(text = levels(df_sort$short_parse))) +
      scale_fill_manual(values = c(breakpoint_pal)) +
      geom_text(data = df_sort[(df_sort$Est. != "Linear"), ], aes(label = paste(Est., "%", sep = ""))) +
      geom_text(data = df_sort[(df_sort$Est. == "Linear"), ], aes(label = Est.)) +
      theme(
        legend.position = "none",
        axis.title.y = element_text(color = "transparent"),
        #axis.title.x = element_text(color = "transparent"),
        strip.text.y = element_blank(),
        strip.placement = "outside",
        axis.text.x = element_text(color = "transparent"),
        axis.ticks.x = element_line(color = "transparent")
      )
    
    return(gg)
  }


make_fig_effects <-
  function(outfile = NA) {
    # This function gathers the two subplots, makes one big figure, and saves it if outfile is specified
    gd <- 
      plot_grid(
      make_effect_plot(),
      make_breakpoint_plot(),
      nrow = 1,
      rel_widths = c(0.7, 0.3),
      labels = c("(a)", "(b)")
    )
    
    if (is.na(outfile)){
      return(gd)
    } else {
      ggsave(plot = gd, 
             filename = outfile,
             height = 9,
             width = 18)
    }
  }
