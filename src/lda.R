###########################################################################################
## CODE TO PERFORM LDA TO EXAMINE DIFFERENCES BETWEEN GENOTYPES
###########################################################################################


prepare_lda_data <-
  function(infile,
           responses,
           leave_out_indicator,
           subset,
           v_just,
           spoke_scale_factor) {
    # This function..
    
    # Read in data
    df <- read.csv(infile, header = T) %>%
      tidyr::drop_na()
    
    # Only include desired response variables in the argument
    df <- df[, c("geno", responses)]
    group_var <- c("geno")
    
    # LDA formula
    lda_form <- as.formula(paste("geno ~ ."))
    
    # Features are everything but the geno
    feats <-
      df %>%
      dplyr::select(-c(group_var))
    
    # Scale in prep for LDA
    feats_scaled <- scale(feats)
    
    # Run LDA and save proportion of variance explained
    ldadat <-
      cbind(df %>% dplyr::select(group_var),
            as.data.frame(feats_scaled))
    lda_results <- lda(formula = lda_form, data = ldadat)
    prop.lda = lda_results$svd ^ 2 / sum(lda_results$svd ^ 2)
    prop.lda
    plda <- predict(object = lda_results, newdata = ldadat)
    
    # Save the loadings
    load_dat <-
      data.frame(varnames = rownames(coef(lda_results)), coef(lda_results))
    load_dat$length <- with(load_dat, sqrt(LD1 ^ 2 + LD2 ^ 2))
    load_dat$angle <- atan2(load_dat$LD2, load_dat$LD1)
    load_dat$x_start <- load_dat$y_start <- 0
    load_dat$x_end <-
      cos(load_dat$angle) * log10(load_dat$length) * spoke_scale_factor # This sets the length of your lines.
    load_dat$y_end <-
      sin(load_dat$angle) * log10(load_dat$length) * spoke_scale_factor # This sets the length of your lines.
    
    # Data containing full / parse-able names for phenotypic measures
    orders <-
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>%
      dplyr::mutate(varnames = as.character(measure))
    load_dat <-
      load_dat %>% left_join(orders, by = "varnames")
    # Add leave out indicator and v_just for plotting
    load_dat <-
      cbind(load_dat, leave_out_indicator, v_just)
    print(load_dat)
    
    # save only top loadings in LD1 and LD2
    load_dat <-
      load_dat %>%
      mutate(h_just = replace(x_end, x_end < 0, 1)) %>%
      mutate(h_just = replace(h_just, x_end >= 0, 0)) %>%
      dplyr::filter(leave_out_indicator != 0)
    
    subset_f <- rep(subset, nrow(plda$x))
    plotdat <-
      cbind(data.frame(pop = df %>%
                         dplyr::select(geno), plda$x),
            subset_f)
    
    centroids <- aggregate(cbind(LD1, LD2) ~ geno, plotdat, mean)
    
    return(list(plotdat, prop.lda, load_dat, centroids, spoke_scale_factor))
  }


make_lda_plot <-
  function(d, vjust, xlim = NA, ylim = NA) {
    # This function..
    
    genotype_colors <-
      c("#d7301f", "#fc8d59", "#3690c0")
    
    G11 <- parse(text = "paste(G11[R])")
    G2 <- parse(text = "paste(G2[R])")
    G5 <- "G5"
    
    gg <-
      ggplot(d[[1]], aes(LD1, LD2)) +
      theme_cowplot() +
      geom_point(
        aes(
          fill = as.factor(geno),
          color = as.factor(geno),
          shape = as.factor(geno)
        ),
        size = 5,
        data = d[[4]]
      ) +
      geom_point(aes(color = as.factor(geno),
                     shape = as.factor(geno)), size = 2.5) +
      scale_color_manual(name = "Genotype",
                         labels = c(G11, G2, G5),
                         values = genotype_colors) +
      scale_fill_manual(name = "Genotype",
                        labels = c(G11, G2, G5),
                        values = genotype_colors) +
      scale_shape_manual(
        name = "Genotype",
        labels = c(G11, G2, G5),
        values = c(23, 24, 21)
      ) +
      labs(
        x = paste("LD1 (", percent(d[[2]][1]), ")", sep = ""),
        y = paste("LD2 (", percent(d[[2]][2]), ")", sep = "")
      ) +
      labs(colour = "Site") +
      geom_spoke(
        aes(
          x_start,
          y_start,
          angle = angle,
          radius = log10(length) * d[[5]]
        ),
        d[[3]],
        color = "gray",
        size = 0.5,
        show.legend = FALSE
      ) +
      geom_label(
        aes(
          y = y_end,
          x = x_end,
          label = short_parse,
          hjust = h_just,
          vjust = v_just
        ),
        #can change alpha=length within aes
        d[[3]],
        alpha = 0.6,
        size = 3,
        colour = "black",
        show.legend = FALSE,
        parse = TRUE
      ) +
      
      # Facet according to phenotypic measure grouping (subset_f)
      facet_grid(
        . ~ subset_f,
        switch = "y",
        # Parse label name
        labeller = as_labeller(effect_names, label_parsed)
      ) +
      
      guides(text = FALSE,
             spoke = FALSE,
             length = FALSE) +
      theme(legend.title = element_blank(),
            legend.position = "none",
            legend.text.align = 0)

    if(!(is.na(xlim))){
      gg <- gg + xlim(-xlim, xlim)
    }
    if(!(is.na(ylim))){
      gg <- gg + ylim(ylim)
      }
    
    
    return(gg)
    
  }


gather_lda_plots <-
  function(outfile = NA) {
    # This function gathers and writes (if indicated) LDA ordination plots
    
    # Growth subset
    d1 <-
      prepare_lda_data(
        infile = "data/all_plants_clean.csv",
        responses =  c("urgr",
                       "uH",
                       "max_rgr",
                       "max_H",
                       "trt"),
        subset = "Growth",
        leave_out_indicator = c(1, 1, 1, 1, 0),
        v_just = c(1,0,1,0.5,0.5),
        spoke_scale_factor = 4
      )
    
    # Cumulative subset
    d2 <-
      prepare_lda_data(
        infile = "data/biomass_plants_clean.csv",
        responses =  c(
          "ssa",
          "srl",
          "SLA",
          "Rv",
          "Rt",
          "Rsa_la",
          "Rsa",
          "Rl",
          "Rd",
          #"LA",
          "DMCv",
          "DMCr",
          "Bv",
          "Br",
          "A_B",
          "trt"
        ),
        subset = "Cumulative",
        leave_out_indicator = c(1, 1, 1, 1, 1,
                                0, 1, 1, 1, 1,
                                1, 1, 1, 1, 0),
        v_just = c(0.5,0.5,0.5,0.5,0.5,
                   0.5,0.5,0.5,0.5,0.5,
                   0.5,0.5,0.5,0.5,0.5),
        spoke_scale_factor = 4
      )
    
    # Instantaneous subset
    d3 <-
      prepare_lda_data(
        infile = "data/phys_plants_clean.csv",
        responses =  c("uWUEi",
                       "ugs",
                       "ufv",
                       "uAnet",
                       "max_WUEi",
                       #"max_gs",
                       "max_fv",
                       #"max_Anet",
                       "trt"),
        subset = "Instantaneous",
        leave_out_indicator = c(1, 1, 1, 1,
                                1, 1, 0),
        v_just = c(0.5,1,1,0,
                   0.5,0.5,0.5),
        spoke_scale_factor = 4
      )
    
    leg <- g_legend(make_lda_plot(d1) + theme(legend.position = "right",
                                              # Spread out legend keys a little bit more
                                              legend.key.size = unit(0.7, "cm")))
    
    # Plot two subplots in appropriate widths
    gd <-
      plot_grid(
        make_lda_plot(d1),
        make_lda_plot(d3),
        make_lda_plot(d2, xlim = 6),
        leg,
        nrow = 1,
        rel_widths = c(1, 1, 1, 0.2),
        labels = c("(a)", "(b)", "(c)", "")
      )
    
    # Write file or return gridded plot
    if (is.na(outfile)) {
      return(gd)
    } else {
      ggsave(
        plot = gd,
        filename = outfile,
        height = 5,
        width = 15
      )
    }
  }

gather_lda_plots()
