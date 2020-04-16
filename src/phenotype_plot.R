###########################################################################################
## CODE FOR PLOTTING REACTION NORMS
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)


make_phenotype_plot <-
  function(measure_, return_plot = T) {
    df <-
      read.csv("output/posterior_output.csv", header = T) %>%
      # Gather parsed names for axis label
      full_join(read.csv("data/measure_order.csv"),
                header = T,
                by = "measure") %>%
      dplyr::filter(measure == measure_ &
                      param == "posterior value")
    
    
    df$trt <- as.numeric(as.character(gsub("Sat'd", 30, df$trt)))
    df$parse_name <- as.character(df$parse_name)
    facet <- c(1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2)
    df <- cbind(df, facet)
    
    G11 <- parse(text = "paste(G11[R])")
    G2 <- parse(text = "paste(G2[R])")
    G5 <- "G5"
    
    genotype_colors <-
      c("#d7301f", "#fc8d59", "#3690c0")
    
    
    gg <-
      ggplot(data = df, aes(x = trt, y = mean, color = geno)) +
      theme_cowplot() +
      geom_errorbar(aes(ymin = `X2.50.`, ymax = `X97.50.`), width = 2) +
      geom_line() +
      geom_point(aes(shape = geno, color = geno),
                 size = 4,
                 fill = "white") +
      xlab("% VWC") +
      ylab(parse(text = df$parse_name[1])) +
      scale_x_continuous(
        breaks = c(10, 15, 20, 25, 30),
        labels = c("10", "15", "20", "25", "Sat'd")
      ) +
      facet_grid(. ~ facet, scales = "free", space = "free") +
      scale_color_manual(name = "Genotype",
                         labels = c(G11, G2, G5),
                         values = genotype_colors) +
      scale_shape_manual(
        name = "Genotype",
        labels = c(G11, G2, G5),
        values = c(23, 24, 21)
      ) +
      theme(
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.direction = "vertical",
        legend.position = element_blank()
      ) +
      theme(legend.position = "right")
    
    if (!(return_plot)) {
      ggsave(
        gg,
        file = paste("figures/phenotypes/", measure_, ".pdf", sep = ""),
        height = 6,
        width = 7
      )
    } else {
      return(gg)
    }
  }


cycle_phenotype_plots <-
  function() {
    # This function cycles through and makes phenotype plots for all given measures
    
    meas_list <- c(
      "urgr_recovery",
      "uH_recovery",
      "max_rgr_recovery",
      "max_H_recovery",
      "B_total_recovery",
      "B_below_recovery",
      "B_above_recovery",
      "A_B_recovery",
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
      "A_B",
      "uWUEi",
      "ugs",
      "ufv",
      "uAnet",
      "max_WUEi",
      "max_gs",
      "max_fv",
      "max_Anet",
      "urgr",
      "uH",
      "max_rgr",
      "max_H"
    )
    
    for (meas in meas_list) {
      make_phenotype_plot(meas, return_plot = F)
    }
    
  }


leg <-
  g_legend(make_phenotype_plot("Bv") + theme(legend.margin = margin(l = 0.6, unit =
                                                                      'cm')))

growth_summary <-
  function() {
    meas_list <- c("max_H",
                   "uH",
                   "max_rgr",
                   "urgr")
    
    p <- list()
    i = 0
    for (meas in meas_list) {
      i = i + 1
      p[[i]] <-
        make_phenotype_plot(meas) + theme(legend.position = "none")
    }
    
    grid <-
      plot_grid(
        p[[1]],
        p[[2]],
        p[[3]],
        p[[4]],
        align = "vh",
        axis = "bl",
        nrow = 2,
        labels = c("(a)", "(b)", "(c)", "(d)")
      )
    
    ggsave(
      gridExtra::grid.arrange(
        grid,
        geno_legend(make_phenotype_plot("Bv") + theme(legend.margin = margin(
          l = 0.6, unit =
            'cm'
        ))),
        nrow = 1,
        widths = c(7, 1)
      ),
      file = paste("figures/growth_summary.pdf", sep = ""),
      height = 9,
      width = 11
    )
    
  }


instantaneous_summary <-
  function() {
    meas_list <- c("uWUEi",
                   "ugs",
                   "ufv",
                   "uAnet",
                   "max_WUEi",
                   "max_gs",
                   "max_fv",
                   "max_Anet")
    
    p <- list()
    i = 0
    for (meas in meas_list) {
      i = i + 1
      p[[i]] <-
        make_phenotype_plot(meas) + theme(legend.position = "none")
    }
    
    grid <-
      plot_grid(
        p[[1]],
        p[[2]],
        p[[3]],
        p[[4]],
        p[[5]],
        p[[6]],
        p[[7]],
        p[[8]],
        align = "vh",
        axis = "bl",
        nrow = 2
      )
    
    ggsave(
      gridExtra::grid.arrange(
        grid,
        geno_legend(make_phenotype_plot("Bv") + theme(legend.margin = margin(
          l = 0.6, unit =
            'cm'
        ))),
        nrow = 1,
        widths = c(10, 1)
      ),
      file = paste("figures/instantaneous_summary.pdf", sep = ""),
      height = 9,
      width = 18
    )
    
  }


cumulative_summary <-
  function() {
    meas_list <- c(
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
      "A_B"
    )
    
    p <- list()
    i = 0
    for (meas in meas_list) {
      i = i + 1
      p[[i]] <-
        make_phenotype_plot(meas) + theme(legend.position = "none")
    }
    
    grid <-
      plot_grid(
        p[[1]],
        p[[2]],
        p[[3]],
        p[[4]],
        p[[5]],
        p[[6]],
        p[[7]],
        p[[8]],
        p[[9]],
        p[[10]],
        p[[11]],
        p[[12]],
        p[[13]],
        p[[14]],
        p[[15]],
        align = "vh",
        axis = "bl",
        nrow = 5
      )
    
    ggsave(
      gridExtra::grid.arrange(
        grid,
        geno_legend(make_phenotype_plot("Bv") + theme(legend.margin = margin(
          l = 0.6, unit =
            'cm'
        ))),
        nrow = 1,
        widths = c(10, 1)
      ),
      file = paste("figures/cumulative_summary.pdf", sep = ""),
      height = 16,
      width = 12
    )
    
  }


recovery_summary <-
  function() {
    meas_list <- c(
      "urgr_recovery",
      "uH_recovery",
      "max_rgr_recovery",
      "max_H_recovery",
      "B_total_recovery",
      "B_below_recovery",
      "B_above_recovery",
      "A_B_recovery"
    )
    
    p <- list()
    i = 0
    for (meas in meas_list) {
      i = i + 1
      p[[i]] <-
        make_phenotype_plot(meas) + theme(legend.position = "none")
    }
    
    grid <-
      plot_grid(
        p[[1]],
        p[[2]],
        p[[3]],
        p[[4]],
        p[[5]],
        p[[6]],
        p[[7]],
        p[[8]],
        align = "vh",
        axis = "bl",
        nrow = 2
      )
    
    ggsave(
      gridExtra::grid.arrange(
        grid,
        geno_legend(make_phenotype_plot("Bv") + theme(legend.margin = margin(
          l = 0.6, unit =
            'cm'
        ))),
        nrow = 1,
        widths = c(10, 1)
      ),
      file = paste("figures/recovery_summary.pdf", sep = ""),
      height = 9,
      width = 18
    )
    
  }
