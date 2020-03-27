###########################################################################################
## CODE FOR PLOTTING REACTION NORMS
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)


make_phenotype_plot <-
  function(measure_){
    df <-
      read.csv("output/posterior_output.csv", header = T) %>%
      # Gather parsed names for axis label
      full_join(read.csv("data/measure_order.csv"), header = T, by = "measure") %>%
      dplyr::filter(measure == measure_ & param == "posterior value")
    
    
    df$trt <- as.numeric(as.character(gsub("Sat'd", 30, df$trt)))
    df$parse_name <- as.character(df$parse_name)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
    df <- cbind(df, facet)
    
    G11 <- parse(text = "paste(G11[R])")
    G2 <- parse(text = "paste(G2[R])")
    G5 <- "G5"
    
    genotype_colors <-
      c("#d7301f", "#fc8d59", "#3690c0")
    
    
    gg <- 
      ggplot(data = df, aes(x = trt, y = mean, color = geno)) +
      theme_cowplot() +
      geom_errorbar(aes(ymin = `X2.5.`, ymax = `X97.5.`), width = 2) +
      geom_line() +
      geom_point(aes(shape = geno, color = geno), size = 4, fill = "white") +
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
      theme(strip.text = element_blank(),
            legend.title = element_blank(),
            legend.text.align = 0,
            legend.direction = "vertical",
            legend.position = element_blank()) +
      theme(legend.position = "right")
    
    ggsave(gg, file = paste("figures/phenotypes/", measure_, ".pdf", sep = ""),
           height = 6,
           width = 7)
  }


cycle_phenotype_plots <- 
  function(){
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
    
    for (meas in meas_list){
      make_phenotype_plot(meas)
    }
    
  }
