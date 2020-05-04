###########################################################################################
# MAKE THE THEORETICAL FIGURE FEATURED IN THE MANUSCRIPT
###########################################################################################
library(ggplot2)
library(viridis)
library(cowplot)


make_subplot <-
  function(df) {
    # This function makes the subplots
    
    gg <-
      ggplot(df, aes(
        x = trt,
        y = value,
        color = as.factor(geno)
      )) +
      ylab("Phenotype") +
      ylim(c(0, 12)) +
      geom_line() +
      geom_point(size = 2) +
      theme_classic() +
      theme(
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank()
      ) +
      scale_color_manual(values = c("#3690c0", "#ef6548")) +
      scale_x_continuous(breaks = c(10, 20), labels = c("dry", "wet"))
    
    return(gg)
  }


make_theor_fig <-
  function(outfile) {
    # This function makes the theoretical figure
    #
    # outfile: sub path containing the name of the output theoretical figure
    
    # Gather data
    df <- read.csv("data/theoretical_fig.csv", header = T)
    d1 <- df[grep("1", df$plot), ]
    d2 <- df[grep("2", df$plot), ]
    d3 <- df[grep("3", df$plot), ]
    d4 <- df[grep("4", df$plot), ]
    
    # Make each subfigure
    # Get rid of x axis labels
    gg1 <- make_subplot(d1) + theme(axis.text.x = element_text(color = "transparent"))
    # Get rid of x axis labels and y label
    gg2 <- make_subplot(d2) + ylab("") + theme(axis.text.x = element_text(color = "transparent"))
    gg3 <- make_subplot(d3)
    # Get rid of y label
    gg4 <- make_subplot(d4) + ylab("")
    
    # Arrange and add labels
    grid <-
      plot_grid(
        gg1,
        gg2,
        gg3,
        gg4,
        align = 'v',
        labels = c("(a)", "(b)", "(c)", "(d)"),
        label_size = 10,
        hjust = -2,
        nrow = 2,
        rel_widths = c(1, 1)
      )
    
    # Make a legend
    grobs <-
      ggplotGrob(gg1 + theme(
        legend.title = element_blank(),
        legend.position = "top", 
        legend.box.just = "left"))$grobs
    legend <-
      grobs[[which(sapply(grobs, function(x)
        x$name) == "guide-box")]]
    
    # Plot and save given outfile arg
    plot_grid(grid,
              legend,
              ncol = 1,
              rel_heights = c(1, .15))
    ggsave(file = outfile,
           height = 3,
           width = 3)
    
    dev.off()
  }
