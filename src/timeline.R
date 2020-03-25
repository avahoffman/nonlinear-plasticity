###########################################################################################
## MAKE AN EXPERIMENTAL TIMELINE
###########################################################################################
library(ggplot2)


make_timeline <-
  function() {
    # This function plots the timeline for the experiments found in the study
    
    # Read in some structural data of dates, phases, etc
    d <- read.csv("data/timeline.csv", header = T)
    
    # Convert date to datetime and not factor
    d$date <- as.Date(d$date)
    
    # Make timeline plot
    ggplot() +
      # Background color for different phases
      geom_col(data = d,
               aes(x = date, y = 1, fill = phase),
               color = "transparent") +
      theme_void() + # Get rid of most axes, etc
      
      # Add vertical lines for measurement points
      geom_segment(
        data = d[(d$subset != "None"), ],
        aes(
          y = bar_height,
          yend = 0.5,
          xend = date,
          x = date
        ),
        color = 'black',
        size = 0.2
      )  +
      
      # Points for different plant subsets
      geom_point(data = d[(d$subset != "None"), ],
                 aes(x = date, y = 0.5, color = subset),
                 size = 3) +
      
      # Measurement names
      geom_text(
        data = d[(d$subset != "None"), ],
        aes(x = date, y = label_pos, label = action) ,
        size = 2.5,
        hjust = 0.1
      )  +
      
      # Horizontal bar for height measurements
      geom_segment(aes(
        y = 1.5,
        yend = 1.5,
        xend = as.Date("2014-11-24"),
        x = as.Date("2014-08-21")
      ), size = 0.2) +
      
      # Vertical bar for height measurements
      geom_segment(aes(
        y = 1.5,
        yend = 2,
        xend = as.Date("2014-09-20"),
        x = as.Date("2014-09-20")
      ), size = 0.2) +
      
      # Horizontal bar for phys measurements
      geom_segment(aes(
        y = -0.5,
        yend = -0.5,
        xend = as.Date("2014-09-07"),
        x = as.Date("2014-10-11")
      ), size = 0.2) +
      
      # Vertical bar for phys measurements
      geom_segment(aes(
        y = -0.5,
        yend = -2.1,
        xend = as.Date("2014-09-30"),
        x = as.Date("2014-09-30")
      ), size = 0.2) +
      
      # Vertical bar for biomass measurements
      geom_segment(aes(
        y = -0.5,
        yend = -0.8,
        xend = as.Date("2014-10-20"),
        x = as.Date("2014-10-20")
      ), size = 0.2) +
      
      # Horizontal bar for biomass measurements 
      geom_segment(aes(
        y = -0.5,
        yend = -0.5,
        xend = as.Date("2014-10-16"),
        x = as.Date("2014-11-25")
      ), size = 0.2) +
      
      # Custom colors for phases 
      scale_fill_manual(
        values = c("#a6bddb", "#ece7f2", "#fdbb84", "#fee8c8"),
        labels = c(
          "Acclimation (misters)",
          "Acclimation (pots)",
          "Drought treatment",
          "Recovery"
        ),
        name = "Experimental phase"
      ) +
      
      # Custom colors for subsets of plants
      scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                         name = "Plant subset") +
      
      # Annotate month (not at x axis line though)
      geom_text(data = d,
                aes(x = date, y = -0.2, label = month),
                size = 4)
    
    # Save figure
    ggsave(file = "figures/timeline.pdf",
           height = 3,
           width = 10)
  }
