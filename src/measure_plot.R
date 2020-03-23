###########################################################################################
## CODE FOR PLOTTING
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)


make_effect_plot <- 
  function(param_,
           outfile){
    orders <- 
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>% 
      dplyr::mutate(measure = as.character(measure))
    
    df <-
      read.csv("output/posterior_output.csv", header = T) %>%
      dplyr::filter(param == param_) %>%
      #dplyr::filter(measure %in% c("Bv", "Rt", "Rl")) %>%
      dplyr::mutate(`X2.5.` = as.numeric(as.character(`X2.5.`)))  %>%
      dplyr::mutate(`X97.5.` = as.numeric(as.character(`X97.5.`)))  %>%
      dplyr::mutate(mean = as.numeric(as.character(mean)))  %>%
      dplyr::mutate(sd = as.numeric(as.character(sd)))  %>%
      dplyr::mutate(`X75.` = as.numeric(as.character(`X75.`)))  %>%
      dplyr::mutate(`X25.` = as.numeric(as.character(`X25.`)))  %>%
      dplyr::mutate(Pr = as.numeric(as.character(Pr))) %>% 
      dplyr::mutate(Pr_yn = (Pr > 0.95)) %>%
      dplyr::full_join(orders, by = "measure")
    
    ggplot(data = df, aes(x = order_desc)) +
      geom_hline(yintercept = 0) +
      coord_flip() +
      scale_x_discrete(labels = c(rev(df$measure))) +
      ylab("Standard deviations") +
      geom_boxplot(aes(
        ymin = `X2.5.`/ sd,
        ymax = `X97.5.`/ sd,
        middle = mean / sd,
        upper = `X75.` / sd,
        lower = `X25.` / sd,
        fill = Pr_yn
      ),
      stat = "identity"
      ) +
      geom_hline(yintercept = 0) +
      theme_cowplot() +
      theme(legend.position = "none",
            axis.title.y = element_blank()) 
    
      ggsave(file = outfile, height = 8, width = 4)
    
  }
