###########################################################################################
## CODE FOR PLOTTING
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)


clean_posterior_data_for_plotting <-
  function(param_){
    # This function ..
    
    orders <- 
      as.tbl(read.csv("data/measure_order.csv", header = T)) %>% 
      dplyr::mutate(measure = as.character(measure))
    
    df <-
      read.csv("output/posterior_output.csv", header = T) %>%
      dplyr::filter(param == param_) %>%
      #dplyr::filter(measure %in% c("Bv", "Rt", "Rl")) %>%
      dplyr::mutate(`X2.5.` = as.numeric(as.character(`X2.50.`)))  %>%
      dplyr::mutate(`X97.5.` = as.numeric(as.character(`X97.50.`)))  %>%
      dplyr::mutate(mean = as.numeric(as.character(mean)))  %>%
      dplyr::mutate(sd = as.numeric(as.character(sd)))  %>%
      dplyr::mutate(`X75.` = as.numeric(as.character(`X75.`)))  %>%
      dplyr::mutate(`X25.` = as.numeric(as.character(`X25.`)))  %>%
      dplyr::mutate(Pr = as.numeric(as.character(Pr))) %>% 
      dplyr::mutate(Pr_yn = (Pr > 0.95)) %>%
      dplyr::full_join(orders, by = "measure")
    
    return(df)
  }


make_effect_plot <- 
  function(param_,
           outfile){
    
    trt_data <- clean_posterior_data_for_plotting("trt_effect")
    geno_data <- clean_posterior_data_for_plotting("geno_effect")
    int_data <- clean_posterior_data_for_plotting("int_effect")
    df <- rbind(trt_data, geno_data, int_data)
    
    # Make a reordered factor to order facets
    df$param_f = factor(df$param, levels=c('trt_effect','geno_effect','int_effect'))
    
    # Make a reordered factor to order facets
    df$facet_left_f = factor(df$facet_left, levels=c('Growth','Instantaneous','Cumulative','Recovery'))
    
    
    ggplot(data = df, aes(x = measure)) +
      geom_hline(yintercept = 0) +
      coord_flip() +
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
            axis.title.y = element_blank()) +
      facet_grid(facet_left_f~param_f, scales = "free_y", space = "free_y")
    
      ggsave(file = outfile, height = 8, width = 4)
    
  }
