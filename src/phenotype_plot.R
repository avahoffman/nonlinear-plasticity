###########################################################################################
## CODE FOR PLOTTING REACTION NORMS
###########################################################################################
library(ggplot2)
library(cowplot)
library(dplyr)


measure_ <- "DMCv"

df <-
  read.csv("output/posterior_output.csv", header = T) %>%
  dplyr::filter(measure == measure_ & param == "posterior value")

df$trt <- as.numeric(as.character(gsub("Sat'd", 30, df$trt)))
facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
df <- cbind(df, facet)
  
  

ggplot(data = df, aes(x=trt,y=mean, color = geno)) +
  geom_errorbar(aes(ymin=`X2.5.`, ymax=`X97.5.`), width=0.3) +
  xlab("% VWC") +
  theme_cowplot() +
  geom_line() +
  geom_point( size = 2)  +
  scale_x_continuous(breaks=c(10,15,20,25,30), labels=c("10","15","20","25","Sat'd")) +
  facet_grid(.~facet, scales="free",space="free") +
  theme(strip.text = element_blank()) 
  


# 
# ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
#   geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
#   xlab("% VWC") +
#   ylab(response.label) +
#   geom_line() +
#   geom_point(size=1) +
#   theme_classic() +
#   theme(legend.position="bottom") +
#   scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
#   scale_x_continuous(breaks=c(10,15,20,25,3), labels=c("10","15","20","25","Sat'd")) +
#   facet_grid(.~facet, scales="free",space="free") +
#   theme(strip.text = element_blank()) +
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = -7, b = 0, l = 0)))
