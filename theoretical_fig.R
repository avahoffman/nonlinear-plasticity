wd <- '/Users/avahoffman/Dropbox/Research/Andropogon_geno_drought_study_2014/Revision_1'
setwd(wd)
df <- read.csv("theoretical_fig.csv",header=T)
d1 <- df[grep("1", df$plot), ]
d2 <- df[grep("2", df$plot), ]
d3 <- df[grep("3", df$plot), ]
d4 <- df[grep("4", df$plot), ]


gg1 <- 
  ggplot( d1, aes(x=trt,y=value,color=as.factor(geno))) +
  ylab("Trait value    ") +
  ylim(c(0,12)) +
  geom_line() +
  geom_point(size=2) +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_color_viridis(discrete=T, begin=0.2, end=0.7, name="Genotype",option="inferno") +
  scale_x_continuous(breaks=c(10,20), labels=c("dry","wet"))


gg2 <- 
  ggplot( d2, aes(x=trt,y=value,color=as.factor(geno))) +
  ylab("Trait value    ") +
  ylim(c(0,12)) +
  geom_line() +
  geom_point(size=2) +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_color_viridis(discrete=T, begin=0.2, end=0.7, name="Genotype",option="inferno") +
  scale_x_continuous(breaks=c(10,20), labels=c("dry","wet"))

gg3 <- 
  ggplot( d3, aes(x=trt,y=value,color=as.factor(geno))) +
  ylab("Trait value    ") +
  ylim(c(0,12)) +
  geom_line() +
  geom_point(size=2) +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_color_viridis(discrete=T, begin=0.2, end=0.7, name="Genotype",option="inferno") +
  scale_x_continuous(breaks=c(10,20), labels=c("dry","wet"))

gg4 <- 
  ggplot( d4, aes(x=trt,y=value,color=as.factor(geno))) +
  ylab("Trait value    ") +
  ylim(c(0,12)) +
  geom_line() +
  geom_point(size=2) +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_color_viridis(discrete=T, begin=0.2, end=0.7, name="Genotype",option="inferno") +
  scale_x_continuous(breaks=c(10,20), labels=c("dry","wet"))

grid=plot_grid( gg1,gg2,gg3,gg4, align = 'v', labels = c("(a)", "(b)", "(c)", "(d)"), label_size=10, hjust = -2, nrow = 2, rel_widths=c(1,1) )
grobs=ggplotGrob(gg1 + theme(legend.position="top",legend.box.just = "left"))$grobs
legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
pdf(file="theoretical_fig.pdf",height=3,width=3)
plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .15)) ; dev.off() ; dev.off() ; dev.off()

