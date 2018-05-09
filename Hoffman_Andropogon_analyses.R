###########################################################################################
##
## R source code to accompany Hoffman & Smith (2018), last updated 18 April 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## If you found this code useful, please use the citation below:
## 
##
###########################################################################################

## SET YOUR WORKING DIRECTORY TO WHEREVER YOU HAVE DOWNLOADED ACCOMPANYING FILES
wd <- '/Users/avahoffman/Dropbox/Research/Andropogon_geno_drought_study_2014/Revision_1'
setwd(wd)
library(reshape2)
library(dplyr)
library(ggplot2)
library(viridis)
library(cowplot)
sessionInfo()

# model1 <- lm(srl~as.factor(geno)*as.factor(trt),data=df)
# summary(model1)

###########################################################################################
## FIRST ORDER OF BUSINESS IS DETERMINING WHETHER TRAITS ARE CORRELATED FOR FURTHER ANALYSIS
###########################################################################################

df <- read.csv("biomass_plants.csv",header=T)
    dir.create(file.path(wd, 'Correlation_plots'), showWarnings = FALSE)
    setwd(file.path(wd, 'Correlation_plots'))
    df <- na.omit(df)
    pdf(file="cor_biomass_plants.pdf",height=8,width=8)
    plot(df[,8:ncol(df)]) ; dev.off() 
    write.csv(cor(df[,8:ncol(df)]),file="cor_biomass_plants.csv"); setwd(wd)
df <- read.csv("phys_plants.csv",header=T)
    setwd(file.path(wd, 'Correlation_plots'))
    df <- na.omit(df)
    pdf(file="cor_phys_plants.pdf",height=30,width=30)
    plot(df[,10:ncol(df)]) ; dev.off()
    write.csv(cor(df[,10:ncol(df)]),file="cor_phys_plants.csv"); setwd(wd)
df <- read.csv("all_plants.csv",header=T)
    setwd(file.path(wd, 'Correlation_plots'))
    df <- na.omit(df)
    pdf(file="cor_growth_plants.pdf",height=15,width=15)
    plot(df[,6:ncol(df)]) ; dev.off()
    write.csv(cor(df[,6:ncol(df)]),file="cor_all_plants.csv"); setwd(wd)
df <- read.csv("recovery_plants.csv",header=T)
    setwd(file.path(wd, 'Correlation_plots'))
    df <- na.omit(df)
    pdf(file="cor_recovery_plants.pdf",height=8,width=8)
    plot(df[,10:ncol(df)]) ; dev.off()
    write.csv(cor(df[,10:ncol(df)]),file="cor_recovery_plants.csv"); setwd(wd)
setwd(wd)
    
###########################################################################################
## PRINCIPAL COMPONENTS WERE USED TO DETERMINE TRAITS OF INTEREST
###########################################################################################

getprcomps <- function(df=read.csv(file="biomass_plants.csv",header=T), limits=c(8:ncol(df))){
  df <- na.omit(df)
  df.responsevars <- df[,limits]
  means <- apply(as.matrix(df.responsevars),2,mean)
  normalized.df <- matrix(nrow=nrow(df.responsevars), ncol=ncol(df.responsevars))
  for(i in 1:ncol(df.responsevars)){ ## RESPONSES MUST BE SCALED
    coldata <- df.responsevars[,i] / means[i]
    normalized.df[,i] <- coldata
  }
  colnames(normalized.df) <- names(df.responsevars)
  pr.cmp2 <- prcomp(normalized.df)
  pr.cmp <- princomp(normalized.df)
  print(summary(pr.cmp2))
  return(pr.cmp2)
}
getprcomps( df=read.csv(file = "biomass_plants.csv",header=T), limits = c(8:23) )
getprcomps( df=read.csv(file = "phys_plants.csv",header=T), limits = c(10:59) )
getprcomps( df=read.csv(file = "all_plants.csv",header=T), limits = c(6:25) )
getprcomps( df=read.csv(file = "recovery_plants.csv",header=T), limits = c(10:12,14,15,17:26) ) #excluding flowering params used in later models 


  ###########################################################################################
## MODEL FOR MORPHOLOGY, PHYSIOLOGICAL TRAITS OF INTEREST
###########################################################################################

  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
        data{
        int<lower=0> N;
        int<lower=0> J;
        vector[N] y;
        matrix[N,J] X;
        }
        parameters{
        vector[J] b;
        real<lower=0> sigma;
        }
        transformed parameters{
        vector[N] mu;
        mu=X*b;
        }
        model{
        sigma ~ cauchy(0,5);
        b ~ normal(0,1000000);
        y ~ normal(mu, sigma);  //likelihood
        }
        generated quantities{
        vector[N] e_y;
        vector[15] Y;
        Y[1] = b[1];
        Y[2] = b[1] + b[2];
        Y[3] = b[1] + b[3];
        Y[4] = b[1] + b[4];
        Y[5] = b[1] + b[5];
        Y[6] = b[1] + b[6];
        Y[7] = b[1] + b[2] + b[6] + b[8];
        Y[8] = b[1] + b[3] + b[6] + b[9];
        Y[9] = b[1] + b[4] + b[6] + b[10];
        Y[10] = b[1] + b[5] + b[6] + b[11];
        Y[11] = b[1] + b[7];
        Y[12] = b[1] + b[2] + b[7] + b[12];
        Y[13] = b[1] + b[3] + b[7] + b[13];
        Y[14] = b[1] + b[4] + b[7] + b[14];
        Y[15] = b[1] + b[5] + b[7] + b[15];
        e_y = y - mu;
        }
        "
  comp <- stan_model(model_code = Stan.model) ##COMPILE
  ##LOAD DATA - CHANGE TO NAME FILES DIFFERENTLY IF DESIRED
  trait.model <- function(infile="biomass_plants.csv",
                        pdfname="rxn_norm_Bv.pdf",
                        pdfname_normtest="normtest_Bv.pdf",
                        response.label=expression(paste("Aboveground Biomass (g)")),
                        response='Bv'){
    setwd(wd)
    df <- read.csv(infile,header=T) ; response.1 <- df[,c(response)]
    df.filt <- df[!(response.1 = is.na(response.1)),] ; nrow(df.filt) ## ALLOWS N TO VARY
    dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df.filt) ## MODEL MATRIX FOR STAN
    response.var <- df.filt[,c(response)] ## NEW DEPENDENT VAR WITH NAs REMOVED
    model.components <- list( 'N' = nrow(df.filt), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
    iter <- 10000
    ##SAMPLE
    fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
    summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
    ## OUTPUT DATA
    dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
    setwd(file.path(wd, 'Parameter_estimates'))
    ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
    write.csv(ests, file=paste(response,"_parameters_",infile,sep="") )
    setwd(wd)
      ##RUN NORMALITY TESTS
      sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
      rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
      dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
      setwd(file.path(wd, 'Normality_tests'))
      pdf(file=pdfname_normtest,height=7,width=7) ##WRITES NORMALITY PLOTS TO DIR 'NORMALITY_TESTS'
    test_norm(rdf$mean); dev.off()
      ##PLOTTING REACTION NORMS
      geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
      trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
      facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2) ## ALLOWS SEPERATION OF X AXIS FOR SAT'D TREATMENT
      plot.params <- cbind(sdf,geno,trt,facet)
    gg <- ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
      geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
      xlab("% VWC") +
      ylab(response.label) +
      geom_line() +
      geom_point(size=1) +
      theme_classic() +
      theme(legend.position="bottom") +
      scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
      scale_x_continuous(breaks=c(10,15,20,25,35), labels=c("10","15","20","25","Sat'd")) +
      facet_grid(.~facet, scales="free",space="free") +
      theme(strip.text = element_blank()) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = -7, b = 0, l = 0)))
    dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
    setwd(file.path(wd, 'Reaction_Norms'))
    ggsave(filename=pdfname,plot=gg,height=4,width=4)
    setwd(wd)
    print(summ_fit[grep("geno", rownames(summ_fit)), ])
    return(gg)
  }
  ##RUN MODEL FOR MORPHOLOGICAL TRAITS
  p_1 <- trait.model(infile="biomass_plants.csv",
              pdfname="rxn_norm_Bv.pdf",
              pdfname_normtest="normtest_Bv.pdf",
              response.label=expression(paste("Aboveground Biomass (g)")),
              response="Bv") + theme(legend.position = "none", axis.title.x = element_blank())
  p_2 <- trait.model(infile="biomass_plants.csv",
              pdfname="rxn_norm_Rl.pdf",
              pdfname_normtest="normtest_Rl.pdf",
              response.label=expression(paste("Root length (mm)")),
              response="Rl") + theme(legend.position = "none", axis.title.x = element_blank())
  p_3 <- trait.model(infile="biomass_plants.csv",
              pdfname="rxn_norm_Rsa.pdf",
              pdfname_normtest="normtest_Rsa.pdf",
              response.label=expression(paste("Root surface area (",mm^2,")")),
              response="Rsa") + theme(legend.position = "none")
  p_4 <- trait.model(infile="biomass_plants.csv",
              pdfname="rxn_norm_Rt.pdf",
              pdfname_normtest="normtest_Rt.pdf",
              response.label=expression(paste("Root tips")),
              response="Rt") + theme(legend.position = "none")
  p_5 <- trait.model(infile="biomass_plants.csv",
                     pdfname="rxn_norm_Rv.pdf",
                     pdfname_normtest="normtest_Rv.pdf",
                     response.label=expression(paste("Root volume (",mm^3,")")),
                     response="Rv") + theme(legend.position = "none")
  p_6 <- trait.model(infile="biomass_plants.csv",
                     pdfname="rxn_norm_LA.pdf",
                     pdfname_normtest="normtest_LA.pdf",
                     response.label=expression(paste("Leaf area (",mm^2,")")),
                     response="LA") + theme(legend.position = "none")
  p_7 <- trait.model(infile="biomass_plants.csv",
                     pdfname="rxn_norm_rsa_la.pdf",
                     pdfname_normtest="normtest_rsa_la.pdf",
                     response.label=expression(paste("Root to leaf area ratio (",mm^2," ",mm^-2,")")),
                     response="Rsa_la") + theme(legend.position = "none")
  p_8 <- trait.model(infile="biomass_plants.csv",
                     pdfname="rxn_norm_srl.pdf",
                     pdfname_normtest="normtest_srl.pdf",
                     response.label=expression(paste("Specific root length (",mm," ",g^-1,")")),
                     response="srl") + theme(legend.position = "none")
  ## CONSTRUCT FIGURE FROM EIGHT INDIVIDUAL FIGURES ABOVE
    grid=plot_grid( p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8, align = 'vh', labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f) ", " (g)", "(h)"), label_size=15, hjust = -3, nrow = 2, rel_widths=c(1,1) )
    grobs=ggplotGrob(p_1 + theme(legend.position="bottom",legend.box.just = "left"))$grobs
    legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    pdf(file="FIGURE_2.pdf",height=5,width=10.5)
    plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off()
  ##RUN MODEL FOR PHYSIOLOGICAL TRAITS
  p_9 <- trait.model(infile="phys_plants.csv",
              pdfname="rxn_norm_Anet6_0913.pdf",
              pdfname_normtest="normtest_Anet6_0913.pdf",
              response.label=expression(paste(A[net]," week 6")),
              response="Anet6_0913") + theme(legend.position = "none", axis.title.x = element_blank())
  p_10 <- trait.model(infile="phys_plants.csv",
              pdfname="rxn_norm_E5_0907.pdf",
              pdfname_normtest="normtest_E5_0907.pdf",
              response.label=expression(paste("Evap. rate week 5")),
              response="E5_0907") + theme(legend.position = "none", axis.title.x = element_blank())
  p_11 <- trait.model(infile="phys_plants.csv",
              pdfname="rxn_norm_gs5_0907.pdf",
              pdfname_normtest="normtest_gs5_0907.pdf",
              response.label=expression(paste(g[s]," week 5")),
              response="gs5_0907") + theme(legend.position = "none", axis.title.x = element_blank())
  p_12 <- trait.model(infile="phys_plants.csv",
              pdfname="rxn_norm_gs6_0913.pdf",
              pdfname_normtest="normtest_gs6_0913.pdf",
              response.label=expression(paste(g[s]," week 6")),
              response="gs6_0913") + theme(legend.position = "none")
  p_13 <- trait.model(infile="phys_plants.csv",
              pdfname="rxn_norm_gs10_1011.pdf",
              pdfname_normtest="normtest_gs10_1011.pdf",
              response.label=expression(paste(g[s]," week 10")),
              response="gs10_1011") + theme(legend.position = "none")
  ##DOESN'T MAKE AN OVERALL FIGURE YET - HAVE TO ADJUST gs AT WEEK 7 (SEE BELOW)

###########################################################################################
## MODEL DIFFERENT FOR THE FOLLOWING; NEEDED TO TRANSFORM gs AT WEEK 7 TO MEET NORMALITY ASSUMPTIONS
###########################################################################################

  #LOAD DATA
  trait.model.modified <- function(infile="phys_plants.csv",
                        pdfname="rxn_norm_gs7_0920.pdf",
                        pdfname_normtest="normtest_gs7_0920.pdf",
                        response.label=expression(paste(g[s]," week 7")),
                        response='gs7_0920'){
    df <- read.csv(infile,header=T) ; response.1 <- df[,c(response)]
    df.filt <- df[!(response.1 = is.na(response.1)),] ; nrow(df.filt)     ## ALLOWS N TO VARY
    dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df.filt)     ## MODEL MATRIX FOR STAN
    response.var <- df.filt[,c(response)]                                 ## NEW DEPENDENT VAR WITH NAs REMOVED
  model.components <- list( 'N' = nrow(df.filt), 'y' = log(sqrt(response.var)), 'X' = dm, 'J' = ncol(dm) )
  iter <- 10000
  ##SAMPLE
  fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
  summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
  ## OUTPUT DATA
  dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
  setwd(file.path(wd, 'Parameter_estimates'))
  ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
  write.csv(ests, file=paste(response,"_parameters_",infile,sep="") )
  setwd(wd)
    ##RUN NORMALITY TESTS
    sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
    rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
    dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
    setwd(file.path(wd, 'Normality_tests'))
    pdf(file=pdfname_normtest,height=7,width=7)
  test_norm(rdf$mean); dev.off() ; setwd(wd)
    ##PLOT REACTION NORMS
    geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
    trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
    plot.params <- cbind(sdf,geno,trt,facet)
    plot.params$mean <- (exp(plot.params$mean))^2 ##BACK TRANSFORM
    plot.params$`2.5%` <- (exp(plot.params$`2.5%`))^2
    plot.params$`97.5%` <- (exp(plot.params$`97.5%`))^2
  gg <- ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
    xlab("% VWC") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
    scale_x_continuous(breaks=c(10,15,20,25,35), labels=c("10","15","20","25","Sat'd")) +
    facet_grid(.~facet, scales="free",space="free") +
    theme(strip.text = element_blank())+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = -7, b = 0, l = 0)))
  dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
  setwd(file.path(wd, 'Reaction_Norms'))
  ggsave(filename=pdfname,plot=gg,height=4,width=4)
  setwd(wd)
  return(gg)
  }
  ## RUN MODEL FOR FINAL PHYSIOLOGICAL TRAIT
  p_14 <- trait.model.modified() + theme(legend.position = "none")
  ## CONSTRUCT FIGURE FROM SIX INDIVIDUAL FIGURES ABOVE
    grid=plot_grid( p_10 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))),
                    p_11 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))),
                    p_9 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))),
                    p_12 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))),
                    p_14 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))),
                    p_13 + theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))), 
                    align = 'vh', labels = c(" (a)", "(b)", "(c)", "(d)", "(e)", "(f)"), label_size=15, hjust = -2.5, nrow = 2)
    grobs=ggplotGrob(p_1 + theme(legend.position="bottom",legend.box.just = "left"))$grobs
    legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    pdf(file="FIGURE_1.pdf",height=5,width=10.5)
    plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off()

###########################################################################################
## MODEL FOR HEIGHT OVER TIME (GROWTH RATE) 
###########################################################################################

  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
      data{
      int<lower=0> N;
      int<lower=0> J;
      vector[N] y;
      matrix[N,J] X;
      }
      parameters{
      vector[J] b;
      real<lower=0> sigma;
      }
      transformed parameters{
      vector[N] mu;
      mu=X*b;
      }
      model{
      sigma ~ cauchy(0,5);
      b ~ normal(0,1000000);
      y ~ normal(mu, sigma);  //likelihood
      }
      generated quantities{
      vector[N] e_y;
      vector[15] Y;
      Y[1] = b[1];
      Y[2] = b[1] + b[2];
      Y[3] = b[1] + b[3];
      Y[4] = b[1] + b[4];
      Y[5] = b[1] + b[5];
      Y[6] = b[1] + b[6];
      Y[7] = b[1] + b[2] + b[6] + b[8];
      Y[8] = b[1] + b[3] + b[6] + b[9];
      Y[9] = b[1] + b[4] + b[6] + b[10];
      Y[10] = b[1] + b[5] + b[6] + b[11];
      Y[11] = b[1] + b[7];
      Y[12] = b[1] + b[2] + b[7] + b[12];
      Y[13] = b[1] + b[3] + b[7] + b[13];
      Y[14] = b[1] + b[4] + b[7] + b[14];
      Y[15] = b[1] + b[5] + b[7] + b[15];
      e_y = y - mu;
      }
      "
  pdfname="FIGURE_3.pdf"
  comp <- stan_model(model_code = Stan.model)
  setwd(wd)
  df <- read.csv("all_plants.csv",header=T)
  dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df)     ## MODEL MATRIX FOR STAN
  height.model <- function(response.var=df$H2_0821,pdfname_normtest="normtest_H2.pdf",date_ID="2",response="H2_0821"){
        model.components <- list( 'N' = nrow(df), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
        iter <- 10000
        ##SAMPLE
        fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3, control=list(adapt_delta=0.99) )
        summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
        ## OUTPUT DATA
        dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
        setwd(file.path(wd, 'Parameter_estimates'))
        ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
        write.csv(ests, file=paste(response,"_parameters_all_plants.csv",sep="") )
        setwd(wd)
        sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
        rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
        dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
        setwd(file.path(wd, 'Normality_tests'))
        pdf(file=pdfname_normtest,height=7,width=7); test_norm(rdf$mean); dev.off() ; setwd(wd) ##PLOT NORMALITY TESTS
        geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
        trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
        facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2) ## SATURATED TREATMENT ON SEPERATED X AXIS
        date <- c(rep(date_ID,15))
        sdf <- cbind(sdf,geno,trt,facet,date)
        return(sdf)
  }
  ##RUN MODELS, PLOT COMES LATER
  df1 <- height.model()
  df2 <- height.model(df$H4_0904,"normtest_H4.pdf","4",response="H4_0904")
  df3 <- height.model(df$H6_0918,"normtest_H6.pdf","6",response="H6_0918")
  df4 <- height.model(df$H7_0926,"normtest_H7.pdf","7",response="H7_0926")
  df5 <- height.model(df$H8_1002,"normtest_H8.pdf","8",response="H8_1002")
  df6 <- height.model(df$H9_1009,"normtest_H9.pdf","9",response="H9_1009")
  df7 <- height.model(df$H10_1015,"normtest_H10.pdf","10",response="H10_1015")

   height.model2 <- function(response.var=df$D4_0904,pdfname_normtest="normtest_D4.pdf",date_ID="week4",response="D4_0904"){
        model.components <- list( 'N' = nrow(df), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
        iter <- 100000
        ##SAMPLE
        fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3, control=list(adapt_delta=0.99) )
        summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
        ## OUTPUT DATA
        dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
        setwd(file.path(wd, 'Parameter_estimates'))
        ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
        write.csv(ests, file=paste(response,"_parameters_all_plants.csv",sep="") )
        setwd(wd)
        sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
        rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
        dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
        setwd(file.path(wd, 'Normality_tests'))
        pdf(file=pdfname_normtest,height=7,width=7); test_norm(rdf$mean); dev.off() ; setwd(wd) ##PLOT NORMALITY TESTS
        geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
        trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
        facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2) ## SATURATED TREATMENT ON SEPERATED X AXIS
        date <- c(rep(date_ID,15))
        sdf <- cbind(sdf,geno,trt,facet,date)
        return(sdf)
  }

  df9 <- height.model2(df$D4_0904,"normtest_D4.pdf","4",response="D4_0904")
  df10 <- height.model2(df$D6_0918,"normtest_D6.pdf","6",response="D6_0918")
  df11 <- height.model2(df$D7_0926,"normtest_D7.pdf","7",response="D7_0926")
  df12 <- height.model2(df$D8_1002,"normtest_D8.pdf","8",response="D8_1002")
  df13 <- height.model2(df$D9_1009,"normtest_D9.pdf","9",response="D9_1009")
  df14 <- height.model2(df$D10_1015,"normtest_D10.pdf","10",response="D10_1015")
  
  height.data <- rbind(df1,df2,df3,df4,df5,df6,df7)
  height.data <- height.data[,c(1,4,8,11,12,13,14)]
  response.label <- expression(paste("Height (cm)"))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  height.data$date <- as.numeric(as.character(height.data$date))
  ##PLOT HEIGHT OVER TIME
  htplot <- ggplot( height.data, aes(x=date,y=mean,color=as.factor(trt))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(2,4,6,7,8,9,10), labels=c("2","4","6","7","8","9","10")) +
    facet_grid(.~geno)
  
  delta.data <- rbind(df9,df10,df11,df12,df13,df14)
  delta.data <- delta.data[,c(1,4,8,11,12,13,14)]
  response.label <- expression(paste("Relative growth rate (ln(cm)"~week^{-1}~")"))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  delta.data$date <- as.numeric(as.character(delta.data$date))
  ##PLOT HEIGHT OVER TIME
  deltaplot <- ggplot( delta.data, aes(x=date,y=mean,color=as.factor(trt))) +
    geom_hline(yintercept=0, color="grey") +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(4,6,7,8,9,10), labels=c("4","6","7","8","9","10")) +
    facet_grid(.~geno)
  
  grid=plot_grid( htplot + theme(legend.position="none"),
                  deltaplot + theme(legend.position="none"), 
                  align = 'v', labels = c("(a)", "(b)"), label_size=15, hjust = -3.5, vjust = 4, nrow = 2)
  grobs=ggplotGrob(htplot + theme(legend.position="bottom",legend.box.just = "left"))$grobs
  legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  pdf(file=pdfname,height=9,width=10.5)
  plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off()

###########################################################################################
## MODELS FOR RECOVERY - Flowering first, different distribution
## NOTE: run theta parameter seperately for genotype-level only flowering 
###########################################################################################
  
  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
        data{
        int<lower=0> N; //total number of observations
        int<lower=0> J; //number of matrix cols/regression parameters
        int<lower=0> G; // number of genotypes
        vector[N] y; //response variable
        matrix[N,J] X; //design matrix
        int geno[N]; //genotype index
        int<lower=0> w[N]; //Did flower or not
        }
        parameters{
        vector[J] b;
        real<lower=0> sigma;
        vector[G] flowerrate;
        }
        transformed parameters{
        vector[N] mu_derived;
        vector[N] mu;
        vector[N] theta;
        vector<lower=0>[N] alpha;
        vector<lower=0>[N] Beta;
        for(n in 1:N){
        theta[n] =  1 ./ (1 + exp( - flowerrate[geno[n]] ) ) ;
        }
        mu =  exp(X*b) ;
        for(n in 1:N){
        mu_derived[n] = (mu[n] + 1) / theta[n] ;
        }
        alpha = mu_derived .* mu_derived / sigma;
        Beta = mu_derived / sigma;
        }
        model{
        sigma ~ cauchy(0,5);
        Beta ~ cauchy(0,5);
        b ~ normal(0,3);
        flowerrate ~ normal(0,10);
        y ~ gamma(alpha,Beta);
        w ~ bernoulli(theta);
        }
        generated quantities{
        vector[N] e_y;
        vector[15] Y;
          Y[1] = exp(b[1]);
          Y[2] = exp(b[1] + b[2]);
          Y[3] = exp(b[1] + b[3]);
          Y[4] = exp(b[1] + b[4]);
          Y[5] = exp(b[1] + b[5]);
          Y[6] = exp(b[1] + b[6]);
          Y[7] = exp(b[1] + b[2] + b[6] + b[8]);
          Y[8] = exp(b[1] + b[3] + b[6] + b[9]);
          Y[9] = exp(b[1] + b[4] + b[6] + b[10]);
          Y[10] = exp(b[1] + b[5] + b[6] + b[11]);
          Y[11] = exp(b[1] + b[7]);
          Y[12] = exp(b[1] + b[2] + b[7] + b[12]);
          Y[13] = exp(b[1] + b[3] + b[7] + b[13]);
          Y[14] = exp(b[1] + b[4] + b[7] + b[14]);
          Y[15] = exp(b[1] + b[5] + b[7] + b[15]);  
          e_y = y - mu;
        }
        
        "
  comp <- stan_model(model_code = Stan.model)
  setwd(wd)
  ##LOAD AND SPECIFY DATA
  infile="recovery_plants.csv"
  pdfname="rxn_norm_recovery_Bf.pdf"
  pdfname_normtest="normtest_recovery_Bf.pdf"
  response.label=expression(paste("Flowering biomass (g)"))
  response='Bf'
    df <- read.csv(infile,header=T) ; response.1 <- df[,c(response)]
    df.filt <- df[!(response.1 = is.na(response.1)),] ; nrow(df.filt)     ## ALLOWS N TO VARY
    dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df.filt)     ## MODEL MATRIX FOR STAN
    response.var <- df.filt[,c(response)]                                 ## NEW DEPENDENT VAR WITH NAs REMOVED
  ##SPECIFY MODEL COMPONENTS
  model.components <- list( 'N' = nrow(df.filt), 
                            'y' = response.var, 
                            'X' = dm, 
                            'J' = ncol(dm),
                            'w' = df.filt$did_flwr, ## WHETHER PLANTS FLOWERED OR NOT, 0 OR 1
                            'geno' = df.filt$geno,
                            'G' = 3)
  iter <- 10000                                                         ## MAY TAKE UP TO 60s ; ADJUST ACCORDINGLY
  ##SAMPLING
  fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
    summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary) ; summ_fit
    ## OUTPUT DATA
    dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
    setwd(file.path(wd, 'Parameter_estimates'))
    ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
    write.csv(ests, file=paste(response,"_parameters_",infile,sep="") )
    setwd(wd)
    thetas <- summ_fit[grep("theta", rownames(summ_fit)), ]
    sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
    rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
    dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
    setwd(file.path(wd, 'Normality_tests'))
    pdf(file=pdfname_normtest,height=7,width=7)
    test_norm(rdf$mean); dev.off() ; setwd(wd) ##PLOT NORMALITY PLOTS
    ##PLOT REACTION NORMS
    geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
    trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
    plot.params <- cbind(sdf,geno,trt,facet)
  gg1 <- ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
    xlab("% VWC") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
    scale_x_continuous(breaks=c(10,15,20,25,35), labels=c("10","15","20","25","Sat'd")) +
    facet_grid(.~facet, scales="free",space="free") +
    theme(strip.text = element_blank())
  gg1
  dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
  setwd(file.path(wd, 'Reaction_Norms'))
  ggsave(filename=pdfname,plot=gg1,height=4,width=4) ;  setwd(wd)
  
###########################################################################################
## MODELS FOR RECOVERY - Rhizome recovery
## NOTE: run theta parameter seperately for genotype-level only recovery
###########################################################################################

  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
        data{
  int<lower=0> N; //total number of observations
  int<lower=0> J; //number of matrix cols/regression parameters
  int<lower=0> G; // number of genotypes
  vector[N] y; //response variable
  matrix[N,J] X; //design matrix
  int geno[N]; //genotype index
  int<lower=0> w[N]; //Did rhizome recover or not
  }
  parameters{
  vector[J] b;
  real<lower=0> sigma;
  vector[G] rhizparam; //rhizome parameter
  }
  transformed parameters{
  vector[N] mu_derived;
  vector[N] mu;
  vector[N] theta;
  vector<lower=0>[N] alpha;
  vector<lower=0>[N] Beta;
  for(n in 1:N){
  theta[n] =  1 ./ (1 + exp( - rhizparam[geno[n]] ) ) ;
  }
  mu =  exp(X*b) ;
  for(n in 1:N){
  mu_derived[n] = (mu[n] + 1) / theta[n] ;
  }
  alpha = mu_derived .* mu_derived / sigma;
  Beta = mu_derived / sigma;
  }
  model{
  sigma ~ cauchy(0,5);
  Beta ~ cauchy(0,5);
  b ~ normal(0,3);
  rhizparam ~ normal(0,10);
  y ~ gamma(alpha,Beta);
  w ~ bernoulli(theta);
  }
  generated quantities{
  vector[N] e_y;
  vector[15] Y;
  Y[1] = exp(b[1]);
  Y[2] = exp(b[1] + b[2]);
  Y[3] = exp(b[1] + b[3]);
  Y[4] = exp(b[1] + b[4]);
  Y[5] = exp(b[1] + b[5]);
  Y[6] = exp(b[1] + b[6]);
  Y[7] = exp(b[1] + b[2] + b[6] + b[8]);
  Y[8] = exp(b[1] + b[3] + b[6] + b[9]);
  Y[9] = exp(b[1] + b[4] + b[6] + b[10]);
  Y[10] = exp(b[1] + b[5] + b[6] + b[11]);
  Y[11] = exp(b[1] + b[7]);
  Y[12] = exp(b[1] + b[2] + b[7] + b[12]);
  Y[13] = exp(b[1] + b[3] + b[7] + b[13]);
  Y[14] = exp(b[1] + b[4] + b[7] + b[14]);
  Y[15] = exp(b[1] + b[5] + b[7] + b[15]);  
  e_y = y - mu;
  }
  
  "
  comp <- stan_model(model_code = Stan.model)
  setwd(wd)
  ##LOAD AND SPECIFY DATA
  infile="recovery_plants.csv"
  pdfname="rxn_norm_recovery_Brh.pdf"
  pdfname_normtest="normtest_recovery_Brh.pdf"
  response.label=expression(paste("Rhizome recovery (cm)"))
  response='Brh'
  df <- read.csv(infile,header=T) ; response.1 <- df[,c(response)]
  df.filt <- df[!(response.1 = is.na(response.1)),] ; nrow(df.filt)     ## ALLOWS N TO VARY
  dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df.filt)     ## MODEL MATRIX FOR STAN
  response.var <- df.filt[,c(response)]                                 ## NEW DEPENDENT VAR WITH NAs REMOVED
  ##SPECIFY MODEL COMPONENTS
  model.components <- list( 'N' = nrow(df.filt), 
                            'y' = response.var, 
                            'X' = dm, 
                            'J' = ncol(dm),
                            'w' = df.filt$did_rh, ## WHETHER PLANTS RECOVERED OR NOT, 0 OR 1
                            'geno' = df.filt$geno,
                            'G' = 3)
  iter <- 10000  ## MAY TAKE UP TO 60s ; ADJUST ACCORDINGLY
  ##SAMPLING
  fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
  summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary) ; summ_fit
  ## OUTPUT DATA
  dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
  setwd(file.path(wd, 'Parameter_estimates'))
  ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
  write.csv(ests, file=paste(response,"_parameters_",infile,sep="") )
  setwd(wd)
  thetas <- summ_fit[grep("theta", rownames(summ_fit)), ]
  sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
  rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
  dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
  setwd(file.path(wd, 'Normality_tests'))
  pdf(file=pdfname_normtest,height=7,width=7)
  test_norm(rdf$mean); dev.off() ; setwd(wd) ##PLOT NORMALITY PLOTS
  ##PLOT REACTION NORMS
  geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
  trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
  facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
  plot.params <- cbind(sdf,geno,trt,facet)
  gg2 <- ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
    xlab("% VWC") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
    scale_x_continuous(breaks=c(10,15,20,25,35), labels=c("10","15","20","25","Sat'd")) +
    facet_grid(.~facet, scales="free",space="free") +
    theme(strip.text = element_blank())
  gg2
  dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
  setwd(file.path(wd, 'Reaction_Norms'))
  ggsave(filename=pdfname,plot=gg1,height=4,width=4) ;  setwd(wd)
  
###########################################################################################
## MODEL FOR RECOVERY GROWTH RATES
###########################################################################################
  
  Stan.model <- "
        data{
        int<lower=0> N;
        int<lower=0> J;
        vector[N] y;
        matrix[N,J] X;
        }
        parameters{
        vector[J] b;
        real<lower=0> sigma;
        }
        transformed parameters{
        vector[N] mu;
        mu=X*b;
        }
        model{
        sigma ~ cauchy(0,5);
        b ~ normal(0,1000000);
        y ~ normal(mu, sigma);  //likelihood
        }
        generated quantities{
        vector[N] e_y;
        vector[15] Y;
        Y[1] = b[1];
        Y[2] = b[1] + b[2];
        Y[3] = b[1] + b[3];
        Y[4] = b[1] + b[4];
        Y[5] = b[1] + b[5];
        Y[6] = b[1] + b[6];
        Y[7] = b[1] + b[2] + b[6] + b[8];
        Y[8] = b[1] + b[3] + b[6] + b[9];
        Y[9] = b[1] + b[4] + b[6] + b[10];
        Y[10] = b[1] + b[5] + b[6] + b[11];
        Y[11] = b[1] + b[7];
        Y[12] = b[1] + b[2] + b[7] + b[12];
        Y[13] = b[1] + b[3] + b[7] + b[13];
        Y[14] = b[1] + b[4] + b[7] + b[14];
        Y[15] = b[1] + b[5] + b[7] + b[15];
        e_y = y - mu;
        }
        "
  comp <- stan_model(model_code = Stan.model)
  ##LOAD AND SPECIFY DATA
  trait.model <- function(infile="recovery_plants.csv",
                          pdfname="rxn_norm_D12.pdf",
                          pdfname_normtest="normtest_D12.pdf",
                          response.label=expression(paste("Growth rate (cm/week)")),
                          response='D12_1104'){
    setwd(wd)
      df <- read.csv(infile,header=T) ; response.1 <- df[,c(response)]
      df.filt <- df[!(response.1 = is.na(response.1)),] ; nrow(df.filt) ## ALLOWS N TO VARY
      dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df.filt) ## MODEL MATRIX FOR STAN
      response.var <- df.filt[,c(response)] ## NEW DEPENDENT VAR WITH NAs REMOVED
    model.components <- list( 'N' = nrow(df.filt), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
    iter <- 10000
    ##SAMPLING
    fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
      summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
      ## OUTPUT DATA
      dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
      setwd(file.path(wd, 'Parameter_estimates'))
      ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
      write.csv(ests, file=paste(response,"_parameters_",infile,sep="") )
      setwd(wd)
      sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
      rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
      dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
      setwd(file.path(wd, 'Normality_tests'))
      pdf(file=pdfname_normtest,height=7,width=7)
      test_norm(rdf$mean); dev.off() ##NORMALITY PLOTS
    ##PLOT REACTION NORMS
    geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
    trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2)
    plot.params <- cbind(sdf,geno,trt,facet)
    gg <- ggplot( plot.params, aes(x=trt,y=mean,color=as.factor(geno))) +
      geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=1) +
      xlab("% VWC") +
      ylab(response.label) +
      geom_line() +
      geom_point(size=1) +
      theme_classic() +
      theme(legend.position="bottom") +
      scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
      scale_x_continuous(breaks=c(10,15,20,25,35), labels=c("10","15","20","25","Sat'd")) +
      facet_grid(.~facet, scales="free",space="free") +
      theme(strip.text = element_blank())
    dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
    setwd(file.path(wd, 'Reaction_Norms'))
    ggsave(filename=pdfname,plot=gg,height=4,width=4)
    setwd(wd)
    return(gg)
  }
  ##RUN MODELS AND GENERATE SEPARATE GRAPHS
  p_1 <- gg1 + theme(legend.position = "none", axis.title.x = element_blank())
  p_2 <- gg2 + theme(legend.position = "none")
  
  ## !! excluding growth rates from final figure
  p_3 <- trait.model(infile="recovery_plants.csv",
                     pdfname="rxn_norm_D12.pdf",
                     pdfname_normtest="normtest_D12.pdf",
                     response.label=expression(paste("Growth rate, week 12 (cm/week)")),
                     response='D12_1104') + theme(legend.position = "none", axis.title.x = element_blank())
  p_4 <- trait.model(infile="recovery_plants.csv",
                     pdfname="rxn_norm_D14.pdf",
                     pdfname_normtest="normtest_D14.pdf",
                     response.label=expression(paste("Growth rate, week 14 (cm/week)")),
                     response="D14_1117") + theme(legend.position = "none", axis.title.x = element_blank())
  
  # ## CONSTRUCT FIGURE FROM INDIVIDUAL FIGURES ABOVE
  grid=plot_grid( p_1,p_2,align = 'v', labels = c("(a)", "(b)"), label_size=15, hjust = -2, nrow = 2, rel_heights = c(1,1) )
  grobs=ggplotGrob(p_1 + theme(legend.position="bottom",legend.box.just = "left"))$grobs
  legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  pdf(file="FIGURE_4.pdf",height=7.2,width=3.5)
  plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off() ; dev.off() ; dev.off()
  
###########################################################################################
## MODELS FOR TRAIT CORRELATIONS
## Not included in final manuscript
###########################################################################################
  # 
  # ##WORKING WITH AVERAGES TO DETERMINE CORRELATIONS
  # setwd(wd)
  # df1=read.csv(file="biomass_plants.csv",header=T)
  # df2=read.csv(file="phys_plants.csv",header=T) ; names(df2)
  # a=aggregate(df1[, "Bv"], list(df1$geno,df1$trt), mean, na.rm=TRUE)
  # b=aggregate(df1[, "Br"], list(df1$geno,df1$trt), mean, na.rm=TRUE)
  # total_bio= a$x + b$x ##ADD ABOVE AND BELOWGROUND BIOMASS
  # ##MEAN PHOTOSYNTHETIC RATE
  # photo=aggregate(df2[, "uAnet"], list(df2$geno,df2$trt), mean, na.rm=TRUE)
  # genos=c(11,2,5,11,2,5,11,2,5,11,2,5,11,2,5) ##SPECIFY GENOTYPES FOR PLOTTING
  # plot.data1=cbind(total_bio,photo$x,genos) ; colnames(plot.data1) <- c("bio","anet","geno") ; plot.data1 <- as.data.frame(plot.data1)
  # p1 <- ggplot(plot.data1, aes(x=anet, y=bio,color=as.factor(geno)))+
  #   geom_point(size=2) +
  #   stat_smooth(method=lm,se=F,color="black") +
  #   theme_classic() +
  #   ylab("Total biomass (g)") +
  #   xlab(expression(Average~A[net]~(mu~mol~m^{-2}~s^{-1}))) +
  #   scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  #   theme(legend.position = "none")
  # ##WATER USE EFFICIENCY
  # iWUE=aggregate(df2[, "uWUEi"], list(df2$geno,df2$trt), mean, na.rm=TRUE)
  # genos=c(11,2,5,11,2,5,11,2,5,11,2,5,11,2,5) ##SPECIFY GENOTYPES FOR PLOTTING
  # plot.data2=cbind(total_bio,iWUE$x,genos) ; colnames(plot.data2) <- c("bio","iWUE","geno") ; plot.data2 <- as.data.frame(plot.data2)
  # p2 <- ggplot(plot.data2, aes(x=iWUE, y=bio,color=as.factor(geno)))+
  #   geom_point(size=2) +
  #   #stat_smooth(method=lm,se=F,color="black") +
  #   theme_classic() +
  #   ylab("Total biomass (g)") +
  #   xlab(expression(Average~iWUE~(mu~mol~mmol^{-1}))) +
  #   scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  #   theme(legend.position = "none")
  # ##PSII EFFICIENCY
  # fvfm=aggregate(df2[, "ufv"], list(df2$geno,df2$trt), mean, na.rm=TRUE)
  # genos=c(11,2,5,11,2,5,11,2,5,11,2,5,11,2,5) ##SPECIFY GENOTYPES FOR PLOTTING
  # plot.data3=cbind(total_bio,fvfm$x,genos) ; colnames(plot.data3) <- c("bio","fvfm","geno") ; plot.data3 <- as.data.frame(plot.data3)
  # p3 <- ggplot(plot.data3, aes(x=fvfm, y=bio,color=as.factor(geno)))+
  #   geom_point(size=2) +
  #   stat_smooth(method=lm,se=F,color="black") +
  #   theme_classic() +
  #   ylab("Total biomass (g)") +
  #   xlab(expression(Average~PSII~efficiency~(mu~Fv/Fm))) +
  #   scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
  #   theme(legend.position = "none")
  # ## CONSTRUCT FIGURE FROM INDIVIDUAL FIGURES ABOVE
  # grid=plot_grid( p1,p3,p2, align = 'h', labels = c("A", "B", "C"), label_size=30, hjust = -0.5, nrow = 1, rel_widths=c(1,1,1) )
  # grobs=ggplotGrob(p1 + theme(legend.position="bottom",legend.box.just = "left"))$grobs
  # legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  # pdf(file="FIGURE_6.pdf",height=3.5,width=9)
  # plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off() ; dev.off() ; dev.off()
  # ##RUN MODEL TO DETERMINE SLOPE OF THESE RELATIONSHIPS
  # Stan.model <- "
  #       data{
  #       int<lower=0> N;
  #       vector[N] y;
  #       vector[N] x;
  #       }
  #       parameters{
  #       real b0; //intercept
  #       real b1; //intercept
  #       real<lower=0> sigma;
  #       }
  #       transformed parameters{
  #       vector[N] mu;
  #       mu = b0 + b1*x;
  #       }
  #       model{
  #       sigma ~ cauchy(0,5);
  #       b0 ~ normal(0,10000);
  #       b1 ~ normal(0,10000);
  #       y ~ normal(mu, sigma);  //likelihood
  #       }
  #       generated quantities{
  #       vector[N] e_y;
  #       e_y = y - mu;
  #       }
  #       "
  # comp <- stan_model(model_code = Stan.model)
  # ##LOAD AND SPECIFY DATA
  # corr.model <- function(in_df=plot.data1,
  #                         pdfname_normtest="normtest_corr_photo_biomass.pdf",
  #                         predictor='anet'){
  #   setwd(wd)
  #   predictor.var <- in_df[,c(predictor)]
  #   model.components <- list( 'N' = nrow(in_df), 'y' = in_df$bio, 'x' = predictor.var)
  #   iter <- 10000 
  #   ##SAMPLING
  #   fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
  #   summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
  #   rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ] ##GATHER RESIDUALS
  #   dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
  #   setwd(file.path(wd, 'Normality_tests'))
  #   pdf(file=pdfname_normtest,height=7,width=7)
  #   test_norm(rdf$mean); dev.off() ##NORMALITY PLOTS
  #   return(summ_fit)
  # }
  #   
  # corr.model(in_df=plot.data1,
  #            pdfname_normtest="normtest_corr_photo_biomass.pdf",
  #            predictor='anet')
  # corr.model(in_df=plot.data3,
  #            pdfname_normtest="normtest_corr_fvfm_biomass.pdf",
  #            predictor='fvfm')
  # corr.model(in_df=plot.data2,
  #            pdfname_normtest="normtest_corr_iWUE_biomass.pdf",
  #            predictor='iWUE')
  
###########################################################################################
## MODEL FOR HEIGHT RECOVERY (GROWTH RATE) 
###########################################################################################
  
  setwd(wd)
  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
  data{
  int<lower=0> N;
  int<lower=0> J;
  vector[N] y;
  matrix[N,J] X;
  }
  parameters{
  vector[J] b;
  real<lower=0> sigma;
  }
  transformed parameters{
  vector[N] mu;
  mu=X*b;
  }
  model{
  sigma ~ cauchy(0,5);
  b ~ normal(0,1000000);
  y ~ normal(mu, sigma);  //likelihood
  }
  generated quantities{
  vector[N] e_y;
  vector[15] Y;
  Y[1] = b[1];
  Y[2] = b[1] + b[2];
  Y[3] = b[1] + b[3];
  Y[4] = b[1] + b[4];
  Y[5] = b[1] + b[5];
  Y[6] = b[1] + b[6];
  Y[7] = b[1] + b[2] + b[6] + b[8];
  Y[8] = b[1] + b[3] + b[6] + b[9];
  Y[9] = b[1] + b[4] + b[6] + b[10];
  Y[10] = b[1] + b[5] + b[6] + b[11];
  Y[11] = b[1] + b[7];
  Y[12] = b[1] + b[2] + b[7] + b[12];
  Y[13] = b[1] + b[3] + b[7] + b[13];
  Y[14] = b[1] + b[4] + b[7] + b[14];
  Y[15] = b[1] + b[5] + b[7] + b[15];
  e_y = y - mu;
  }
  "
  
pdfname="FIGURE_5.pdf"
  comp <- stan_model(model_code = Stan.model)
  df <- read.csv("recovery_plants.csv",header=T)
  dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df)     ## MODEL MATRIX FOR STAN
  height.model <- function(response.var=df$H11_1023,pdfname_normtest="normtest_H11.pdf",date_ID="11", response="H11_1023"){
    model.components <- list( 'N' = nrow(df), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
    iter <- 10000
    ##SAMPLE
    fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3, control=list(adapt_delta=0.99) )
    summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
    ## OUTPUT DATA
    dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
    setwd(file.path(wd, 'Parameter_estimates'))
    ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4,8) ])
    write.csv(ests, file=paste(response,"_parameters_recovery_plants.csv",sep="") )
    setwd(wd)
    sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
    rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
    dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
    setwd(file.path(wd, 'Normality_tests'))
    pdf(file=pdfname_normtest,height=7,width=7); test_norm(rdf$mean); dev.off() ; setwd(wd) ##PLOT NORMALITY TESTS
    geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
    trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2) ## SATURATED TREATMENT ON SEPERATED X AXIS
    date <- c(rep(date_ID,15))
    sdf <- cbind(sdf,geno,trt,facet,date)
    return(sdf)
  }
  ##RUN MODELS, PLOT COMES LATER
  df1 <- height.model()
  df2 <- height.model(df$H12_1104,"normtest_H12.pdf","12",response="H12_1104")
  df3 <- height.model(df$H13_1111,"normtest_H13.pdf","13",response="H13_1111")
  df4 <- height.model(df$H14_1117,"normtest_H14.pdf","14",response="H14_1117")
  df5 <- height.model(df$H15_1124,"normtest_H15.pdf","15",response="H15_1124")
  df6 <- height.model(df$D11_1023,"normtest_D11.pdf","11",response="D11_1023")
  df7 <- height.model(df$D12_1104,"normtest_D12.pdf","12",response="D12_1104")
  df8 <- height.model(df$D13_1111,"normtest_D13.pdf","13",response="D13_1111")
  df9 <- height.model(df$D14_1117,"normtest_D14.pdf","14",response="D14_1117")
  df10 <- height.model(df$D15_1124,"normtest_D15.pdf","15",response="D15_1124")
  height.data.recov <- rbind(df1,df2,df3,df4,df5)
  height.data.recov <- height.data.recov[,c(1,4,8,11,12,13,14)]
  response.label <- expression(paste("Height (cm)"))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  height.data.recov$date <- as.numeric(as.character(height.data.recov$date))
  ##PLOT HEIGHT OVER TIME
  recovp1 <- ggplot( height.data.recov, aes(x=date,y=mean,color=as.factor(trt))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week (recovery)") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(11,12,13,14,15), labels=c("11","12","13","14","15")) +
    facet_grid(.~geno)
  delta.data.recov <- rbind(df6,df7,df8,df9,df10)
  delta.data.recov <- delta.data.recov[,c(1,4,8,11,12,13,14)]
  delta.data.recov$date <- as.numeric(as.character(delta.data.recov$date))
  response.label <- expression(paste("Relative growth rate (ln(cm)"~week^{-1}~")"))
  recovp2 <- ggplot( delta.data.recov, aes(x=date,y=mean,color=as.factor(trt))) +
    geom_hline(yintercept=0, color="grey") +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week (recovery)") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(11,12,13,14,15), labels=c("11","12","13","14","15")) +
    facet_grid(.~geno)
  
  grid=plot_grid( recovp1 + theme(legend.position="none"),
                  recovp2 + theme(legend.position="none"), 
                  align = 'v', labels = c("(a)", "(b)"), label_size=15, hjust = -3.5, vjust = 4, nrow = 2)
  grobs=ggplotGrob(htplot + theme(legend.position="bottom",legend.box.just = "left"))$grobs
  legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
  pdf(file=pdfname,height=9,width=10.5)
  plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) ; dev.off()
  
  
###########################################################################################
## MODEL AND PLOTS FOR GROWTH RATE CORRELATIONS
## Not included in final manuscript
###########################################################################################
# 
#   setwd(wd)
#   df=read.csv(file="mechanism_cor.csv",header=T)
#   names(df)
#   df$geno <- as.factor(df$geno)
#   G11data <- df[(df$geno == 11),]
#   G2data <- df[(df$geno == 2),]
#   p <- ggplot(df, aes(x=avg_cond, y=Dht_0926, color=as.factor(geno), shape=as.factor(geno))) +
#     geom_point(aes(shape=geno)) +
#     scale_shape_manual(values=c(1,2,4),name="Genotype") +
#     scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
#     stat_smooth(data=G11data,method=lm,se=F,show.legend = F,color="#BBDF27FF") + #DON'T WANT TO SHOW NS DATA, USE CODE viridis(n=3, end=0.9) TO GET HEX CODE
#     theme_classic() +
#     theme(legend.position = "none") +
#     xlab(expression(mu~g[s]))+
#     ylab("growth rate week 7")
#   p
#   q <- ggplot(df, aes(x=avg_cond, y=Dht_1002, color=geno, shape=geno))+
#     geom_point(aes(shape=geno)) +
#     scale_shape_manual(values=c(1,2,4),name="Genotype") +
#     scale_color_viridis(discrete=T, end=0.9, name="Genotype") +
#     stat_smooth(data=G11data,method=lm,se=F,show.legend = F,color="#BBDF27FF") + #don't want to show a stat line for non-significant data
#     stat_smooth(data=G2data,method=lm,se=F,show.legend = F,color="#440154FF") + #don't want to show a stat line for non-significant data
#     theme_classic() +
#     theme(legend.position = "none") +
#     xlab(expression(mu~g[s]))+
#     ylab("growth rate week 8")
#   q
#   get_legend<-function(a.gplot){
#     tmp <- ggplot_gtable(ggplot_build(a.gplot))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     return(legend)}
#     mylegend<-get_legend(q + theme(legend.position = "right"))
#   figure7<- plot_grid(p,q,labels = "AUTO", align = 'h', nrow = 1, ncol=2,  label_size = 18)
#   finalplot_7 <- plot_grid(figure7, mylegend, ncol=2, nrow = 1, rel_widths = c(1, .15))
#   finalfigure7 <- "FIGURE_7.pdf"
#   pdf(file=finalfigure7,height=3,width=6)
#   finalplot_7
#   dev.off();dev.off()
#   ##RUN MODEL TO DETERMINE SLOPE OF THESE RELATIONSHIPS
#   Stan.model <- "
#   data{
#   int<lower=0> N;
#   vector[N] y;
#   vector[N] x;
#   int geno[N];
#   }
#   parameters{
#   real b0; //intercept
#   vector[3] b1; //SLOPES
#   real<lower=0> sigma;
#   }
#   transformed parameters{
#   vector[N] mu;
#   mu = b0 + b1[geno].*x; //DIFFERENT SLOPE FOR EACH GENOTYPE
#   }
#   model{
#   sigma ~ cauchy(0,5);
#   b0 ~ normal(0,10000);
#   b1 ~ normal(0,10000);
#   y ~ normal(mu, sigma);  //likelihood
#   }
#   generated quantities{
#   vector[N] e_y;
#   e_y = y - mu;
#   }
#   "
#   comp <- stan_model(model_code = Stan.model)
#   ##LOAD AND SPECIFY DATA
#   corr.model <- function(in_df=df,
#                          pdfname_normtest="normtest_corr_ugs_growth7.pdf",
#                          response='Dht_0926'){
#     setwd(wd)
#     response.var <- in_df[,c(response)]
#     model.components <- list( 'N' = nrow(in_df), 'y' = response.var, 'x' = in_df$avg_cond, 'geno' = as.numeric(in_df$factor_geno) ) 
#     iter <- 10000 
#     ##SAMPLING
#     fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3)
#     summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
#     rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ] ##GATHER RESIDUALS
#     dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
#     setwd(file.path(wd, 'Normality_tests'))
#     pdf(file=pdfname_normtest,height=7,width=7)
#     test_norm(rdf$mean); dev.off() ##NORMALITY PLOTS
#     return(summ_fit)
#   }
#   c1 <- corr.model(in_df=df,
#             pdfname_normtest="normtest_corr_ugs_growth7.pdf",
#             response='Dht_0926');c1
#   c2 <- corr.model(in_df=df,
#                    pdfname_normtest="normtest_corr_ugs_growth8.pdf",
#                    response='Dht_1002');c2
#   
###########################################################################################
## PHYSIOLOGICAL TRAITS OVER TIME (plots only)
###########################################################################################
  
  library(rstan)
  library(LambertW)
  options(mc.cores = parallel::detectCores())
  ##MODEL
  Stan.model <- "
  data{
  int<lower=0> N;
  int<lower=0> J;
  vector[N] y;
  matrix[N,J] X;
  }
  parameters{
  vector[J] b;
  real<lower=0> sigma;
  }
  transformed parameters{
  vector[N] mu;
  mu=X*b;
  }
  model{
  sigma ~ cauchy(0,5);
  b ~ normal(0,1000000);
  y ~ normal(mu, sigma);  //likelihood
  }
  generated quantities{
  vector[N] e_y;
  vector[15] Y;
  Y[1] = b[1];
  Y[2] = b[1] + b[2];
  Y[3] = b[1] + b[3];
  Y[4] = b[1] + b[4];
  Y[5] = b[1] + b[5];
  Y[6] = b[1] + b[6];
  Y[7] = b[1] + b[2] + b[6] + b[8];
  Y[8] = b[1] + b[3] + b[6] + b[9];
  Y[9] = b[1] + b[4] + b[6] + b[10];
  Y[10] = b[1] + b[5] + b[6] + b[11];
  Y[11] = b[1] + b[7];
  Y[12] = b[1] + b[2] + b[7] + b[12];
  Y[13] = b[1] + b[3] + b[7] + b[13];
  Y[14] = b[1] + b[4] + b[7] + b[14];
  Y[15] = b[1] + b[5] + b[7] + b[15];
  e_y = y - mu;
  }
  "
  comp <- stan_model(model_code = Stan.model) ##COMPILE
  ##LOAD DATA - CHANGE TO NAME FILES DIFFERENTLY IF DESIRED
  pdfname="FIGURE_S3a.pdf"
  df <- read.csv("phys_plants.csv",header=T)
  dm <- model.matrix(~as.factor(trt)*as.factor(geno), data=df)     ## MODEL MATRIX FOR STAN
  phys.time.model <- function(response.var=df$Anet5_0907,date_ID="week05"){
    model.components <- list( 'N' = nrow(df), 'y' = response.var, 'X' = dm, 'J' = ncol(dm) )
    iter <- 10000
    ##SAMPLE
    fit <- sampling( comp, data = model.components, iter = iter, warmup = iter/2, thin = 1, chains = 3, control=list(adapt_delta=0.99) )
    summ_fit <- summary(fit) ; summ_fit <- as.data.frame(summ_fit$summary)
    sdf <- summ_fit[grep("Y", rownames(summ_fit)), ]
    rdf <- summ_fit[grep("e_y", rownames(summ_fit)), ]
    setwd(wd)
    geno <- c(11,11,11,11,11,2,2,2,2,2,5,5,5,5,5)
    trt <- c(10,15,20,25,35,10,15,20,25,35,10,15,20,25,35)
    facet <- c(1,1,1,1,2,1,1,1,1,2,1,1,1,1,2) ## SATURATED TREATMENT ON SEPERATED X AXIS
    date <- c(rep(date_ID,15))
    sdf <- cbind(sdf,geno,trt,facet,date)
    return(sdf)
  }

  
  ##RUN MODELS, PLOT COMES LATER
  df1 <- phys.time.model()
  df2 <- phys.time.model(df$Anet6_0913,"week06")
  df3 <- phys.time.model(df$Anet7_0920,"week07")
  df4 <- phys.time.model(df$Anet8_0925,"week08")
  df5 <- phys.time.model(df$Anet9_1004,"week09")
  df6 <- phys.time.model(df$Anet10_1011,"week10")
  phys.time.data <- rbind(df1,df2,df3,df4,df5,df6)
  phys.time.data <- phys.time.data[,c(1,4,8,11,12,13,14)]
  response.label <- expression(A[net]~(mu~mol~m^{-2}~s^{-1}))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  pdf(file=pdfname,height=4,width=10)
  ##PLOT HEIGHT OVER TIME
  ggplot( phys.time.data, aes(x=as.numeric(date),y=mean,color=as.factor(trt))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("5","6","7","8","9","10")) +
    #facet_grid(.~facet, scales="free",space="free") +
    facet_grid(.~geno)
  dev.off()
  
  pdfname="FIGURE_S3b.pdf"
  df7 <- phys.time.model(df$gs5_0907,"week05")
  df8 <- phys.time.model(df$gs6_0913,"week06")
  df9 <- phys.time.model(df$gs7_0920,"week07")
  df10 <- phys.time.model(df$gs8_0925,"week08")
  df11 <- phys.time.model(df$gs9_1004,"week09")
  df12 <- phys.time.model(df$gs10_1011,"week10")
  phys.time.data <- rbind(df7,df8,df9,df10,df11,df12)
  phys.time.data <- phys.time.data[,c(1,4,8,11,12,13,14)]
  response.label <- expression(g[s]~(mol~H[2]*O~m^{-2}~s^{-1}))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  pdf(file=pdfname,height=4,width=10)
  ##PLOT HEIGHT OVER TIME
  ggplot( phys.time.data, aes(x=as.numeric(date),y=mean,color=as.factor(trt))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("5","6","7","8","9","10")) +
    #facet_grid(.~facet, scales="free",space="free") +
    facet_grid(.~geno)
  dev.off()
  
  pdfname="FIGURE_S3c.pdf"
  df13 <- phys.time.model(df$E5_0907,"week05")
  df14 <- phys.time.model(df$E6_0913,"week06")
  df15 <- phys.time.model(df$E7_0920,"week07")
  df16 <- phys.time.model(df$E8_0925,"week08")
  df17 <- phys.time.model(df$E9_1004,"week09")
  df18 <- phys.time.model(df$E10_1011,"week10")
  phys.time.data <- rbind(df13,df14,df15,df16,df17,df18)
  phys.time.data <- phys.time.data[,c(1,4,8,11,12,13,14)]
  response.label <- expression(E~(mmol~H[2]*O~m^{-2}~s^{-1}))
  COLOR.PAL <- c("#d7191c","#d95f02","#4daf4a","#1b9e77","#1f78b4") ##COLOR PALETTE FOR TREATMENTS OVER TIME
  pdf(file=pdfname,height=4,width=10)
  ##PLOT HEIGHT OVER TIME
  ggplot( phys.time.data, aes(x=as.numeric(date),y=mean,color=as.factor(trt))) +
    geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=0.5) +
    xlab("Week") +
    ylab(response.label) +
    geom_line() +
    geom_point(size=1) +
    theme_classic() +
    theme(legend.position="bottom") +
    scale_color_manual(values=COLOR.PAL, name="% VWC Treatment", labels=c("10","15","20","25","Sat'd")) +
    scale_x_continuous(breaks=c(1,2,3,4,5,6), labels=c("5","6","7","8","9","10")) +
    #facet_grid(.~facet, scales="free",space="free") +
    facet_grid(.~geno)
  dev.off()