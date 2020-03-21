###########################################################################################
## MODEL FOR RECOVERY - Flowering first, different distribution
## NOTE: run theta parameter seperately for genotype-level only flowering
###########################################################################################
library(rstan)
library(LambertW)
options(mc.cores = parallel::detectCores())

# Declare model
Stan_model <- "
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


run_flwr_model <-
  function(comp) {
    # This function...
    
    # Summarize model
    summ <-
      run_flwr_mcmc(comp)
    
    
  }


## OUTPUT DATA
dir.create(file.path(wd, 'Parameter_estimates'), showWarnings = FALSE)
setwd(file.path(wd, 'Parameter_estimates'))
ests <- (summ_fit[grep("b", rownames(summ_fit)), c(4, 8)])
write.csv(ests, file = paste(response, "_parameters_", infile, sep = ""))
setwd(wd)
thetas <- summ_fit[grep("theta", rownames(summ_fit)),]
sdf <- summ_fit[grep("Y", rownames(summ_fit)),]
rdf <- summ_fit[grep("e_y", rownames(summ_fit)),]
dir.create(file.path(wd, 'Normality_tests'), showWarnings = FALSE)
setwd(file.path(wd, 'Normality_tests'))
pdf(file = pdfname_normtest,
    height = 7,
    width = 7)
test_norm(rdf$mean)
dev.off()
setwd(wd) ##PLOT NORMALITY PLOTS
##PLOT REACTION NORMS
geno <- c(11, 11, 11, 11, 11, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5)
trt <- c(10, 15, 20, 25, 35, 10, 15, 20, 25, 35, 10, 15, 20, 25, 35)
facet <- c(1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2)
plot.params <- cbind(sdf, geno, trt, facet)
gg1 <-
  ggplot(plot.params, aes(x = trt, y = mean, color = as.factor(geno))) +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 1) +
  xlab("% VWC") +
  ylab(response.label) +
  geom_line() +
  geom_point(size = 1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_viridis(discrete = T,
                      end = 0.9,
                      name = "Genotype") +
  scale_x_continuous(
    breaks = c(10, 15, 20, 25, 35),
    labels = c("10", "15", "20", "25", "Sat'd")
  ) +
  facet_grid(. ~ facet, scales = "free", space = "free") +
  theme(strip.text = element_blank())
gg1
dir.create(file.path(wd, 'Reaction_Norms'), showWarnings = FALSE)
setwd(file.path(wd, 'Reaction_Norms'))
ggsave(
  filename = pdfname,
  plot = gg1,
  height = 4,
  width = 4
)
setwd(wd)