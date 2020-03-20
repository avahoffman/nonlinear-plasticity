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


trait.model <-
  function(infile,
           response) {
    summ <-
      run_mcmc(comp = comp,
               infile = "data/biomass_plants.csv",
               response = "Bv")
    
    write.csv(gather_posterior_data(summ), 
              file = "output/posterior_output.csv",
              append = T)
    
    make_normality_plots(summ,
                         response = "Bv")
    
  }



## RUN MODEL FOR MORPHOLOGICAL TRAITS

p_1 <- 
  trait.model(
  infile = "data/biomass_plants.csv",
  response = "Bv"
)

p_2 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_Rl.pdf",
  pdfname_normtest = "normtest_Rl.pdf",
  response.label = expression(paste("Root length (mm)")),
  response = "Rl"
) + theme(legend.position = "none", axis.title.x = element_blank())
p_3 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_Rsa.pdf",
  pdfname_normtest = "normtest_Rsa.pdf",
  response.label = expression(paste("Root surface area (", mm ^
                                      2, ")")),
  response = "Rsa"
) + theme(legend.position = "none")
p_4 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_Rt.pdf",
  pdfname_normtest = "normtest_Rt.pdf",
  response.label = expression(paste("Root tips")),
  response = "Rt"
) + theme(legend.position = "none")
p_5 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_Rv.pdf",
  pdfname_normtest = "normtest_Rv.pdf",
  response.label = expression(paste("Root volume (", mm ^
                                      3, ")")),
  response = "Rv"
) + theme(legend.position = "none")
p_6 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_LA.pdf",
  pdfname_normtest = "normtest_LA.pdf",
  response.label = expression(paste("Leaf area (", mm ^
                                      2, ")")),
  response = "LA"
) + theme(legend.position = "none")
p_7 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_rsa_la.pdf",
  pdfname_normtest = "normtest_rsa_la.pdf",
  response.label = expression(paste(
    "Root to leaf area ratio (", mm ^ 2, " ", mm ^ -2, ")"
  )),
  response = "Rsa_la"
) + theme(legend.position = "none")
p_8 <- trait.model(
  infile = "biomass_plants.csv",
  pdfname = "rxn_norm_srl.pdf",
  pdfname_normtest = "normtest_srl.pdf",
  response.label = expression(paste("Specific root length (", mm, " ", g ^
                                      -1, ")")),
  response = "srl"
) + theme(legend.position = "none")