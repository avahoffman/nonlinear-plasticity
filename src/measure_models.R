###########################################################################################
## MODEL FOR MORPHOLOGY, PHYSIOLOGICAL TRAITS OF INTEREST
###########################################################################################
library(rstan)
library(LambertW)
options(mc.cores = parallel::detectCores())

# Declare model
Stan_model <- "
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


run_measure_model <-
  function(comp,
           infile = "data/biomass_plants.csv",
           response) {
    # This function...
    
    # Summarize model
    summ <-
      run_mcmc(comp = comp,
               infile = infile,
               response = response)
    
    # Append parameter results
    if (response == "Bv"){
      write_col = NA
    } else {
      write_col = F
    }
    write.table(gather_posterior_data(summ_fit = summ, response = response),
              file = "output/posterior_output.csv",
              sep = ",",
              col.names = write_col,
              row.names = T,
              append = T)
    
    # Run posterior predictive checks
    make_normality_plots(summ,
                         response = response)
  }


do_measure_mcmc_sampling <-
  function() {
    comp <- 
      stan_model(model_code = Stan_model)
    
    # Morphology
    run_measure_model(comp, response = "Bv")
    run_measure_model(comp, response = "Rl")
    run_measure_model(comp, response = "Rsa")
    run_measure_model(comp, response = "Rt")
    run_measure_model(comp, response = "Rv")
    run_measure_model(comp, response = "LA")
    run_measure_model(comp, response = "Rsa_la")
    run_measure_model(comp, response = "srl")
    
    # Physiology
    run_measure_model(comp, 
                      infile = "data/phys_plants.csv",
                      response = "uAnet")
    run_measure_model(comp,
                      infile = "data/phys_plants.csv",
                      response = "ugs")
    run_measure_model(comp,
                      infile = "data/phys_plants.csv",
                      response = "uE")
    
    # 
  }