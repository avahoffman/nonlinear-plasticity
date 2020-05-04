###########################################################################################
## MODEL FOR THETA PARAMETER : PROBABILITY OF FLOWERING OR PROBABILITY OF RE-SPROUTING
###########################################################################################
library(rstan)
library(LambertW)
options(mc.cores = parallel::detectCores())


Stan_model_theta <- "
        data { 
          int<lower=0> N; // Number of observations
          int<lower=0> G; // number of genotypes
          int<lower=0,upper=1> w[N]; // Did / didn't flower data OR did / didn't resprout data
          int<lower=0,upper=G> geno[N]; // Genotype index
        } 
        parameters {
          vector[G] theta; // Probability of flowering / resprouting
        } 
        model {
          theta ~ beta(1,1);
          
          for (n in 1:N) 
            w[n] ~ bernoulli(theta[geno[n]]);
        }
      "


run_flwr_rh_model <-
  function(comp,
           response) {
    # This function wraps the sampling
    
    # Summarize model
    summ <-
      run_flwr_rh_mcmc(comp,
                       response = response)
    
    gather_flwr_rh_posterior_data(summ, response = response)
    
    # Append parameter results
    if (response == "did_flwr") {
      write_col = NA
    } else {
      write_col = F
    }
    write.table(
      gather_flwr_rh_posterior_data(summ, 
                                    response = response),
      file = "output/posterior_output_theta.csv",
      sep = ",",
      col.names = write_col,
      row.names = T,
      append = T
    )
    
  }


do_flwr_rh_mcmc_sampling <-
  function() {
    # This function runs the posterior sampling of flowering and re-sprouting
    # probability
    
    # Compile Stan model
    comp_theta <-
      stan_model(model_code = Stan_model_theta, verbose = T)
    
    # Flowering probability
    run_flwr_rh_model(comp_theta, response = "did_flwr")
    
    # Resprouting probability
    run_flwr_rh_model(comp_theta, response = "did_rh")
    
  }