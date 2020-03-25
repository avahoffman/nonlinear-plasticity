###########################################################################################
## MODEL FOR MORPHOLOGY, PHYSIOLOGICAL TRAITS OF INTEREST
###########################################################################################
library(rstan)
library(LambertW)
options(mc.cores = parallel::detectCores())


# Declare model -
Stan_model_normal <- "
        data{
            int<lower=0> N; // Number of observations
            int<lower=0> J; // Parameters; or columns of model matrix
            vector[N] y; // Response variable
            matrix[N,J] X; // Model matrix (G x T)
        }
        parameters{
            vector[J] b; // Regression parameters
            real<lower=0> sigma; // Variance parameter
        }
        transformed parameters{
            vector[N] mu; // Expected, modelled values
            mu = X * b;
        }
        model{
            b[1] ~ cauchy(0,100); // Prior for the intercept following Gelman 2008
  
            for(i in 2:J)
               b[i] ~ cauchy(0,2.5); // Prior for the slopes following Gelman 2008
            
            y ~ normal(mu, sigma);  // likelihood
        }
        generated quantities{
            vector[N] e_y;
            vector[15] Y;
            vector[3] G;
            vector[5] T;
            vector[1] T_E;
            vector[1] G_E;
            vector[1] I_E;
    
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
    
            G[1] = ( b[1] + (b[1] + b[2]) + (b[1] + b[3]) + (b[1] + b[4]) + (b[1] + b[5])) / 5 ;
            G[2] = ( (b[1] + b[6]) + (b[1] + b[2] + b[6] + b[8]) + (b[1] + b[3] + b[6] + b[9]) + (b[1] + b[4] + b[6] + b[10]) + (b[1] + b[5] + b[6] + b[11])) / 5 ;
            G[3] = ( (b[1] + b[7]) + (b[1] + b[2] + b[7] + b[12]) + (b[1] + b[3] + b[7] + b[13]) + (b[1] + b[4] + b[7] + b[14]) + (b[1] + b[5] + b[7] + b[15])) / 5 ;
    
            T[1] = ( b[1] + (b[1] + b[6]) + (b[1] + b[7])) / 3 ;
            T[2] = ( (b[1] + b[2]) + (b[1] + b[2] + b[6] + b[8]) + (b[1] + b[2] + b[7] + b[12]) ) / 3 ;
            T[3] = ( (b[1] + b[3]) + (b[1] + b[3] + b[6] + b[9]) + (b[1] + b[3] + b[7] + b[13]) ) / 3 ;
            T[4] = ( (b[1] + b[4]) + (b[1] + b[4] + b[6] + b[10]) + (b[1] + b[4] + b[7] + b[14]) ) / 3 ;
            T[5] = ( (b[1] + b[5]) + (b[1] + b[5] + b[6] + b[11]) + (b[1] + b[5] + b[7] + b[15]) ) / 3 ;
    
            G_E[1] = ( b[6] + b[7] ) / 2 ;
    
            T_E[1] = ( b[2] + b[3] + b[4] + b[5] ) / 4 ;
    
            I_E[1] = ( b[8] + b[9] + b[10] + b[11] + b[12] + b[13] + b[14] + b[15] ) / 8 ;
    
            e_y = y - mu;
        }
      "

Stan_model_gamma <- "
        data{
            int<lower=0> N; // Number of observations
            int<lower=0> J; // Parameters; or columns of model matrix
            vector<lower=0>[N] y; // Response variable
            matrix[N,J] X; // Model matrix (G x T)
        }
        parameters{
            vector[J] b; // Regression parameters
            real<lower=0> phi; // Gamma variance parameter
        }
        transformed parameters{
            vector[N] mu; // Expected, modelled values
            vector[N] alpha; // Shape for gamma distribution
            vector[N] beta; // Rate for gamma distribution
            
            mu = exp( X * b); 
            alpha = ( mu .* mu ) / phi;
            beta = mu / phi;
        }
        model{
            b[1] ~ cauchy(0,10); // Prior for the intercept following Gelman 2008
  
            for(i in 2:J)
               b[i] ~ cauchy(0,2.5); // Prior for the slopes following Gelman 2008
               
            y ~ gamma(alpha, beta); // Likelihood
        }
        generated quantities{
            vector[N] e_y; // Draws of the difference between y and mu
            vector[15] Y; // Modeled posteriors of the G x T groups
            vector[3] G; // Genotype distribution averaged over treatments
            vector[5] T; // Treatment distribution averaged over genotypes
            vector[1] T_E; // Average treatment effect
            vector[1] G_E; // Average genotype effect
            vector[1] I_E; // Average interaction effect
        
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
            
            G[1] = exp( b[1] + (b[1] + b[2]) + (b[1] + b[3]) + (b[1] + b[4]) + (b[1] + b[5])) / 5 ;
            G[2] = exp( (b[1] + b[6]) + (b[1] + b[2] + b[6] + b[8]) + (b[1] + b[3] + b[6] + b[9]) + (b[1] + b[4] + b[6] + b[10]) + (b[1] + b[5] + b[6] + b[11])) / 5 ;
            G[3] = exp( (b[1] + b[7]) + (b[1] + b[2] + b[7] + b[12]) + (b[1] + b[3] + b[7] + b[13]) + (b[1] + b[4] + b[7] + b[14]) + (b[1] + b[5] + b[7] + b[15])) / 5 ;
        
            T[1] = exp( b[1] + (b[1] + b[6]) + (b[1] + b[7])) / 3 ;
            T[2] = exp( (b[1] + b[2]) + (b[1] + b[2] + b[6] + b[8]) + (b[1] + b[2] + b[7] + b[12]) ) / 3 ;
            T[3] = exp( (b[1] + b[3]) + (b[1] + b[3] + b[6] + b[9]) + (b[1] + b[3] + b[7] + b[13]) ) / 3 ;
            T[4] = exp( (b[1] + b[4]) + (b[1] + b[4] + b[6] + b[10]) + (b[1] + b[4] + b[7] + b[14]) ) / 3 ;
            T[5] = exp( (b[1] + b[5]) + (b[1] + b[5] + b[6] + b[11]) + (b[1] + b[5] + b[7] + b[15]) ) / 3 ;
            
            G_E[1] = ( b[6] + b[7] ) / 2 ; 
            
            T_E[1] = ( b[2] + b[3] + b[4] + b[5] ) / 4 ; 
            
            I_E[1] = ( b[8] + b[9] + b[10] + b[11] + b[12] + b[13] + b[14] + b[15] ) / 8 ;
            
            e_y = y - mu;
        }
      "


run_measure_model <-
  function(comp,
           infile = "data/biomass_plants_clean.csv",
           response,
           recovery = F) {
    # This function...
    
    # Summarize model
    summ <-
      run_mcmc(comp = comp,
               infile = infile,
               response = response)
    
    # Append parameter results
    if (response == "Bv") {
      write_col = NA
    } else {
      write_col = F
    }
    write.table(
      gather_posterior_data(
        summ_fit = summ,
        response = response,
        recovery = recovery
      ),
      file = "output/posterior_output.csv",
      sep = ",",
      col.names = write_col,
      row.names = T,
      append = T
    )
    
    # Run posterior predictive checks
    make_normality_plots(summ,
                         response = response)
  }


do_measure_mcmc_sampling <-
  function() {
    comp_gamma <-
      stan_model(model_code = Stan_model_gamma, verbose = T)
    comp_normal <-
      stan_model(model_code = Stan_model_normal, verbose = T)
    
    # Morphology
    # All of these values are positive
    run_measure_model(comp_gamma, response = "Bv")
    run_measure_model(comp_gamma, response = "Br")
    run_measure_model(comp_gamma, response = "DMCv")
    run_measure_model(comp_gamma, response = "DMCr")
    run_measure_model(comp_gamma, response = "A_B")
    run_measure_model(comp_gamma, response = "SLA")
    run_measure_model(comp_gamma, response = "LA")
    run_measure_model(comp_gamma, response = "Rl")
    run_measure_model(comp_gamma, response = "Rsa")
    run_measure_model(comp_gamma, response = "Rd")
    run_measure_model(comp_gamma, response = "Rv")
    run_measure_model(comp_gamma, response = "Rt")
    run_measure_model(comp_gamma, response = "Rsa_la")
    run_measure_model(comp_gamma, response = "srl")
    run_measure_model(comp_gamma, response = "ssa")
    
    # Physiology
    # All of these values are positive
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "uAnet")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "ugs")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "ufv")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "uWUEi")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "max_Anet")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "max_gs")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "max_fv")
    run_measure_model(comp_gamma,
                      infile = "data/phys_plants_clean.csv",
                      response = "max_WUEi")
    
    # All plants
    # RGR is sometimes negative, use normal model
    run_measure_model(comp_gamma,
                      infile = "data/all_plants_clean.csv",
                      response = "uH")
    run_measure_model(comp_normal,
                      infile = "data/all_plants_clean.csv",
                      response = "urgr")
    run_measure_model(comp_gamma,
                      infile = "data/all_plants_clean.csv",
                      response = "max_H")
    run_measure_model(comp_normal,
                      infile = "data/all_plants_clean.csv",
                      response = "max_rgr")
    
    # Recovery
    # RGR is sometimes negative, use normal model
    # run_measure_model(comp_gamma,
    #                   infile = "data/recovery_plants_clean.csv",
    #                   response = "Bv",
    #                   recovery = T)
    # run_measure_model(comp_gamma,
    #                   infile = "data/recovery_plants_clean.csv",
    #                   response = "Bf",
    #                   recovery = T)
    # run_measure_model(comp_gamma,
    #                   infile = "data/recovery_plants_clean.csv",
    #                   response = "Br",
    #                   recovery = T)
    # run_measure_model(comp_gamma,
    #                   infile = "data/recovery_plants_clean.csv",
    #                   response = "Brh",
    #                   recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "B_above",
                      recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "B_below",
                      recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "B_total",
                      recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "A_B",
                      recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "uH",
                      recovery = T)
    run_measure_model(comp_normal,
                      infile = "data/recovery_plants_clean.csv",
                      response = "urgr",
                      recovery = T)
    run_measure_model(comp_gamma,
                      infile = "data/recovery_plants_clean.csv",
                      response = "max_H",
                      recovery = T)
    run_measure_model(comp_normal,
                      infile = "data/recovery_plants_clean.csv",
                      response = "max_rgr",
                      recovery = T)
  }