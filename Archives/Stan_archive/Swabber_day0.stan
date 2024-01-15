data{
    int<lower=0> Ntot;                           // Number of PCR data points (all swabs)
    int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD (non-censored swabs)
    int<lower=0> n_id;                           // Number of individuals
    int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
    real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
    real log10_cens_vl[Ntot];                    // censoring value for censored observation (vl at CT = 40)
    real<lower=0> Time_max;                      // Max follow-up time
    real<lower=0,upper=Time_max> obs_day[Ntot];  // Time since randomisation for sample
    real RNaseP[Ntot];                           // Scaled RNaseP CT values (mean 0)

    int<lower=0> n_swabber;
    int<lower=0> swabber_id[Ntot];
    
    int<lower=0> swab_ID[Ntot];

    // priors
   // real alpha_0_prior_mean; // prior mean intercept alpha_0
  //  real alpha_0_prior_sd;   // prior sd intercept alpha_0
//    real beta_0_prior_mean;  // prior mean slope beta_0
  //  real beta_0_prior_sd;    // prior sd slope beta_0
   // real sigma_logvl_mean;   // prior mean of sd in error model sigma_logvl
   // real sigma_logvl_sd;     // prior sd of sd in error model sigma_logvl
}


parameters{
  // hyperparameters
  
  // Population parameters
  real mu_0;
  real gamma_rnasep;                        // Adjustment for RNaseP
  real beta_swabID;

  // Random effects
  real theta_rand[n_id];            // individual random effects vector
  // Measurement error
  //real rand_sigma_logvl[n_swabber];
  real<lower=0> sigma_logvl[n_swabber];
  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}


transformed parameters{
    real pred_log10_vl[Ntot];
  //  real sigma_logvl[n_swabber];
  
  // calculate predicted log viral load under the model parameters
    for(i in 1:Ntot){
      pred_log10_vl[i] = mu_0 + theta_rand[id[i]] + gamma_rnasep*RNaseP[i] + swab_ID[i]*beta_swabID;
    }
    
  //  for(i in 1:n_swabber)
  //  sigma_logvl[i] = sigma0_logvl * exp(rand_sigma_logvl[i]);

}

model{
  //***** Likelihood *****
  // Non censored observations
  for(i in 1:N_obs){
      target += student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl[swabber_id[i]]);
  }
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl[swabber_id[i]]);
  }

  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ exponential(1);
 // rand_sigma_logvl ~ normal(0,1);
  // Population parameters
  mu_0 ~ normal(0,5);
  
  theta_rand ~ normal(0,1);
}
