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

    // priors
    real alpha_0_prior_mean; // prior mean intercept alpha_0
    real alpha_0_prior_sd;   // prior sd intercept alpha_0
    real beta_0_prior_mean;  // prior mean slope beta_0
    real beta_0_prior_sd;    // prior sd slope beta_0
    real sigma_logvl_mean;   // prior mean of sd in error model sigma_logvl
    real sigma_logvl_sd;     // prior sd of sd in error model sigma_logvl
}

transformed data {
  vector[2] zeros2;
  for(i in 1:2) zeros2[i] = 0;
}

parameters{
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects
  
  // Population parameters
  real alpha_0;                             // population intercept
  real beta_0;                              // population slope
  real gamma_rnasep;                        // Adjustment for RNaseP

  // Random effects
  vector[2] theta_rand_id[n_id];            // individual random effects vector
  // Measurement error
  real<lower=0> sigma_logvl;
  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}


transformed parameters{
    real pred_log10_vl[Ntot];
    real alpha[Ntot];
    real beta[Ntot];
  
  // calculate predicted log viral load under the model parameters
    for(i in 1:Ntot){
      alpha[i] =  alpha_0 + theta_rand_id[id[i]][1] + gamma_rnasep*RNaseP[i];
      beta[i] =  beta_0 * exp(theta_rand_id[id[i]][2]);
      pred_log10_vl[i] = alpha[i] + beta[i]*obs_day[i];
    }
}

model{
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_logvl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }

  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ normal(sigma_logvl_mean,sigma_logvl_sd) T[0,];
  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2); // covariance matrix - random effects for individs
  // individual random effects
  for(i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));
  
}


generated quantities {
  real preds[Ntot]; // For plotting
  vector[Ntot] log_lik;

  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
}
