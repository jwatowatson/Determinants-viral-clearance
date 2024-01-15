//   Copyright 2022 James Watson, Mahidol Oxford Tropical Medicine Research Unit
//
//   Licensed under a CC-BY License
/**
Model 1 with the following characteristics:
- analysis on the copies per ml scale (batch effect adjustement done without uncertainty propagation)
- Adjustment for RNaseP
- Covariate adjustment (optional)

**/

// taken from https://github.com/milkha/Splines_in_Stan/blob/master/splines_in_stan.Rmd
functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]); 
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) / 
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) / 
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines 
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) + 
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}

data {
  // Patient data
  int<lower=0> Ntot;                           // Number of PCR data points
  int<lower=0,upper=Ntot> N_obs;               // Number of PCR data points with viral load less than LOD
  int<lower=0> n_id;                           // Number of individuals
  int<lower=1,upper=n_id> id[Ntot];            // Patient identifier for each PCR sample
  int<lower=1,upper=Ntot> ind_start[n_id];     // Starting index for each patient: used in generated quantities block
  real<lower=0> Time_max;                      // Max follow-up time
  real<lower=0,upper=Time_max> obs_day[Ntot];  // Time since randomisation for sample
  real log_10_vl[Ntot];                        // log base 10 viral load in copies per mL
  real log10_cens_vl[Ntot];                    // censoring value for censored observation
  real RNaseP[Ntot];                           // Scaled RNaseP CT values (mean 0)
  int<lower=1> K_trt;                          // Number of treatment arms
  matrix[Ntot,K_trt] trt_mat;                  // Trt matrix
  int<lower=0> K_cov_intercept;                // number of columns in covariate design matrix for intercept
  matrix[Ntot,K_cov_intercept] x_intercept;    // covariate design matrix for regression onto intercept
  int<lower=0> K_cov_slope;                    // number of columns in covariate design matrix for slope
  matrix[Ntot,K_cov_slope] x_slope;            // covariate design matrix for regression onto slope
  
  // priors
  real alpha_0_prior_mean; // prior mean intercept
  real alpha_0_prior_sd;   // prior sd intercept
  
  real beta_0_prior_mean;  // prior mean slope
  real beta_0_prior_sd;    // prior sd slope
  
  real trt_effect_sd;   // prior sd on treatment effect
  
  real sigma_logvl_mean;   // prior mean of sd in error model
  real sigma_logvl_sd;     // prior sd of sd in error model
  
  real slope_coefs_sd;     // prior sd for prior on covariate effect on slope
  real intercept_coefs_sd; // prior sd for prior on covariate effect on intercept
  
  // temporal spline
  int num_knots;            // num of knots
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
}

transformed data {
  vector[2] zeros2;
  
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, n_id] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));
  B[num_knots + spline_degree - 1, num_data] = 1; 
  
  for(i in 1:2) zeros2[i] = 0;

}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects
  
  // Measurement error
  real<lower=0> sigma_logvl;
  
  // Population parameters
  real alpha_0;                             // population intercept
  real beta_0;                              // population slope
  real gamma_rnasep;                        // Adjustment for RNaseP
  vector[K_trt-1] trt_effect;               // Estimates of the treatment effect
  vector[K_cov_slope] slope_coefs;          // Covariate effects on slope
  vector[K_cov_intercept] intercept_coefs;  // Covariate effects on intercept
  
  // Random effects
  vector[2] theta_rand_id[n_id];            // individual random effects vector
  
  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[Ntot] trt_slope;
  vector[Ntot] beta_cov;
  vector[Ntot] alpha_cov;
  {
    vector[K_trt] trt_effect_prime;
    
    
    trt_effect_prime = append_row(0, trt_effect);
    trt_slope = trt_mat * trt_effect_prime;
    
    // make individual covariate transform
    beta_cov = x_slope*slope_coefs;
    beta_time_spline = 
    alpha_cov = x_intercept*intercept_coefs;
    
    // calculate predicted log viral load under the model parameters
    for(i in 1:Ntot){
      pred_log10_vl[i] =
      alpha_0 + theta_rand_id[id[i]][1] + alpha_cov[i] + // intercept for patient id[i]
      gamma_rnasep*RNaseP[i] + // RNaseP adjustment
      beta_0*exp(trt_slope[i]+theta_rand_id[id[i]][2]+beta_cov[i])*obs_day[i]; // slope
    }
  }
}

model {
  //***** Prior *****
  // error model - degrees of freedom
  t_dof ~ exponential(1);
  // error model variance
  sigma_logvl ~ normal(sigma_logvl_mean,sigma_logvl_sd) T[0,];
  
  // random effects
  sigmasq_u[1] ~ exponential(1);
  sigmasq_u[2] ~ exponential(1);
  L_Omega ~ lkj_corr_cholesky(2); // covariance matrix - random effects for individs
  // individual random effects
  for(i in 1:n_id) theta_rand_id[i] ~ multi_normal_cholesky(zeros2, diag_pre_multiply(sigmasq_u, L_Omega));
  
  // Population parameters
  alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);
  slope_coefs ~ normal(0,slope_coefs_sd);
  intercept_coefs ~ normal(0,intercept_coefs_sd);
  trt_effect ~ normal(0,trt_effect_sd); // Treatment effect
  
  //***** Likelihood *****
  // Non censored observations
  log_10_vl[1:N_obs] ~ student_t(t_dof, pred_log10_vl[1:N_obs], sigma_logvl);
  
  // Censored observations
  for(i in (N_obs+1):Ntot){
    target += student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  
}

generated quantities {
  real preds[Ntot]; // For plotting
  vector[Ntot] log_lik;
  vector[n_id] slope;
  
  for(i in 1:N_obs){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lpdf(log_10_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in (N_obs+1):Ntot){
    preds[i] = pred_log10_vl[i]-gamma_rnasep*RNaseP[i];
    log_lik[i] = student_t_lcdf(log10_cens_vl[i] | t_dof, pred_log10_vl[i], sigma_logvl);
  }
  for(i in 1:n_id){
    int j = ind_start[i];
    slope[i] = beta_0*exp(trt_slope[j]+theta_rand_id[id[j]][2]+beta_cov[j]);
  }
}
