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
  real study_time[n_id];
  
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
  int num_knots_alpha;            // num of knots
  vector[num_knots_alpha] knots_alpha;  // the sequence of knots
  int spline_degree_alpha;        // the degree of spline (is equal to order - 1)
  
  int num_knots_beta;            // num of knots
  vector[num_knots_beta] knots_beta;  // the sequence of knots
  int spline_degree_beta;        // the degree of spline (is equal to order - 1)
}

transformed data {
  // For individual random effects
  vector[2] zeros2;
  
  // B-splines for Intercepts
  int num_basis_alpha = num_knots_alpha + spline_degree_alpha - 1; // total number of B-splines
  matrix[num_basis_alpha, n_id] B_alpha;  // matrix of B-splines
  vector[spline_degree_alpha + num_knots_alpha] ext_knots_temp_alpha;
  vector[2*spline_degree_alpha + num_knots_alpha] ext_knots_alpha; // set of extended knots
  
  // B-splines for Slope
  int num_basis_beta = num_knots_beta + spline_degree_beta - 1; // total number of B-splines
  matrix[num_basis_beta, n_id] B_beta;  // matrix of B-splines
  vector[spline_degree_beta + num_knots_beta] ext_knots_temp_beta;
  vector[2*spline_degree_beta + num_knots_beta] ext_knots_beta; // set of extended knots

  // For individual random effects
  for(i in 1:2) zeros2[i] = 0;
  
  // B-splines for Intercepts
  ext_knots_temp_alpha = append_row(rep_vector(knots_alpha[1], spline_degree_alpha), knots_alpha);
  ext_knots_alpha = append_row(ext_knots_temp_alpha, rep_vector(knots_alpha[num_knots_alpha], spline_degree_alpha));
  for (ind in 1:num_basis_alpha){
    B_alpha[ind,:] = to_row_vector(build_b_spline(study_time, to_array_1d(ext_knots_alpha), ind, spline_degree_alpha + 1));
  }
  B_alpha[num_knots_alpha + spline_degree_alpha - 1, n_id] = 1; 
  
  // B-splines for Slope
  ext_knots_temp_beta = append_row(rep_vector(knots_beta[1], spline_degree_beta), knots_beta);
  ext_knots_beta = append_row(ext_knots_temp_beta, rep_vector(knots_beta[num_knots_beta], spline_degree_beta));
  for (ind in 1:num_basis_beta){
    B_beta[ind,:] = to_row_vector(build_b_spline(study_time, to_array_1d(ext_knots_beta), ind, spline_degree_beta + 1));
  }
  B_beta[num_knots_beta + spline_degree_beta - 1, n_id] = 1; 
}

parameters {
  // hyperparameters
  cholesky_factor_corr[2] L_Omega;     // correlation matrix
  vector<lower=0>[2] sigmasq_u;        // variance of random effects
  
  // Measurement error
  real<lower=0> sigma_logvl;
  
  // Population parameters
 // real alpha_0;                             // population intercept
  real beta_0;                              // population slope
  real gamma_rnasep;                        // Adjustment for RNaseP
  vector[K_trt-1] trt_effect;               // Estimates of the treatment effect
  vector[K_cov_slope] slope_coefs;          // Covariate effects on slope
  vector[K_cov_intercept] intercept_coefs;  // Covariate effects on intercept
  
  // Random effects
  vector[2] theta_rand_id[n_id];            // individual random effects vector
  
  // Degrees of freedom for the t-distribution error model
  real<lower=0> t_dof;
  
  // Splines for intercepts
  row_vector[num_basis_alpha] a_raw_alpha;
  real<lower=0> tau_alpha;
  real a0_alpha; // intercept at time 0
  
  // Splines for slope
  row_vector[num_basis_beta] a_raw_beta;
  real<lower=0> tau_beta;
  real a0_beta; // intercept at time 0
}

transformed parameters {
  real pred_log10_vl[Ntot];
  vector[Ntot] trt_slope;
  vector[Ntot] beta_cov;
  vector[Ntot] alpha_cov;

  row_vector[num_basis_alpha] a_alpha;
  vector[n_id] alpha_hat;
  row_vector[num_basis_beta] a_beta;
  vector[n_id] beta_hat;
  
  // Splines for intercept
  a_alpha[1] = a_raw_alpha[1];
  for (i in 2:num_basis_alpha){
      a_alpha[i] = a_alpha[i-1] + a_raw_alpha[i]*tau_alpha;
  }
  alpha_hat = a0_alpha*to_vector(study_time) + to_vector(a_alpha*B_alpha);
  
  // Splines for slope 
  a_beta[1] = a_raw_beta[1];
  for (i in 2:num_basis_beta){
      a_beta[i] = a_beta[i-1] + a_raw_beta[i]*tau_beta;
  }
  beta_hat = a0_beta*to_vector(study_time) + to_vector(a_beta*B_beta);
  
  
  
  {
    vector[K_trt] trt_effect_prime;
    
    trt_effect_prime = append_row(0, trt_effect);
    trt_slope = trt_mat * trt_effect_prime;
    
    // make individual covariate transform
    beta_cov = x_slope*slope_coefs;
    alpha_cov = x_intercept*intercept_coefs;
    
    // calculate predicted log viral load under the model parameters
    for(i in 1:Ntot){
      pred_log10_vl[i] =
      alpha_hat[id[i]] + theta_rand_id[id[i]][1] + alpha_cov[i] +
      gamma_rnasep*RNaseP[i] +
      beta_hat[id[i]]*exp(trt_slope[i]+theta_rand_id[id[i]][2]+beta_cov[i])*obs_day[i];
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
  //alpha_0 ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  //beta_0 ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  gamma_rnasep ~ normal(0,1);
  slope_coefs ~ normal(0,slope_coefs_sd);
  intercept_coefs ~ normal(0,intercept_coefs_sd);
  trt_effect ~ normal(0,trt_effect_sd); // Treatment effect
  
  
  // Splines for intercepts
  a_raw_alpha ~ normal(0, 1);
  a0_alpha ~ normal(alpha_0_prior_mean,alpha_0_prior_sd);
  tau_alpha ~ normal(0, 1);
  
   // Splines for intercepts
  a_raw_beta ~ normal(0, 1);
  a0_beta ~ normal(beta_0_prior_mean,beta_0_prior_sd);
  tau_beta ~ normal(0, 1);
  
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
    slope[i] = beta_hat[i]*exp(trt_slope[j]+theta_rand_id[i][2]+beta_cov[j]);
  }
}
