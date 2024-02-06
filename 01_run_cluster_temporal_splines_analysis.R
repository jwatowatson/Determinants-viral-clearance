# script for running all the stan models with all settings on the BMRC cluster
args = commandArgs(trailingOnly = FALSE) # comes from the SGE_TASKID in *.sh file
job_i = as.numeric(args[6])
print(paste0("job(i) = ", job_i)) # this will print out in the *.o file

## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)
library(dplyr)

source('functions.R')
source('priors.R')

load('Rout/model_run_setup_Unblinded_meta.RData')

Max_job = nrow(model_settings)


  if(job_i > Max_job) stop('no model setting corresponding to job ID')
  
  writeLines('Doing the following job:')
  print(model_settings[job_i, ])
  
  Dmax <- model_settings$Dmax[job_i]
  
  # Analysis data
  platcov_dat_analysis <- platcov_dat_analysis_list[[model_settings$dataset[job_i]]]
  stan_input_job <- stan_inputs[[model_settings$dataset[job_i]]]
  
  
  options(mc.cores = model_settings$Nchain[job_i])
  stopifnot(model_settings$Nchain[job_i]>getDoParWorkers()) # check worker number assigned
  
  mod = stan_model(file = as.character(model_settings$mod[job_i])) # compile 
  
  analysis_data_stan = stan_input_job$analysis_data_stan
  analysis_data_stan$trt_mat = stan_input_job$Trt_matrix
  analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)
  
  x_intercept = stan_input_job$cov_matrices$X_int[[model_settings$cov_matrices[job_i]]]
  if(ncol(x_intercept)==0) x_intercept = array(0, dim=c(nrow(x_intercept),1))
  analysis_data_stan$x_intercept = x_intercept
  analysis_data_stan$K_cov_intercept= ncol(x_intercept)
  
  x_slope = stan_input_job$cov_matrices$X_slope[[model_settings$cov_matrices[job_i]]]
  if(ncol(x_slope)==0) x_slope = array(0, dim=c(nrow(x_slope),1))
  analysis_data_stan$x_slope = x_slope
  analysis_data_stan$K_cov_slope=ncol(x_slope)
  
  num_knots_alpha <- model_settings$num_knots_alpha[job_i]
  num_knots_beta <- model_settings$num_knots_beta[job_i]
  
  study_time <- as.vector(unique(platcov_dat_analysis[,c("ID", "Study_time")])[,2])
  study_time <- study_time[[1]][,1]
  
  knots_alpha <- unname(quantile(study_time, probs=seq(from=0, to=1, length.out = num_knots_alpha)))
  knots_beta <- unname(quantile(study_time, probs=seq(from=0, to=1, length.out = num_knots_beta)))
  
  analysis_data_stan$num_knots_alpha <- num_knots_alpha
  analysis_data_stan$knots_alpha <- knots_alpha
  analysis_data_stan$spline_degree_alpha <- model_settings$spline_degree_alpha[job_i]
  
  analysis_data_stan$num_knots_beta <- num_knots_beta
  analysis_data_stan$knots_beta <- knots_beta
  analysis_data_stan$spline_degree_beta <- model_settings$spline_degree_beta[job_i]
  
  analysis_data_stan$study_time <- study_time
  
  # sample posterior
  out = sampling(mod, 
                 data=c(analysis_data_stan,
                        all_priors[[model_settings$prior[job_i]]]),
                 iter=model_settings$Niter[job_i],
                 chain=model_settings$Nchain[job_i],
                 thin=model_settings$Nthin[job_i],
                 warmup=model_settings$Nwarmup[job_i],
                 save_warmup = FALSE,
                 seed=job_i,
                 pars=c("a_alpha"), # we don't save this as it takes up lots of memory!
                 include=FALSE)
  
  
  save(out, file = paste0('Rout/01_Rout_temporal_splines_analysis_epoch/model_fits_',job_i,'.RData'))# save output
  
  writeLines('Finished job')


