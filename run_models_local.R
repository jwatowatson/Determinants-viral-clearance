#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)


load(paste0('Rout/model_run_setup_',args[1],'.RData'))
source('functions.R')
writeLines(sprintf('Doing analysis for %s......', args[1]))

for(i in 1:nrow(model_settings)){
  
  writeLines('Doing the following job:')
  print(model_settings[i, ])
  
  options(mc.cores = model_settings$Nchain[i])
  stopifnot(model_settings$Nchain[i]>getDoParWorkers()) # check worker number assigned
  
  mod = stan_model(file = as.character(model_settings$mod[i])) # compile 
  
  stan_input_job = stan_inputs
  
  analysis_data_stan = stan_input_job$analysis_data_stan
  analysis_data_stan$trt_mat = stan_input_job$Trt_matrix
  analysis_data_stan$K_trt = ncol(analysis_data_stan$trt_mat)
  
  analysis_data_stan = make_cov_matrices_stan(stan_input_job, analysis_data_stan,i)
  
  # sample posterior
  out = sampling(mod, 
                 data=c(analysis_data_stan,
                        all_priors[[model_settings$prior[i]]]),
                 iter=model_settings$Niter[i],
                 chain=model_settings$Nchain[i],
                 thin=model_settings$Nthin[i],
                 warmup=model_settings$Nwarmup[i],
                 save_warmup = FALSE,
                 seed=i,
                 pars=c('L_Omega'), # we don't save this as it takes up lots of memory!
                 include=FALSE)
  
  
  save(out, file = paste0('Rout/model_fits_',i,'_',args[1],'.RData'))# save output
  
  writeLines('Finished job')
}
