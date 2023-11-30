#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)
############################################################################
## Packages needed
library(rstan)
library(matrixStats)
library(doParallel)
library(tidyverse)

############################################################################
source('functions.R')

############################################################################
covs_base = c('Variant')

############################################################################
load('Rout/model_settings_bootstraps.RData') # change here for ineffective drugs

Max_job = nrow(model_settings)

for(job_i in 1:Max_job){
if(job_i > Max_job) stop('no model setting corresponding to job ID')

writeLines('Doing the following job:')
print(model_settings[job_i, ])

options(mc.cores = model_settings$Nchain[job_i])
stopifnot(model_settings$Nchain[job_i]>getDoParWorkers()) # check worker number assigned

mod = stan_model(file = as.character(model_settings$mod[job_i])) # compile 
############################################################################
platcov_dat_analysis <- data_list[[model_settings$data_ID[job_i]]]
Dmax = model_settings$Dmax[job_i]

ref_arm = model_settings$ref_arm[job_i]
trts = model_settings$intervention[job_i]

# only going to use data from main thai site (best quality data)
platcov_dat_analysis = platcov_dat_analysis %>%
  filter(Trt %in% c(ref_arm, trts), 
         Timepoint_ID <= Dmax, # timepoint is the day of follow-up
         Time < Dmax+1, # sample has to be taken at most 24 hours after last day
         mITT, Site=='th001') %>%
  mutate(Trt = factor(Trt, levels=c(ref_arm, trts)),
         Variant = as.factor(Variant),
         Site = as.factor(Site),
         RnaseP_scaled = scale(40 - CT_RNaseP,scale = F),
         Mean_age = mean(Age[!duplicated(ID)]),
         SD_age = sd(Age[!duplicated(ID)]),
         Age_scaled = (Age-Mean_age)/SD_age,
         Symptom_onset = ifelse(is.na(Symptom_onset),2,Symptom_onset)) 


## Bootstrap data (by bootstrapping patients)
ids_boot = sort(sample(x = unique(platcov_dat_analysis$ID), 
                       size = length(unique(platcov_dat_analysis$ID)),
                       replace = T))
platcov_dat_analysis  = platcov_dat_analysis %>%
  filter(ID %in% ids_boot) %>%
  group_by(ID) %>%
  mutate(ntimes = sum(ids_boot==ID[1])) %>% ungroup() %>%
  uncount(weights = ntimes, .id = 'ID_boot')
platcov_dat_analysis$ID = apply(platcov_dat_analysis[, c('ID','ID_boot')],1,function(x) paste(x,collapse = '_'))
platcov_dat_analysis = platcov_dat_analysis %>%
  arrange(ID, Time) %>%
  arrange(log10_viral_load==log10_cens_vl) 

stan_input_job = make_stan_inputs(input_data_fit = platcov_dat_analysis,
                                  int_covs_base = covs_base,
                                  int_covs_full = covs_base,
                                  slope_covs_base = covs_base,
                                  slope_covs_full = covs_base,
                                  trt_frmla = formula('~ Trt'),
                                  epoch = F,
                                  Dmax = Dmax+1)

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
               pars=c('trt_effect'), # only save trt effect parameter 
               include=T)

save(out, file = paste0('Rout/02 Rout bootstraps_analysis/model_fits_bootstraps_',job_i,'.RData'))# save output # change here for ineffective drugs

writeLines('Finished job')
}









