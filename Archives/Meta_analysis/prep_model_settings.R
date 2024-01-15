source('../functions.R')
source('../priors.R')
#######################################################################
library(RColorBrewer)
library(dplyr)
#######################################################################
intervention = c("Unblinded_meta") # prefix of analysis file
ref_arm = 'No study drug'
#######################################################################
data_list <- list()
for(i in 1:length(intervention)){
  
  f_name = paste0("../Analysis_Data/", intervention[i],'_analysis.csv')
  platcov_dat = read.csv(f_name)
  platcov_dat$Rand_date = as.POSIXct(platcov_dat$Rand_date)
  trt_intervention = unique(platcov_dat$Trt)
  trts = trt_intervention[trt_intervention!=ref_arm] # get interventions
  #######################################################################
  platcov_dat = platcov_dat %>% group_by(ID) %>%
    mutate(
      mITT = any(Per_protocol_sample==1 & Timepoint_ID>=3) &
        !all(CT_NS==40))
  #######################################################################
  platcov_dat = platcov_dat %>% ungroup() %>%
    mutate(Study_time = as.numeric(difftime(Rand_date,min(Rand_date),units = 'weeks')),
           Study_time = scale(Study_time) ) %>%
    group_by(ID, Timepoint_ID) %>%
    mutate(daily_VL = mean(log10_viral_load),
           Sex = as.factor(ifelse(Sex==1,'Male','Female')),
           Site = as.factor(Site),
           Trt = factor(Trt, levels=c(ref_arm, trts)),
           Vaccinated = as.factor(ifelse(N_dose>0,'Yes','No')),
           Variant = as.factor(Variant)#normalise
    )  %>%
    ungroup() %>%
    mutate(trt_color = brewer.pal(name = 'Dark2',8)[c(1,7)][as.numeric(Trt)]) 
  #######################################################################
  data_list[[i]] <- platcov_dat
}

Dmax <- c(5,7)
Data_ID <- 1:length(intervention)
mod <- c("../Stan_models/Linear_model_RNaseP.stan", "../Stan_models/Nonlinear_model_RNaseP.stan")

model_settings <-  unique(do.call(expand.grid, list("Dmax" = Dmax,
                                                    "Data_ID" = Data_ID, 
                                                    "mod" = mod)))

model_settings$intervention <- sapply(model_settings$Data_ID, FUN = function(x) intervention[x])
model_settings$ref_arm <- ref_arm
model_settings$prior <- 1
model_settings$cov_matrices <- 1
model_settings$Niter <- 2000
model_settings$Nwarmup <- 1000
model_settings$Nthin <- 4
model_settings$Nchain <- 4
#######################################################################
save(data_list,
     model_settings,
     all_priors,
     file = "../Rout/model_settings_Unblinded.RData")
#######################################################################