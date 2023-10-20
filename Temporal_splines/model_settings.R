source('../functions.R')
source('../priors.R')
#######################################################################
library(RColorBrewer)
library(tidyverse)
library(dplyr)
#######################################################################
platcov_dat = read.csv("../Analysis_Data/interim_all_analysis.csv")
platcov_dat <- platcov_dat %>% filter(Trt %in% c("No study drug", "Ivermectin",
                                                 "Regeneron",
                                                 "Remdesivir",
                                                 "Favipiravir", 
                                                 "Molnupiravir", 
                                                 "Nirmatrelvir + Ritonavir"))
platcov_dat$Rand_date = as.POSIXct(platcov_dat$Rand_date)
trt_intervention = unique(platcov_dat$Trt)
ref_arm = "No study drug"
trts = trt_intervention[trt_intervention!=ref_arm] # get interventions
#######################################################################
platcov_dat = platcov_dat %>% group_by(ID) %>%
  mutate(
    mITT = any(Per_protocol_sample==1 & Timepoint_ID>=3) &
      !all(CT_NS==40))
#######################################################################
platcov_dat = platcov_dat %>% ungroup() %>%
  mutate(Study_time = as.numeric(difftime(Rand_date,min(Rand_date),units = 'weeks')),
         Study_time_normal = Study_time,
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
Dmax <- c(7)
mod <- c(#"../Stan_models/Temporal_spline_mod1.stan")
         "../Stan_models/Temporal_spline_mod2_w_slope.stan"
         )
num_knots_alpha <- seq(4,10,1)
spline_degree_alpha <- seq(1,4,1)
num_knots_beta <- seq(4,10,1)
spline_degree_beta <- seq(1,4,1)
model_settings <-  unique(do.call(expand.grid, list("Dmax" = Dmax,
                                                    "mod" = mod,
                                                    "num_knots_alpha" = num_knots_alpha,
                                                    "spline_degree_alpha" = spline_degree_alpha,
                                                    "num_knots_beta" = num_knots_beta,
                                                    "spline_degree_beta" = spline_degree_beta
)))
model_settings$ref_arm <- ref_arm
model_settings$prior <- 1
model_settings$cov_matrices <- 1
model_settings$Niter <- 2000
model_settings$Nwarmup <- 1000
model_settings$Nthin <- 4
model_settings$Nchain <- 4
#######################################################################
save(platcov_dat,
     model_settings,
     all_priors,
     file = "Rout/model_settings_temporal_splines_mod2.RData")
#######################################################################