library(rstan)
getwd()


mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
data <- read.csv("Analysis_Data/drug_test_NS_NR_F.csv")
# choose drug
data <- data[data$Trt == "Nirmatrelvir + Ritonavir",] #"Nirmatrelvir + Ritonavir", ]
data$ID_code <- data$ID
data$ID <- as.numeric(as.factor(data$ID))
# choose the first 30
data <- data[data$ID %in% 1:30 ,] #,] #c(1:24, 26:31)
data$ID <- as.numeric(as.factor(data$ID))


data <- data[order(c(data$censor), decreasing = T),]

#Preparing data
Ntot <- nrow(data)
N_obs <- sum(data$censor == "none")
n_id <- length(unique(data$ID))
id <- data$ID
log_10_vl <- data$log10_viral_load
log10_cens_vl <- data$log10_cens_vl
Time_max <- max(data$Timepoint_ID)
obs_day <- data$Timepoint_ID
RNaseP <- as.vector(scale(data$CT_RNaseP))

alpha_0_prior_mean <- 5
alpha_0_prior_sd <- 2
beta_0_prior_mean <- -0.5
beta_0_prior_sd <- 1
sigma_logvl_mean <- 1
sigma_logvl_sd <- 1


data_for_stan <- list(
  #data
  Ntot = Ntot,
  N_obs = N_obs,
  n_id = n_id,
  id = id,
  log_10_vl = log_10_vl,
  log10_cens_vl = log10_cens_vl,
  Time_max = Time_max,
  obs_day = obs_day,
  RNaseP = RNaseP,
  
  #priors
  alpha_0_prior_mean = alpha_0_prior_mean,
  alpha_0_prior_sd = alpha_0_prior_sd,
  beta_0_prior_mean = beta_0_prior_mean,
  beta_0_prior_sd = beta_0_prior_sd,
  sigma_logvl_mean = sigma_logvl_mean,
  sigma_logvl_sd = sigma_logvl_sd
  
)
#############################################################################
# Running stan code
model <- stan_model("Stan_models/Biexp_for_ineffective_arm.stan",verbose = T)

fit <- sampling(
  object = model,         # Stan model
  data = data_for_stan,            # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  thin = 2,
  cores =4,              # number of cores (could use one per chain)
  refresh = 100
)

fit
save(fit, file = "Fit/Biexp_NR_only.RData")

names(fit)
traceplot(fit, inc_warmup  = T,  pars = c("logalpha_0", "logbeta_0", "A0_0", "B0_0"))
