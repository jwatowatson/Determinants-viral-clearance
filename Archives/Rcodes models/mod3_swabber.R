library(rstan)
getwd()


mainDir <- "/Users/pwongnak/Project/Determinants-viral-clearance"  # "D:/Determinants-viral-clearance"
setwd(mainDir)
data <- read.csv("Analysis_Data/swabber_analysis.csv")
data <- data[!data$Swabber == "-",]
data$Swabber_key <- as.factor(data$Swabber)
data$Swabber <- as.numeric(data$Swabber_key)

data$ID_code <- data$ID
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

n_swabber <- length(unique(data$Swabber))
swabber_id <- data$Swabber

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
  n_swabber = n_swabber,
  swabber_id = swabber_id,
  
  #priors
  #alpha_0_prior_mean = alpha_0_prior_mean,
  #alpha_0_prior_sd = alpha_0_prior_sd,
  #beta_0_prior_mean = beta_0_prior_mean,
  #beta_0_prior_sd = beta_0_prior_sd#,
  sigma_logvl_mean = sigma_logvl_mean,
  sigma_logvl_sd = sigma_logvl_sd
  
)
#############################################################################
# Running stan code
model <- stan_model("Stan_models/TwoCmpt_for_ineffective_arm_swabber.stan",verbose = T)

fit <- sampling(
  object = model,         # Stan model
  data = data_for_stan,            # named list of data
  chains = 4,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 4000,            # total number of iterations per chain
  thin = 2,
  cores =4,              # number of cores (could use one per chain)
  refresh = 200
)

summary(fit, par = "sigma_logvl", prob = c(0.025, 0.5, 0.975))


save(fit, file = "Fit/TwoCmpt_nRNaseP_swabber.RData")

traceplot(fit, pars = c("sigma_logvl")) + ylim(0, 2) +
  geom_hline(yintercept = 1, col = "red", linetype = "dashed", linewidth = 0.75)


