library(rstan)
getwd()
#####################################################################
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
#####################################################################
Arms <- c("Nirmatrelvir + Ritonavir", "No study drug")
Mods <- c("Linear", "TwoCmpt", "UpDown", "Biexp")

Pairs <- as.data.frame(expand.grid(Arms, Mods))
colnames(Pairs) <- c("Arms", "Mods")

i <- 2
for(i in 1:8){
  Arm <- Pairs[i,]$Arms
  if(Arm == "Nirmatrelvir + Ritonavir"){ID_inc <- 1:30; Arm_lab <- "NR"}else{ID_inc <- c(1:24, 26:31); Arm_lab <- "NS"}
  #####################################################################
  data <- read.csv("Analysis_Data/drug_test_NS_NR_F.csv")
  data <- data[data$Trt == Arm,]
  data$ID <- as.numeric(as.factor(data$ID))
  data <- data[data$ID %in% ID_inc ,]
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
  #####################################################################
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
  #####################################################################
  Mod <- Pairs[i,]$Mods
  model <- stan_model(paste0("Stan_models/", Mod, "_for_ineffective_arm.stan"),verbose = T)

  print(paste0("0", i, "_", Mod, "_", Arm_lab))

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

  save(fit, file = paste0("Fit/0", i, "_", Mod, "_", Arm_lab, ".RData"))
}








