load("Rout/model_settings_temporal_splines_mod2.RData")
################################################################################
library(rstan)
library(dplyr)
################################################################################
source('../functions.R')
source('../priors.R')
################################################################################
Dmax <- 7

#platcov_dat <- platcov_dat[platcov_dat$Trt %in% c("No study drug", "Ivermectin"),]

# Analysis data
platcov_dat_analysis = 
  platcov_dat %>% ungroup() %>%
  filter(Time <= Dmax, mITT) %>%
  arrange(log10_viral_load==log10_cens_vl) %>%
  mutate(Variant = as.factor(Variant),
         Site = as.factor(Site),
         RnaseP_scaled = scale(40 - CT_RNaseP,scale = F),
         Mean_age = mean(Age[!duplicated(ID)]),
         SD_age = sd(Age[!duplicated(ID)]),
         Age_scaled = (Age-Mean_age)/SD_age,
         Symptom_onset = ifelse(is.na(Symptom_onset),2,Symptom_onset)) 

covs_base = c('Site') #'Study_time'
covs_full=c(covs_base, 'Age_scaled','Symptom_onset','N_dose')

stan_input_job = 
  make_stan_inputs(input_data_fit = platcov_dat_analysis,
                   int_covs_base = covs_base,
                   int_covs_full = covs_full,
                   slope_covs_base = covs_base,
                   slope_covs_full = covs_full,
                   trt_frmla = formula('~ Trt'),
                   Dmax = Dmax)

n_id <- stan_input_job$analysis_data_stan$n_id
ind_start <- stan_input_job$analysis_data_stan$ind_start
################################################################################
res_dir <- "Rout_mod2_splines"
list_files <- list.files(res_dir)

i <- 4

Study_time_data_ALL <- NULL
Baseline_vl_data_ALL <- NULL

#for(i in 1:length(list_files)){
  load(paste0(res_dir, "/model_fits_",i,".RData"))
  
  #-----------------------------------------------------------------------
  #Alpha_hat
  post_alpha_hat <- apply(rstan::extract(out, "alpha_hat")[[1]],2,median)
  Study_time_data <- unique(platcov_dat_analysis[,c("ID", "Rand_date", "Study_time", "Study_time_normal")])
  Study_time_data$alpha_hat <- post_alpha_hat
  
  Study_time_data$num_knots_alpha <- model_settings$num_knots_alpha[i] 
  Study_time_data$spline_degree_alpha <- model_settings$spline_degree_alpha[i] 
  #Baseline viral loads
  Baseline_vl_data <- platcov_dat_analysis[platcov_dat_analysis$Timepoint_ID == 0,]
  Baseline_vl_data$num_knots_alpha <- model_settings$num_knots_alpha[i]
  Baseline_vl_data$spline_degree_alpha <- model_settings$spline_degree_alpha[i] 
  
  Study_time_data_ALL <- rbind(Study_time_data_ALL, Study_time_data)
  Baseline_vl_data_ALL <- rbind(Baseline_vl_data_ALL, Baseline_vl_data)
  #-----------------------------------------------------------------------
  post_beta_hat <- rstan::extract(out, "beta_hat")[[1]]
  post_trt_slope <- rstan::extract(out, "trt_slope")[[1]]
  post_trt_theta_rand_id <- rstan::extract(out, "theta_rand_id")[[1]]
  post_trt_beta_cov <- rstan::extract(out, "beta_cov")[[1]]
  
  slope_ALL <- matrix(NA, ncol = n_id, nrow = 1000)
  slope_trt <- matrix(NA, ncol = n_id, nrow = 1000)
  
  
  for(j in 1:n_id){
    slope_ALL[,j]  <- post_beta_hat[,j] * exp(post_trt_slope[,ind_start[j]] +  post_trt_theta_rand_id[,j,2] +  post_trt_beta_cov[,ind_start[j]])
    slope_trt[,j] <-  post_beta_hat[,j] * exp(post_trt_slope[,ind_start[j]])
  }
  
  
  slope_summarize <- apply(slope_ALL, 2, quantile, c(0.025, 0.5, 0.975))
  slope_trt_summarize <- apply(slope_trt, 2, quantile, c(0.025, 0.5, 0.975))
  post_beta_hat_summarize <- apply(post_beta_hat, 2, quantile, c(0.025, 0.5, 0.975))
  
  
  
  data_for_plot_slope <- platcov_dat_analysis[ind_start,]
  data_for_plot_slope$slope_low <- slope_summarize[1,]
  data_for_plot_slope$slope_med <- slope_summarize[2,]
  data_for_plot_slope$slope_up <- slope_summarize[3,]
  
  data_for_plot_slope$slope_trt_low <- slope_trt_summarize[1,]
  data_for_plot_slope$slope_trt_med <- slope_trt_summarize[2,]
  data_for_plot_slope$slope_trt_up <- slope_trt_summarize[3,]
  
  data_for_plot_slope$beta_hat <- post_beta_hat_summarize[2,]
  
#}
################################################################################
library(ggplot2)

Baseline_vl_data_ALL$num_knots_alpha <- as.factor(Baseline_vl_data_ALL$num_knots_alpha)
levels(Baseline_vl_data_ALL$num_knots_alpha) <- paste0("K = ", levels(Baseline_vl_data_ALL$num_knots_alpha))

Study_time_data_ALL$num_knots_alpha <- as.factor(Study_time_data_ALL$num_knots_alpha)
levels(Study_time_data_ALL$num_knots_alpha) <- paste0("K = ", levels(Study_time_data_ALL$num_knots_alpha))

num_knots_plot <- levels(Baseline_vl_data_ALL$num_knots_alpha)[c(1,3,5,7)]

Baseline_vl_data_Plot <- Baseline_vl_data_ALL[Baseline_vl_data_ALL$num_knots_alpha %in% num_knots_plot, ]
Study_time_data_Plot <- Study_time_data_ALL[Study_time_data_ALL$num_knots_alpha %in% num_knots_plot, ]

A <- ggplot(Baseline_vl_data_Plot, aes(x = Rand_date, y = log10_viral_load)) + 
  geom_point(size = 3, alpha = 0.2, aes(col = Variant)) +
  theme_bw() +
 # facet_wrap(num_knots_alpha ~ ., ncol = 2, nrow = 2) +
  geom_line(data = Study_time_data_Plot, aes(x = Rand_date, y = alpha_hat), 
            linewidth = 1.25, alpha = 0.7, col = "red") +
#  scale_color_manual(values = c("#974EC3", "#F94C10", "#0D1282", "#C70039"), name = "Spline degree") +
  xlab("Randomisation date") +
  ylab("Log10 SARS-CoV-2 genomes/mL") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle("Baseline viral loads")

A
################################################################################
data_for_plot_slope

B <- ggplot(data_for_plot_slope, aes(x = Rand_date, y = slope_med)) +
  geom_point(size = 3, alpha = 0.5, col = "grey"#,
             #aes(col = Variant)
             ) +
  theme_bw() +
  geom_errorbar(aes(ymin = slope_low, ymax = slope_up), width = 0,
                alpha = 0.5, col = "grey") +
  geom_line(aes(x = Rand_date, y = beta_hat), 
            linewidth = 1.25, alpha = 0.7, col = "red") +
  xlab("Randomisation date") +
  ylab("Viral reduction per day (log10 units)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines"))

B
################################################################################
library(ggpubr)
ggarrange(A, B, ncol = 2, nrow = 1, labels = "AUTO", align = "hv")  
################################################################################
data_for_plot_slope$hl_med <- 24*log10(2)/(-(data_for_plot_slope$slope_med))
data_for_plot_slope$hl_up <- 24*log10(2)/(-(data_for_plot_slope$slope_up))
data_for_plot_slope$hl_low <- 24*log10(2)/(-(data_for_plot_slope$slope_low))

data_for_plot_slope$hl_hat <- 24*log10(2)/(-(data_for_plot_slope$beta_hat))




G3 <- ggplot(data_for_plot_slope, aes(x = Rand_date, y = hl_med)) +
  geom_point(size = 3, alpha = 0.25, col = "grey",
            ) +
  theme_bw() +
  geom_errorbar(aes(ymin = hl_low, ymax = hl_up), width = 0,
                alpha = 0.25, col = "grey") +
  geom_line(aes(x = Rand_date, y = hl_hat), 
            linewidth = 1.25, alpha = 0.7, col = "red") +
  xlab("Randomisation date") +
  ylab("Clearance half-life (hours)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ylim(0,60)


png("splines_half_life.png", width = 6, height = 5, units = "in", res = 400)
G3
dev.off()


write.csv(data_for_plot_slope, "estimated_slope.csv", row.names = F)

png("baseline_viral_loads.png", width = 8, height = 5, units = "in", res = 400)
A
dev.off()



