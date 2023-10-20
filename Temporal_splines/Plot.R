load("Rout/model_settings_temporal_splines.RData")
################################################################################
library(rstan)
library(dplyr)
################################################################################
Dmax <- 7

platcov_dat <- platcov_dat[platcov_dat$Trt %in% c("No study drug", "Ivermectin"),]

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
################################################################################
res_dir <- "Rout_mod1_spline_intercept"
list_files <- list.files(res_dir)

i <- 1

Study_time_data_ALL <- NULL
Baseline_vl_data_ALL <- NULL

for(i in 1:length(list_files)){
  load(paste0(res_dir, "/model_fits_",i,".RData"))
  
  post_alpha_hat <- apply(rstan::extract(out, "alpha_hat")[[1]],2,median)
  Study_time_data <- unique(platcov_dat_analysis[,c("ID", "Rand_date", "Study_time", "Study_time_normal")])
  Study_time_data$alpha_hat <- post_alpha_hat
  
  Study_time_data$num_knots_alpha <- model_settings$num_knots_alpha[i] 
  Study_time_data$spline_degree_alpha <- model_settings$spline_degree_alpha[i] 
  
  
  Baseline_vl_data <- platcov_dat_analysis[platcov_dat_analysis$Timepoint_ID == 0,]
  Baseline_vl_data$num_knots_alpha <- model_settings$num_knots_alpha[i]
  Baseline_vl_data$spline_degree_alpha <- model_settings$spline_degree_alpha[i] 
  
  Study_time_data_ALL <- rbind(Study_time_data_ALL, Study_time_data)
  Baseline_vl_data_ALL <- rbind(Baseline_vl_data_ALL, Baseline_vl_data)
  
}
################################################################################
library(ggplot2)

Baseline_vl_data_ALL$num_knots_alpha <- as.factor(Baseline_vl_data_ALL$num_knots_alpha)
levels(Baseline_vl_data_ALL$num_knots_alpha) <- paste0("K = ", levels(Baseline_vl_data_ALL$num_knots_alpha))

Study_time_data_ALL$num_knots_alpha <- as.factor(Study_time_data_ALL$num_knots_alpha)
levels(Study_time_data_ALL$num_knots_alpha) <- paste0("K = ", levels(Study_time_data_ALL$num_knots_alpha))

num_knots_plot <- levels(Baseline_vl_data_ALL$num_knots_alpha)[c(1,3,5,7)]

Baseline_vl_data_Plot <- Baseline_vl_data_ALL[Baseline_vl_data_ALL$num_knots_alpha %in% num_knots_plot, ]
Study_time_data_Plot <- Study_time_data_ALL[Study_time_data_ALL$num_knots_alpha %in% num_knots_plot, ]

ggplot(Baseline_vl_data_Plot, aes(x = Rand_date, y = log10_viral_load)) + geom_point(size = 3, alpha = 0.1, col = "grey") +
  theme_bw() +
  facet_wrap(num_knots_alpha ~ ., ncol = 2, nrow = 2) +
  geom_line(data = Study_time_data_Plot, aes(x = Rand_date, y = alpha_hat, col = as.factor(spline_degree_alpha)), linewidth = 1.25, alpha = 0.7) +
  scale_color_manual(values = c("#974EC3", "#F94C10", "#0D1282", "#C70039"), name = "Spline degree") +
  xlab("Randomisation date") +
  ylab("Log10 SARS-CoV-2 genomes/mL") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ggtitle("Baseline viral loads")


