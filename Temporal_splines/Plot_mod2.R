load("../Rout/model_settings_all_analysis.RData")
################################################################################
library(rstan)
library(dplyr)
library(ggpubr)
library(grid) 
library(ggplot2)
library(DescTools)
library(scales)

################################################################################
back_transfrom_date <- function(platcov_dat, scaled_time){
  mean_time <- mean(platcov_dat$Study_time_normal)
  sd_time <- sd(platcov_dat$Study_time_normal)

  min_date <- as.numeric(min(as.Date(platcov_dat$Rand_date)))
  
  backward_time <- (scaled_time * sd_time) + mean_time
  
  date <- as.POSIXct.Date(backward_time * 7 + min_date)
  return(date)
}

new_variant <- function(platcov_dat){
  platcov_dat$Variant2 <- as.character(platcov_dat$Variant)
  platcov_dat$Variant2[platcov_dat$Variant2 %in% c("BA.5.2", "BA.5.5", "BQ.1")] <- "BA.5"
  platcov_dat$Variant2[platcov_dat$Variant2 %in% c("BN.1.2", "BN.1.3", "CH.1.1")] <- "BA.2.75"
  platcov_dat$Variant2[platcov_dat$Variant2 %in% c("XBB1.5-like with F456L")] <- "XBB.1.5-like"
  platcov_dat$Variant2 <- as.factor(platcov_dat$Variant2)
  
  platcov_dat
}

platcov_dat <- new_variant(platcov_dat)
################################################################################
source('../functions.R')
source('../priors.R')
################################################################################
Dmax <- 7


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


# platcov_dat_analysis$Symptom_onset[platcov_dat_analysis$ID == "PLT-TH1-627"] <- 2
################################################################################
res_dir <- "Rout"
list_files <- list.files(res_dir)

i <- 8
model_settings[i,]

Study_time_data_ALL <- NULL
Baseline_vl_data_ALL <- NULL

#for(i in 1:length(list_files)){
  load(paste0(res_dir, "/model_fits_",i,".RData"))
  
  #-----------------------------------------------------------------------
  #Alpha_hat
  post_alpha_hat <- apply(rstan::extract(out, "alpha_hat")[[1]],2,median)
  Study_time_data <- unique(platcov_dat_analysis[,c("ID", "Rand_date", "Study_time", "Study_time_normal", "N_dose")])
  Study_time_data$alpha_hat <- post_alpha_hat
  Study_time_data$alpha_hat_low <- apply(rstan::extract(out, "alpha_hat")[[1]],2,quantile, 0.025)
  Study_time_data$alpha_hat_high <- apply(rstan::extract(out, "alpha_hat")[[1]],2,quantile, 0.975)
  
  
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
  
  data_for_plot_slope <- data_for_plot_slope[data_for_plot_slope$beta_hat != max(data_for_plot_slope$beta_hat),]
  
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

Study_time_data_Plot <- Study_time_data_Plot[!Study_time_data_Plot$alpha_hat<1,]


library(RColorBrewer)

Baseline_vl_data_Plot$Variant2

Baseline_vl_data_Plot$Variant2 <- as.factor(Baseline_vl_data_Plot$Variant2)
Baseline_vl_data_Plot$Variant2 <- factor(Baseline_vl_data_Plot$Variant2, levels = as.character(Baseline_vl_data_Plot$Variant2[!duplicated(Baseline_vl_data_Plot$Variant2)]))
vars <- levels(Baseline_vl_data_Plot$Variant2)
Baseline_vl_data_Plot$Variant2 <- factor(Baseline_vl_data_Plot$Variant2, levels = c(vars[which(vars != "Other")], "Other"))

nb.cols <- length(levels(Baseline_vl_data_Plot$Variant2))
mycolors <-  c("#E41A1C", "#377EB8", "#4DAF4A", "black", "#984EA3", "#FF8400", "#2B3499", "#E95793", "#999999")  #colorRampPalette(palette.colors(palette = "set1"))(nb.cols) #hcl.colors(nb.cols, palette = "Sunset")      #colorRampPalette(palette.colors(palette = "Spectral"))(nb.cols)

A <- ggplot() + 
  geom_point(data = Baseline_vl_data_Plot, size = 2.5, 
             aes(x = Rand_date, y = log10_viral_load, col = Variant2, alpha = Variant2, shape = Variant2)) +
  theme_bw() +
  geom_ribbon(data = Study_time_data_Plot, aes(x = Rand_date, ymin = alpha_hat_low, ymax = alpha_hat_high), 
              alpha = 0.45, fill = "#27374D") +
  geom_line(data = Study_time_data_Plot, aes(x = Rand_date, y = alpha_hat), 
            linewidth = 1.5, alpha = 1, col = "black") +
  scale_color_manual(values = mycolors, name = "Variants") +
  xlab("Randomisation date") +
  ylab("Baseline SARS-CoV-2 genomes/mL") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom") +
  scale_alpha_manual(values = rep(0.4, length(vars)), guide = "none") +
  scale_shape_manual(values = rep(19, length(vars)), name = "Variants") +
  guides(col=guide_legend(nrow=2),
         shape = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(labels=label_math(), breaks = seq(0,9,1)) 
  
A
################################################################################

B <- ggplot(data_for_plot_slope, aes(x = Rand_date, y = slope_med)) +
  geom_point(size = 3, alpha = 0.5, col = "grey"#,
             #aes(col = Variant)
             ) +
  theme_bw() +
  geom_errorbar(aes(ymin = slope_low, ymax = slope_up), width = 0,
                alpha = 0.5, col = "grey") +
  geom_line(aes(x = Rand_date, y = beta_hat), 
            linewidth = 1.75, alpha = 0.7, col = "red") +
  xlab("Randomisation date") +
  ylab(expression(bold(paste("Viral clearance rates, ",alpha["0-7"]," (log"["10"]," genomes mL"^-1, " day"^-1,")")))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  ylim(-3,0)

B

png("slope_by_time_7d.png", width = 8, height = 6, units = "in", res = 400)
B
dev.off()

################################################################################
library(ggpubr)
ggarrange(A, B, ncol = 2, nrow = 1, labels = "AUTO", align = "hv")  
################################################################################
data_for_plot_slope$hl_med <- 24*log10(2)/(-(data_for_plot_slope$slope_med))
data_for_plot_slope$hl_up <- 24*log10(2)/(-(data_for_plot_slope$slope_up))
data_for_plot_slope$hl_low <- 24*log10(2)/(-(data_for_plot_slope$slope_low))

data_for_plot_slope$hl_trt_med <- 24*log10(2)/(-(data_for_plot_slope$slope_trt_med))



data_for_plot_slope$hl_hat <- 24*log10(2)/(-(data_for_plot_slope$beta_hat))
################################################################################
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
################################################################################
png("splines_half_life.png", width = 6, height = 5, units = "in", res = 400)
G3
dev.off()
################################################################################
write.csv(data_for_plot_slope, "estimated_slope.csv", row.names = F)
################################################################################
png("baseline_viral_loads.png", width = 8, height = 6, units = "in", res = 400)
A
dev.off()
################################################################################
data_for_plot_slope$Trt <- factor(
  data_for_plot_slope$Trt,
  levels = c(
    "Ivermectin",
    "Favipiravir",
    "Regeneron",
    "Remdesivir",
    "Molnupiravir",
    "Nirmatrelvir + Ritonavir",
    "No study drug"
  )
)

Trts <- unique(data_for_plot_slope$Trt)
Trts <- Trts[-which(Trts %in% c("No study drug"))]
Trts <- sort(Trts)
col <- c("#C70039", "#A03C78", "#FA7070", "#FF9130", "#1C6DD0", "#219C90")

plot_list <- list()
plot_list_hl <- list()

for(i in 1:length(Trts)){
  
  trt <- Trts[i]
  
  lab <- as.character(trt)
  if(lab == "Nirmatrelvir + Ritonavir"){lab <- "Nirmatrelvir"}
  if(lab == "Regeneron"){lab <- "Casirivimab/imdevimab"}
  
  
  G <- ggplot() +
    theme_bw() +
    geom_point(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, y = slope_med),
               size = 2.5, alpha = 0.5, col = "grey") +
    geom_errorbar(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, ymin = slope_low, ymax = slope_up), width = 0,
                  alpha = 0.5, col = "grey") +
    geom_line(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, y = slope_trt_med), 
              linewidth = 1.25, alpha = 0.7, col = "black") +
    
    
    
    geom_point(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, y = slope_med),
               size = 2.5, alpha = 0.25, col = col[i]) +
    geom_errorbar(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, ymin = slope_low, ymax = slope_up), width = 0,
                  alpha = 0.25, col = col[i]) +
    geom_line(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, y = slope_trt_med), 
              linewidth = 1.25, alpha = 0.9, col = col[i]) +
    
    
    xlab("Randomisation date") +
    ylab("Viral reduction per day (log10 units)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.title = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(1, "lines")) +
    ggtitle(lab) +
    ylim(-3,0)
  
  G2 <- ggplot() +
    theme_bw() +
    geom_point(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, y = hl_med),
               size = 2.5, alpha = 0.5, col = "grey") +
    geom_errorbar(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, ymin = hl_low, ymax = hl_up), width = 0,
                  alpha = 0.5, col = "grey") +
    geom_line(data = data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ], aes(x = Rand_date, y = hl_trt_med), 
              linewidth = 1.25, alpha = 0.7, col = "black") +
    
    
    
    geom_point(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, y = hl_med),
               size = 2.5, alpha = 0.25, col = col[i]) +
    geom_errorbar(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, ymin = hl_low, ymax = hl_up), width = 0,
                  alpha = 0.25, col = col[i]) +
    geom_line(data = data_for_plot_slope[data_for_plot_slope$Trt == trt, ], aes(x = Rand_date, y = hl_trt_med), 
              linewidth = 1.25, alpha = 0.9, col = col[i]) +
    
    
    xlab("Randomisation date") +
    ylab("Viral clearance half-life (hours)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(strip.text = element_text(size = 10, face = "bold"),
          axis.title = element_blank(),
          plot.title = element_text(size = 14, face = "bold"),
          panel.spacing = unit(1, "lines")) +
    ggtitle(lab) +
    ylim(0,50)
  
  plot_list[[i]] <- G
  plot_list_hl[[i]] <- G2 
  
}

################################################################################
GG <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2)
GG2 <- annotate_figure(GG, left = textGrob(expression(bold(paste("Viral clearance rates, ",alpha["0-5"]," (log"["10"]," genomes mL"^-1, " day"^-1,")"))), 
                                           rot = 90, vjust = 0.5, gp = gpar(cex = 1.1)),
                       bottom = textGrob(expression(bold("Randomisation date")), gp = gpar(cex = 1.1)))

png("slope_by_time_trt_5d.png", width = 10, height = 6, units = "in", res = 350)
GG2
dev.off()


################################################################################
GG_hl <- ggarrange(plotlist = plot_list_hl, ncol = 3, nrow = 2)
GG2_hl <- annotate_figure(GG_hl, left = textGrob("Viral clearance half-life (hours)", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Randomisation date", gp = gpar(cex = 1.3)))

png("hl_by_time_trt.png", width = 10, height = 6, units = "in", res = 350)
GG2_hl
dev.off()

################################################################################
library(ggpmisc)

rand_intercept <- apply(post_trt_theta_rand_id[,,1],2,median)
rand_slope <- apply(post_trt_theta_rand_id[,,2],2,median)

random_effects <- data.frame(rand_intercept, rand_slope)

ran_cor <- cor.test(rand_intercept, rand_slope)
ran_cor_med <- round(ran_cor$estimate,4)
ran_cor_low <- round(ran_cor$conf.int[1],4)
ran_cor_up <- round(ran_cor$conf.int[2],4)

lab_cor <- paste0("r = ", ran_cor_med, " [", ran_cor_low, "; ", ran_cor_up, "]")

png("random_effects.png", width = 8, height = 6, units = "in", res = 350)
ggplot(random_effects, aes(x = rand_intercept, y = rand_slope)) +
  geom_point( size = 3, alpha = 0.4) +
  theme_bw() +
  xlab("Random effects on baseline viral loads") +
  ylab("Random effects on average viral clearance rate") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10)) +
  geom_text(x=-3.1, y = 1, label = lab_cor, fontface = "plain") +
  stat_poly_line(col = "#E41A1C") +
  stat_poly_eq() 
dev.off()

  ?geom_text

