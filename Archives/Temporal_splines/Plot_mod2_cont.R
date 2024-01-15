library(mgcv)
library(tidymv)
library(splines)

source("Plot_mod2.R")

Baseline_vl_data_Plot$Symptom_onset_round <- round(Baseline_vl_data_Plot$Symptom_onset)
Baseline_vl_data_Plot$Symptom_onset_round[Baseline_vl_data_Plot$Symptom_onset_round == 0] <- 1


Baseline_vl_data_Plot$Symptom_onset_round_factor <- as.factor(Baseline_vl_data_Plot$Symptom_onset_round)
levels(Baseline_vl_data_Plot$Symptom_onset_round_factor)  <- c("1 day", "2 days", "3 days", "4 days")

mean(unique(Baseline_vl_data_Plot[,c("ID", "Symptom_onset")])$Symptom_onset)


formatC(10^mean(Baseline_vl_data_Plot$log10_viral_load), format = "e", digits = 2)

cor.test(Baseline_vl_data_Plot$Symptom_onset, Baseline_vl_data_Plot$log10_viral_load)
lm_baseline_mod <- lm(log10_viral_load ~ Symptom_onset, data = Baseline_vl_data_Plot)
summary(lm_baseline_mod)$r.squared

mycolors <-  c("#E41A1C", "#377EB8", "#4DAF4A", "black", "#984EA3", "#FF8400", "#2B3499", "#E95793", "#999999")  #colorRampPalette(palette.colors(palette = "set1"))(nb.cols) #hcl.colors(nb.cols, palette = "Sunset")      #colorRampPalette(palette.colors(palette = "Spectral"))(nb.cols)



########################################################################################
# Temporal trend of time since symptom onset
mod_time_onset <- gam(Symptom_onset_round ~ s(Study_time), 
                      data = Baseline_vl_data_Plot) 

pred_time_onset <- predict_gam(mod_time_onset)
pred_time_onset

pred_time_onset$Rand_date <-  back_transfrom_date(platcov_dat, pred_time_onset$Study_time)

omicron_date <- min(Baseline_vl_data_Plot$Rand_date[Baseline_vl_data_Plot$Variant2 == "BA.1"])
omicron_date_2 <- max(Baseline_vl_data_Plot$Rand_date[Baseline_vl_data_Plot$Variant2 == "BA.1"])

A <-pred_time_onset %>%
  ggplot(aes(Rand_date, fit)) +
  #geom_rect(xmin = omicron_date, xmax = omicron_date_2, ymin = -Inf, ymax = Inf, fill = "#190482", alpha = 0.005) +
  geom_vline(xintercept = omicron_date, linetype = "dashed", linewidth = 0.75) +
  geom_point(data = Baseline_vl_data_Plot, aes(x = Rand_date, y = Symptom_onset_round),
             size = 2, alpha = 0.06, col = "#435585") +
  geom_smooth_ci(linewidth = 1.5, col = "#D80032", ci_alpha = 0.25) +
  theme_bw() +
  xlab("Randomisation date") +
  ylab("Time since symptomp onset (days)") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8)) 
A
########################################################################################
# Effects of time since symptom onset on baseline viral loads
mod_onset_baseline <- gam(log10_viral_load ~ s(Symptom_onset_round, k=4), 
                      data = Baseline_vl_data_Plot) 

pred_onset_baseline <- predict_gam(mod_onset_baseline)
pred_onset_baseline

B <- pred_onset_baseline %>%
  ggplot(aes(Symptom_onset_round, fit)) +
  
 # scale_fill_manual(values = c("#D0A2F7", "#D0A2F7","#E5D4FF", "#F1EAFF"), guide = "none") +
   geom_jitter(data = Baseline_vl_data_Plot, aes(x = Symptom_onset_round_factor, y = log10_viral_load),
                size = 2, alpha = 0.3, width = 0.3, col = "#435585") +
  
  geom_violin(data = Baseline_vl_data_Plot, aes(x = Symptom_onset_round_factor, y = log10_viral_load),
              trim=T, fill='#8CC0DE', linewidth = 0.75, width = 0.75, alpha = 0.5) +
  geom_boxplot(data = Baseline_vl_data_Plot, aes(x = Symptom_onset_round_factor, y = log10_viral_load),
               width=0.1, size = 0.75,outlier.shape = NA,  coef = 0, alpha = 0.5) +
  geom_smooth_ci(linewidth = 1.5, col = "#D80032", ci_alpha = 0.25) +
  theme_bw() +
  xlab("Time since symptom onset (days)") +
  ylab("Baseline SARS-CoV-2 genomes/mL") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8)) +
  scale_y_continuous(labels=label_math(), breaks = seq(0,9,2)) 
B


#################################################################################################
png("time_onset_dynamics.png", width = 10, height = 5, units = "in", res = 350)
ggarrange(A, B, ncol = 2, labels = "AUTO", align = "hv")
dev.off()
#################################################################################################

mod_all <- gam(log10_viral_load ~ s(Symptom_onset_round, k=4) + 
                                  s(Study_time) +
                                  s(Study_time, by = Symptom_onset_round), 
                                  data = Baseline_vl_data_Plot) 

summary(mod_all)

preds <- predict_gam(mod_all, values = list(Symptom_onset_round = c(1,2,3,4)))
preds$Rand_date <-  back_transfrom_date(platcov_dat, preds$Study_time)

preds$Symptom_onset_round <- as.factor(preds$Symptom_onset_round)
levels(preds$Symptom_onset_round)  <- c("1 day", "2 days", "3 days", "4 days")

Baseline_vl_data_Plot$Symptom_onset_round <- as.factor(Baseline_vl_data_Plot$Symptom_onset_round)
levels(Baseline_vl_data_Plot$Symptom_onset_round)  <- c("1 day", "2 days", "3 days", "4 days")

Baseline_vl_data_Plot$Rand_date

library(ggnewscale)

mycolors <-  c("#E41A1C", "#377EB8", "#4DAF4A", "black", "#984EA3", "#FF8400", "#2B3499", "#E95793", "#999999")  #colorRampPalette(palette.colors(palette = "set1"))(nb.cols) #hcl.colors(nb.cols, palette = "Sunset")      #colorRampPalette(palette.colors(palette = "Spectral"))(nb.cols)


C <- preds %>%
  ggplot(aes(Rand_date, fit)) +
  geom_point(data = Baseline_vl_data_Plot, 
             aes(x = Rand_date, y = log10_viral_load, col = Variant2, alpha = Variant2),
             size = 2) +
  scale_color_manual(values = mycolors, name = "Variants") +
  facet_wrap(.~ Symptom_onset_round, ncol = 4) +
  guides(col = guide_legend(override.aes = list(size = 4))) +
  ggnewscale::new_scale_color() +
  geom_smooth_ci(Symptom_onset_round,  linewidth = 1, col = "black", ci_alpha = 0.25) +
  #scale_color_viridis(discrete = T, option="magma") +
  scale_linetype_manual(values = rep("solid",5), guide = "none") +
  theme_bw() +
  xlab("Randomisation date") +
  ylab("Baseline SARS-CoV-2 genomes/mL") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8),
        legend.position = "bottom",
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold")) +
  scale_alpha_manual(values = rep(0.3, length(vars)), guide = "none") +
  ggtitle("Temporal dynamics of baseline viral loads by time since symptom onset") +
  scale_y_continuous(labels=label_math(), breaks = seq(0,9,2)) 

C
#################################################################################################
png("baseline_dynamics_by_onset.png", width = 8, height = 6, units = "in", res = 350)
C
dev.off()
#################################################################################################

D <- ggarrange(ggarrange(B, NULL, A, ncol = 3, labels = c("A","", "B"), align = "hv", widths = c(1,0.05,1)),
               NULL,
               C,
               nrow = 3, labels = c("", "", "C"), heights = c(1, 0.1, 1.2)
)

png("onset_dynamics_all.png", width = 10, height = 8, units = "in", res = 350)
D
dev.off()



