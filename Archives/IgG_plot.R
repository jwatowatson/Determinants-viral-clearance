IDs <- unique(platcov_dat$ID)

sero_unblinded_bl <- output_df_unk_unique %>% 
                     filter(ID %in% IDs, Day == 0) 

rand_date <- unique(platcov_dat[,c("ID", "Rand_date", "VariantClass_new", "N_dose", "Time_since_last_dose",
                                   "Study_time")])


sero_unblinded_bl_merged <- merge(rand_date, sero_unblinded_bl, by = "ID")


mod <- gam(Med_log_IgG ~ s(Study_time), 
            data = sero_unblinded_bl_merged) 

pred_mod <- predict_gam(mod)
pred_mod$Rand_date <-  back_transfrom_date(platcov_dat, pred_mod$Study_time)

library(ggplot2)
vars <- levels(sero_unblinded_bl_merged$VariantClass_new)

sero_unblinded_bl_merged$VariantClass_new <- factor(sero_unblinded_bl_merged$VariantClass_new, levels = vars[c(6, 1:5,8,9,7)])
mycolors <-  c("#E41A1C", "#377EB8", "#4DAF4A", "black", "#984EA3", "#FF8400", "#2B3499", "#E95793", "#999999")



G <- pred_mod %>%
  ggplot(aes(Rand_date, fit)) +
  geom_point(data = sero_unblinded_bl_merged, aes(x = Rand_date, y = Med_log_IgG, col = VariantClass_new,
                                                  alpha = VariantClass_new, shape = VariantClass_new), 
             size = 2.5) +
  # geom_errorbar(data = sero_unblinded_bl_merged, aes(x = Rand_date, ymin = Low_log_IgG, ymax = Up_log_IgG), 
  #               width = 0, alpha = 0.6, linewidth = 0.1) +
  geom_smooth_ci(linewidth = 1.5, col = "black", ci_alpha = 0.25)  +
  theme_bw() +
  scale_color_manual(values = mycolors, name = "Variants")  +
  scale_y_continuous(labels=label_math(), breaks = seq(-5,5,1)) +
  xlab("Randomisation date") +
  ylab("Median baseline IgG estimate (log10)") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "bottom") +
  scale_alpha_manual(values = rep(0.6, length(vars)), guide = "none") +
  scale_shape_manual(values = rep(19, length(vars)), name = "Variants") +
  guides(col=guide_legend(nrow=2),
         shape = guide_legend(override.aes = list(size = 4)))
  
png("Plots/Figx_baseline_IgG.png", width = 8, height = 6, units = "in", res = 400)
G
dev.off()







ggplot(sero_unblinded_bl_merged, aes(x = Rand_date, y = Med_log_IgG, col = VariantClass_new)) +
  geom_errorbar(aes(ymin = Low_log_IgG, ymax = Up_log_IgG), width = 0, alpha = 0.6, linewidth = 0.1) +
  geom_point(size = 2.5, alpha = 0.6) +
  theme_bw() +
  scale_color_manual(values = mycolors, name = "Variants")  +
  scale_y_continuous(labels=label_math(), breaks = seq(-5,5,1)) +
  xlab("Randomisation date") +
  ylab("Baseline SARS-CoV-2 genomes/mL") 
  





IDs[which(! IDs %in% sero_unblinded_bl$ID)]
sero_unblinded_bl
unique(sero_unblinded_bl$ID)
