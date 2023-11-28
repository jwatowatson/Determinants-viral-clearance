mod_time_dose <- gam(N_dose ~ s(Study_time, k = 8), 
                      data = Baseline_vl_data_Plot %>% filter(Site=='th001')) 

pred_time_dose <- predict_gam(mod_time_dose)
pred_time_dose

pred_time_dose$Rand_date <-  back_transfrom_date(platcov_dat, pred_time_dose$Study_time)


A <- pred_time_dose %>%
  ggplot(aes(Rand_date, fit)) +
  geom_point(data = Baseline_vl_data_Plot  %>% filter(Site=='th001'), aes(x = Rand_date, y = N_dose),
             size = 3, alpha = 0.03) +
  geom_smooth_ci(linewidth = 1, col = "#D80032", ci_alpha = 0.25) +
  theme_bw() +
  xlab("Randomisation date") +
  ylab("Number of vaccination doses") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8))
A
#######################################################################################
mod_viral_dose <- gam(log10_viral_load ~ s(N_dose, k = 5), 
                     data = Baseline_vl_data_Plot %>% filter(Site=='th001')) 

pred_viral_dose <- predict_gam(mod_viral_dose)

#pred_viral_dose$Rand_date <-  back_transfrom_date(platcov_dat, pred_viral_dose$Study_time)


B<- pred_viral_dose %>%
  ggplot(aes(N_dose, fit)) +
  geom_point(data = Baseline_vl_data_Plot  %>% filter(Site=='th001'), aes(x = N_dose, y = log10_viral_load),
             size = 3, alpha = 0.08) +
  geom_smooth_ci(linewidth = 1, col = "#D80032", ci_alpha = 0.25) +
  theme_bw() +
  ylab("Baseline SARS-CoV-2 genomes/mL") +
  xlab("Number of vaccination doses") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8))


ggarrange(A, B, ncol = 2)




