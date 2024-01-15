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
col <- c("#C70039", "#A03C78", "#FA7070", "#FF9130", "#1C6DD0", "#219C90", "black")
################################################################################
mod_viral_dose <- gam(slope_med ~ s(N_dose, k=6) + 
                        s(Study_time) +
                        s(Study_time, by = N_dose), 
                      data = data_for_plot_slope %>% filter(Site=='th001')) 

pred_viral_dose <- predict_gam(mod_viral_dose, values = list(N_dose = 0:5))

pred_viral_dose$Rand_date <-  back_transfrom_date(platcov_dat, pred_viral_dose$Study_time)

#pred_viral_dose$Rand_date <-  back_transfrom_date(platcov_dat, pred_viral_dose$Study_time)

data_for_plot_slope$Rand_date

pred_viral_dose %>%
  ggplot(aes(Rand_date, fit)) +
  facet_wrap(.~ N_dose, ncol = 3) +
  geom_point(data = data_for_plot_slope  %>% filter(Site=='th001', N_dose %in% 0:5), aes(x = Rand_date, y = slope_med, col = Trt),
             size = 3, alpha = 0.5) +
  scale_color_manual(values = col) +
  geom_smooth_ci(linewidth = 1, col = "#D80032", ci_alpha = 0.25) +
  theme_bw() +
  ylab("Viral reduction per day (log10 units)") +
  xlab("Randomisation date") +
  theme(axis.title = element_text(size = 10, face = "bold"),
        axis.text = element_text(size = 8)) +
  ylim(-3,0) +
  geom_hline(yintercept = 0, linetype = "dashed")
