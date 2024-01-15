data_for_plot_slope_NSD <-  data_for_plot_slope[data_for_plot_slope$Trt == "No study drug", ]

A <- ggplot() +
  theme_bw() +
  geom_point(data = data_for_plot_slope_NSD, aes(x = Time_since_last_dose, y = slope_med, color = as.Date(Rand_date)),
              size = 3, alpha = 0.75) +
  scale_colour_viridis_c(trans = "date", name = "Randomisation date") +
  ylim(-2.5, 0) +
  ylab("Viral reduction per day (log10 units)") +
  xlab("Time since the last vaccination (days)") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed")
  
  
cor.test(data_for_plot_slope_NSD$Time_since_last_dose, data_for_plot_slope_NSD$slope_med, method = "spearman")



# B <- ggplot() +
#   theme_bw() +
#   geom_point(data = data_for_plot_slope_NSD, aes(x = Time_since_last_dose, y = slope_med, color = N_dose),
#              size = 3, alpha = 0.75) +
#   scale_colour_viridis_b(name = "Number of vaccinations \n(doses)", option = "plasma") +
#   ylim(-2.5, 0) +
#   ylab("Viral reduction per day (log10 units)") +
#   xlab("Time since the last vaccination (days)") +
#   theme(strip.text = element_text(size = 10, face = "bold"),
#         axis.title = element_text(size = 12, face = "bold"),
#         plot.title = element_text(size = 12, face = "bold"),
#         panel.spacing = unit(1, "lines")) +
#   geom_hline(yintercept = 0, linetype = "dashed")

##################################################################################################
B <- ggplot() +
  theme_bw() +
  geom_point(data = data_for_plot_slope_NSD, aes(x = Rand_date, y = Time_since_last_dose, col = N_dose),
             size = 3, alpha = 0.75) +
  scale_colour_viridis_b(name = "Number of vaccinations \n(doses)", option = "plasma") +
  xlab("Randomisation date") +
  ylab("Time since the last vaccination (days)") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        panel.spacing = unit(1, "lines"),
        legend.position = "right") +
  geom_hline(yintercept = 0, linetype = "dashed") 


png("Plots/Time_since_last_vaccines.png", width = 12, height = 5, units = "in", res = 350)
ggarrange(B, A, labels = "AUTO")
dev.off()
