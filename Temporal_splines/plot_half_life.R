Baseline_vl_data_ALL

Baseline_vl_avg <- Baseline_vl_data_ALL %>%
  group_by(ID) %>%
  summarise_at(vars(log10_viral_load), list(log10_viral_load = mean))

data_for_plot_slope2 <- merge(data_for_plot_slope, Baseline_vl_avg, by = "ID")

plot(data_for_plot_slope2$log10_viral_load.y, data_for_plot_slope$slope_med)

library(ggplot2)

ggplot(data_for_plot_slope2, aes(x = log10_viral_load.y, y = hl_med)) +
  geom_point(size = 3, alpha = 0.5, aes(col = Study_time_normal)) +
  theme_bw() +
  facet_wrap(Trt~.) + 
  scale_colour_continuous(type = "viridis", name = "Study time") +
  xlab("Baseline viral loads") +
  ylab("Viral clearance half-life (hours)") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines"))




data_for_plot_slope2 <- transform(data_for_plot_slope2, vl_group = cut(log10_viral_load.y, 
                                               c(seq(1, 9, 1))))
data_for_plot_slope2$vl_group <- as.character(data_for_plot_slope2$vl_group)
data_for_plot_slope2$vl_group[data_for_plot_slope2$vl_group %in% c("(1,2]")] <- "(2,3]"
data_for_plot_slope2$vl_group[data_for_plot_slope2$vl_group %in% c("(8,9]")] <- "(7,8]"
data_for_plot_slope2$vl_group <- as.factor(data_for_plot_slope2$vl_group)
levels(data_for_plot_slope2$vl_group) <- c("<3", "3-4", "4-5",
                                           "5-6", "6-7", ">7")

data_for_plot_slope3 <- data_for_plot_slope2[data_for_plot_slope2$Trt == "No study drug",]

quantile <- aggregate(list("hl_med" = data_for_plot_slope3$hl_med), by = list("vl_group" = data_for_plot_slope3$vl_group), quantile, c(0.25, 0.5, 0.75))
quantile$Q1 <- quantile$hl_med[,1]
quantile$Q2 <- quantile$hl_med[,2]
quantile$Q3 <- quantile$hl_med[,3]

ggplot() +
  geom_jitter(data = data_for_plot_slope3, aes(x = vl_group, y = hl_med), size = 3, alpha = 0.1, width = 0.2) +
  theme_bw() +
  geom_point(data = quantile, aes(x = vl_group, y = Q2), col = "#CD1818", size = 3.5) +
  geom_errorbar(data = quantile, aes(x = vl_group, ymin = Q1, ymax = Q3), width = 0, col = "#CD1818",
                linewidth = 1.2) +
  xlab("Baseline viral loads (log10 units)") +
  ylab("Viral clearance half-life (hours)") +
  ggtitle("No study drug") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines"))
  
  
quantile$vl_group


?geom_errorbar

quantile$Q1
quantile$vl_group
quantile$Q3
