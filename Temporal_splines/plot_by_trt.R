library(ggpubr)
library(grid) 


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

for(i in 1:length(Trts)){
  
  trt <- Trts[i]
  
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
    ggtitle(trt) +
    ylim(-3,0)
  
  plot_list[[i]] <- G
  
}


GG <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2)


GG2 <- annotate_figure(GG, left = textGrob("Viral reduction per day (log10 units)", rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
                       bottom = textGrob("Randomisation date", gp = gpar(cex = 1.3)))

png("slope_by_time_trt.png", width = 10, height = 6, units = "in", res = 350)
GG2
dev.off()




