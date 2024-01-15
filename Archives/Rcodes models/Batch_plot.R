#####################################################################
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
#####################################################################
library(stringr)
library(ggplot2)
library(ggforce)
#####################################################################
Mods <- c("Linear", "TwoCmpt", "UpDown", "Biexp")
Arms <- c("NR", "NS")

Pairs <- as.data.frame(expand.grid(Arms, Mods))
colnames(Pairs) <- c("Arms", "Mods")

fileList <- list.files("Fit")

i <- 2

dataPLOT <- NULL
for(i in 1:8){
  modID_check <- str_detect(fileList, paste0("0", as.character(rownames(Pairs)[i])))
  Arm <- Pairs$Arms[i]
  Mod <- Pairs$Mods[i]
  if(Arm == "NR"){ID_inc <- 1:30;Arm_lab <- "Nirmatrelvir + Ritonavir"} else {ID_inc <- c(1:24, 26:31);Arm_lab <- "No study drug"}


  load(paste0("Fit/", fileList[modID_check]))

  data <- read.csv("Analysis_Data/drug_test_NS_NR_F.csv")
  data <- data[data$Trt == Arm_lab,]
  data$ID_code <- data$ID
  data$ID <- as.numeric(as.factor(data$ID))
  data <- data[data$ID %in% ID_inc ,]
  data$ID <- as.numeric(as.factor(data$ID))
  data <- data[order(c(data$censor), decreasing = T),]

  samples <- fit@sim$samples[[4]]
  preds_ind <- which(str_detect(string = names(samples), pattern = "preds"))
  data$preds <- as.vector(sapply(samples[preds_ind], median))

  data_plot <- data.frame("ID_code" = data$ID_code, 
                          "log10_viral_load" = data$log10_viral_load, 
                          "Timepoint_ID" = data$Timepoint_ID,
                          "preds" = data$preds, 
                          "Trt" = data$Trt, 
                          "Mod" = Mod)
  dataPLOT <- rbind(dataPLOT, data_plot)
}


ggplot(data = dataPLOT[dataPLOT$Trt == "Nirmatrelvir + Ritonavir",]) +
  geom_point(mapping = aes(x = Timepoint_ID, y = log10_viral_load)) +
  geom_point(mapping = aes(x = Timepoint_ID, y = preds, col = Mod), size = 2) +
  geom_line(mapping = aes(x = Timepoint_ID, y = preds, col = Mod), linewidth = 0.5) +
  theme_bw() +
  facet_wrap_paginate(~ ID_code, ncol = 6, nrow = 5, page = 1) +
  scale_y_continuous(breaks = seq(-9,9,3), limits = c(-9,9)) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Time (days)") +
  ylab("Log 10 of viral loads") +
  ggtitle("Nirmatrelvir + Ritonavir") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = c("#FFB84C", "#F266AB", "#A459D1", "#2CD3E1"), name = "")



ggplot(data = dataPLOT[dataPLOT$Trt != "Nirmatrelvir + Ritonavir",]) +
  geom_point(mapping = aes(x = Timepoint_ID, y = log10_viral_load)) +
  geom_point(mapping = aes(x = Timepoint_ID, y = preds, col = Mod), size = 2) +
  geom_line(mapping = aes(x = Timepoint_ID, y = preds, col = Mod), linewidth = 0.5) +
  theme_bw() +
  facet_wrap_paginate(~ ID_code, ncol = 6, nrow = 5, page = 1) +
  scale_y_continuous(breaks = seq(-9,9,3), limits = c(-9,9)) +
  geom_hline(yintercept = 0, col = "red", linetype = "dashed") +
  xlab("Time (days)") +
  ylab("Log 10 of viral loads") +
  ggtitle("No study drug") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        title = element_text(size = 14, face = "bold")) +
  scale_color_manual(values = c("#FFB84C", "#F266AB", "#A459D1", "#2CD3E1"), name = "")
  