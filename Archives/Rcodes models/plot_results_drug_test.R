library(purrr)
library(stringr)
library(ggplot2)
library(ggforce)
library(rstan)
library(dplyr)

getwd()
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)

load(file = "Fit/TwoCmpt_nRNaseP_NS_only_fixed_eta2.RData")
data <- read.csv("Analysis_Data/drug_test_NS_NR_F.csv")
# choose drug
data <- data[data$Trt == "No study drug",] #  "Nirmatrelvir + Ritonavir", ]
data$ID_code <- data$ID
data$ID <- as.numeric(as.factor(data$ID))
# choose the first 30
data <- data[data$ID %in% c(1:24, 26:31)  ,] #,] #1:30
data <- data[order(c(data$censor), decreasing = T),]

samples <- fit@sim$samples[[4]]
preds_ind <- which(str_detect(string = names(samples), pattern = "preds"))
data$preds <- as.vector(sapply(samples[preds_ind], median))

#data$preds2 <- data$preds * data$CT_RNaseP

#par(mfrow = c(1,2))
#hist(samples[preds_ind][[1]], main = "Posterior of gamma_rnasep", xlab = "gamma_rnasep", xlim = c(-0.3,0.3))
#abline(v = median(samples[preds_ind][[1]]), col = "red")

#hist(data$preds2, main = "Predicted RNaseP effects", xlab = "gamma_rnasep * CT_RNaseP", xlim = c(-3.5,-2))
#abline(v = median(data$preds2), col = "red")

#traceplot(fit, pars = c("gamma_rnasep"))

G <- ggplot(data = data) +
  geom_point(mapping = aes(x = Timepoint_ID, y = log10_viral_load)) +
  geom_point(mapping = aes(x = Timepoint_ID, y = preds), col = "red", size = 2) +
  geom_line(mapping = aes(x = Timepoint_ID, y = preds), col = "red") +
  theme_bw() +
  facet_wrap_paginate(~ ID_code, ncol = 6, nrow = 5, page = 1) +
  scale_y_continuous(breaks = seq(-4,9,2), limits = c(-6,9))
G

estim_lambda2_ind <- which(str_detect(string = names(samples), pattern = "estim_lambda2"))
lambda2 <- (lapply(samples[estim_lambda2_ind],  quantile, c(0.025, 0.5, 0.975)))
lambda2 <- as.data.frame(bind_rows(lambda2, .id = "column_label"))
colnames(lambda2) <- c("ID", "Low", "Med", "Up")
lambda2$ID_code <- as.factor(unique(data$ID_code))
lambda2$half_life <- log10(2)/lambda2$Med
lambda2$lab <- paste0("Rate = ", sprintf("%.4f", round(lambda2$Med, 4)))
lambda2$lab2 <- paste0("Half-life = ", sprintf("%.2f", round(lambda2$half_life, 2)), " days")

G2 <- ggplot(data = data) +
  geom_point(mapping = aes(x = Timepoint_ID, y = log10_viral_load)) +
  geom_point(mapping = aes(x = Timepoint_ID, y = preds), col = "red", size = 2) +
  geom_line(mapping = aes(x = Timepoint_ID, y = preds), col = "red") +
  theme_bw() +
  facet_wrap_paginate(~ ID_code, ncol = 6, nrow = 5, page = 1) +
  scale_y_continuous(breaks = seq(-4,9,2), limits = c(-4,9)) +
  xlab("Time (days)") +
  ylab("Log10 of viral loads") +
  geom_text(data = lambda2, aes(x = 0, y = -1, label = lab), size = 3, hjust = 0) +
  geom_text(data = lambda2, aes(x = 0, y = -3, label = lab2), size = 3, hjust = 0) +
  theme(axis.title = element_text(size = 14, face = "bold"))

pdf("Results/TwoCmpt_NR_only.pdf", width = 12, height = 8)
G2
dev.off()



