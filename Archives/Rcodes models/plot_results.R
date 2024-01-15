library(purrr)
library(stringr)
library(ggplot2)
library(ggforce)
library(rstan)
getwd()
load(file = "Swab variations/Fit/TwoCmpt_nRNaseP.RData")
data <- read.csv("Swab variations/Data/swabber_analysis.csv")
data$ID_code <- data$ID
data$ID <- as.numeric(as.factor(data$ID))
data <- data[order(c(data$censor), decreasing = T),]

samples <- fit@sim$samples[[1]]
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
  scale_y_continuous(breaks = seq(-4,9,2), limits = c(-4,9))

pdf("Swab variations/Results/TwoCmpt_model_1.pdf", width = 12, height = 8)
G
dev.off()



