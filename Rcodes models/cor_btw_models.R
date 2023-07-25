library(rstan)
library(stringr)
########################################################################
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
########################################################################
load("Fit/Linear_naive_nRNaseP.RData")
linear_fit <- fit
########################################################################
load("Fit/TwoCmpt_nRNaseP.RData")
twocmpt_fit <- fit
########################################################################

posterior_twocmpt <- twocmpt_fit@sim$samples[[1]]
lambda2_ind <- which(str_detect(names(posterior_twocmpt), "lambda2") & !str_detect(names(posterior_twocmpt), "log"))
lambda2 <- NULL
for(i in lambda2_ind){
  lambda2 <- c(lambda2, (median(unlist(posterior_twocmpt[i]))))
}


posterior_linear <- linear_fit@sim$samples[[1]]
beta_ind <- which(str_detect(names(posterior_linear), "beta") & !str_detect(names(posterior_linear), "_0"))
beta <- NULL
for(i in beta_ind){
  beta <- c(beta, (-1* median(unlist(posterior_linear[i]))))
}
beta

dat <- data.frame(beta, lambda2)

cor.test(lambda2, beta)

library(ggplot2)
ggplot(dat, aes(y = beta, x = lambda2)) + geom_point(size = 4) +
  theme_bw() +
  xlab("Viral clearance rate from two-compartment model") +
  ylab("Viral vlearance rate from log-linear decay model") +
  theme(axis.title = element_text(size = 14, face = "bold"))
