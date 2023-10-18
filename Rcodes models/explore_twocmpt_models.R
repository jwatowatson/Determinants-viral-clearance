library(rstan)
library(stringr)
#####################################################################
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
#############################################################
Mods <- c("Linear", "TwoCmpt", "UpDown", "Biexp")
Arms <- c("NR", "NS")
#############################################################
Pairs <- as.data.frame(expand.grid(Arms, Mods))
colnames(Pairs) <- c("Arms", "Mods")

fileList <- list.files("Fit")
i <- 3

modID_check <- str_detect(fileList, paste0("0", as.character(rownames(Pairs)[i])))
(Arm <- Pairs$Arms[i])
(Mod <- Pairs$Mods[i])

load(paste0("Fit/", fileList[modID_check]))
#############################################################
vl_cal_m4 <- function(logA0, logB0, lambda1, lambda2){
  t <- seq(0, 14, 0.1)
  vl_t <- log10(  ((lambda1/(lambda2-lambda1)) * exp(logA0) * (exp(-lambda1 * t) - exp(-lambda2 * t))) + (exp(logB0) * (exp(-lambda2 * t))))
  out <- vl_t
  return(out)
}
#############################################################
sample <- fit@sim$samples[[1]]
ii <- 1
theta1 <- sample[names(sample) == paste0("theta_rand_id[", ii, ",1]")][[1]]
theta2 <- sample[names(sample) == paste0("theta_rand_id[", ii, ",2]")][[1]]
theta3 <- sample[names(sample) == paste0("theta_rand_id[", ii, ",3]")][[1]]
theta4 <- sample[names(sample) == paste0("theta_rand_id[", ii, ",4]")][[1]]

lambda1 <- exp(sample$loglambda1_0 + theta1)
lambda2 <- exp(sample$loglambda2_0 + theta2)
logA0 <- sample$logA0_0 + theta3
logB0 <- sample$logB0_0 + theta4

par(mfrow = c(2,5))

ALL <- NULL
for(x1 in c(0.5,1,2)){
  for(x2 in c(0.5,1,2)){
    
res <- mapply(vl_cal_m4, logA0, logB0, lambda1*x1, lambda2*x2)
res2 <- apply(res,1,median)
up <- apply(res,1,quantile, 0.975)
low <- apply(res,1,quantile, 0.025)

lam1 <- median(lambda1) * x1
lam2 <- median(lambda2) * x2

out <- data.frame(x1, x2, res2, up, low, lam1, lam2)
ALL <- rbind(ALL, out)
#par(mfrow = c(2,3))
#plot(t, up, type = "l", lwd = 1, ylim=c(-5,8), lty = 2, col = "red",
#     ylab = "Log 10 of viral loads", 
#     xlab = "Time (days)",
#     main = paste0(xxx,"x lambda2"))
#lines(t, res2, lwd = 3)
#lines(t, low, lwd = 1, lty = 2, col = "red")
#abline(h = 0, lty = 3)
#hist(log10(lambda1))
#hist(log10(lambda2))
#hist(log10(exp(logA0)))
#hist(log10(exp(logB0)))

}
}


ALL
t <- seq(0, 14, 0.1)
ALL$time <- t
library(ggplot2)
ALL$x1 <- as.factor(ALL$x1)
ALL$x2 <- as.factor(ALL$x2)

levels(ALL$x1) <- c("0.5x lambda1", "1.0x lambda1", "2.0x lambda1")
levels(ALL$x2) <- c("0.5x lambda2", "1.0x lambda2", "2.0x lambda2")

ALL$lam1 <- paste0("lambda1 = ", round(ALL$lam1,2))
ALL$lam2 <- paste0("lambda2 = ", round(ALL$lam2,2))

lab <- unique(ALL[,c("x1", "x2", "lam1", "lam2")])


ggplot(ALL, aes(x = time)) +
  geom_ribbon(aes(ymin=low, ymax=up), alpha = 0.15, fill = "#009292") +
  geom_line(aes(y = res2), linewidth = 1) +
  facet_grid(x2~x1) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) + 
  coord_cartesian(ylim = c(-7, 7), xlim = c(0,14)) +
  scale_x_continuous(breaks = seq(0,14,2)) +
  xlab("Time (days)") +
  ylab("Log 10 of viral loads") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold")) +
  geom_text(data = lab, aes(x = 7, y = 6, label = lam1), size = 3, hjust = 0) +
  geom_text(data = lab, aes(x = 7, y = 4, label = lam2), size = 3, hjust = 0) 








