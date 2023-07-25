library(rstan)
###########################################################################
mainDir <- "D:/Determinants-viral-clearance"
setwd(mainDir)
###########################################################################
data <- read.csv("Analysis_Data/interim_all_analysis.csv")
data$Rand_date <- as.Date(data$Rand_date)
data <- data[data$Site == "th001",]
data <- data[data$Timepoint_ID %in% c(0:7, 14),]

data$Sex <- as.factor(data$Sex)
levels(data$Sex) <- c("Female", "Male")

data$Epoch <- as.factor(data$Epoch)
levels(data$Epoch) <- c(paste0("Epoch ", 0:5))
###########################################################################
control_arm <- "No study drug"
positive_arm <- "Nirmatrelvir + Ritonavir"
test_arm <- "Fluoxetine"
###########################################################################
# Define test period
rand_date_data <- aggregate(Rand_date  ~ Trt, data = data, FUN = function(x) c(min = min(x), max = max(x)))
rand_date_data[,2:3] <- data.frame(rand_date_data[,2][,1], rand_date_data[,2][,2])
colnames(rand_date_data) <- c("Trt", "min_date", "max_date") 
rand_date_data$min_date <- as.Date(rand_date_data$min_date)
rand_date_data$max_date <- as.Date(rand_date_data$max_date)

control_duration <- seq.Date(rand_date_data$min_date[rand_date_data$Trt == control_arm], 
                             rand_date_data$max_date[rand_date_data$Trt == control_arm],
                             "1 day")
positive_duration <- seq.Date(rand_date_data$min_date[rand_date_data$Trt == positive_arm], 
                             rand_date_data$max_date[rand_date_data$Trt == positive_arm],
                             "1 day")
test_duration <- seq.Date(rand_date_data$min_date[rand_date_data$Trt == test_arm], 
                              rand_date_data$max_date[rand_date_data$Trt == test_arm],
                              "1 day")
common_duration <- as.Date(Reduce(intersect, list(control_duration,positive_duration,test_duration)))
min_study_date <- min(common_duration)
max_study_date <- max(common_duration)
###########################################################################
# Subsetting data
subdata <- data[data$Trt %in% c(control_arm, positive_arm, test_arm) & data$Rand_date %in% common_duration,]

subdata$Weight_class <- NA
subdata$Weight_class[subdata$Weight >= 60] <- "Weight >= 60 kg"
subdata$Weight_class[subdata$Weight < 60] <- "Weight < 60 kg"
subdata$Weight_class <- as.factor(subdata$Weight_class)

# visualize
ggplot(data = subdata, aes(x=Rand_date, y = Trt, col = Variant)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab("") +
  ylab("") +
  geom_vline(xintercept = seq.Date(as.Date("2022-01-01"), as.Date("2024-01-01"), "6 months"),
             col = "red", linetype = "dashed") +
  scale_x_date(date_labels =  "%b %y") +
  theme(axis.text.x = element_text(angle = -90, vjust = 0.5, hjust=0),
        strip.text = element_text(size = 12, face = "bold")) +
  scale_color_brewer(palette="Set1")

# number by group
table(unique(data.frame(subdata$ID, subdata$Trt))[,2])
###########################################################################
subdata_median <-  aggregate(log10_viral_load~Timepoint_ID+Trt, 
                             data = subdata, FUN = median)

G1 <- ggplot() +
  geom_jitter(data = subdata, aes(x = Timepoint_ID, y = log10_viral_load, 
                                  shape = censor, col = Trt), width = 0.2, alpha = 0.4, size = 1) +
  scale_shape_manual(values = c(2,1), name = "Censor") +
  geom_line(data = subset(subdata_median, as.numeric(Timepoint_ID) <= 7), 
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 1, linetype = 1 ) + 
  geom_line(data =  subset(subdata_median, as.numeric(Timepoint_ID) >= 7),
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 0.75, linetype = 2 ) + 
  geom_point(data = subdata_median, aes(x = Timepoint_ID, y = log10_viral_load, col = Trt), size = 2.5) + 
  scale_color_manual(values = c("#00235B", "#E21818", "#EBB02D"), name = "Treatment") +
  theme_bw() +
  scale_y_continuous(breaks = 0:10) +
  scale_x_continuous(breaks = 0:14) +
  xlab("Time (days)") +
  ylab("Log10 of viral loads") +
  theme(axis.title = element_text(size = 12, face = "bold"))

pdf("Results/Observed_vl_NS_NR_F.pdf", width = 8, height = 5)
G1
dev.off()
###########################################################################
subdata2 <- subdata[!(is.na(subdata$Weight_class)),]
subdata_median_Sex_weight <-  aggregate(log10_viral_load~Timepoint_ID+Trt+Weight_class, 
                             data = subdata2, FUN = median)

G2 <- ggplot() +
  facet_grid(.~Weight_class)+
  geom_jitter(data = subdata2, aes(x = Timepoint_ID, y = log10_viral_load, 
                                  shape = censor, col = Trt), width = 0.2, alpha = 0.4, size = 1) +
  scale_shape_manual(values = c(2,1), name = "Censor") +
  geom_line(data = subset(subdata_median_Sex_weight, as.numeric(Timepoint_ID) <= 7), 
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 1, linetype = 1 ) + 
  geom_line(data =  subset(subdata_median_Sex_weight, as.numeric(Timepoint_ID) >= 7),
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 0.75, linetype = 2 ) + 
  geom_point(data = subdata_median_Sex_weight, aes(x = Timepoint_ID, y = log10_viral_load, col = Trt), size = 2.5) + 
  scale_color_manual(values = c("#00235B", "#E21818", "#EBB02D"), name = "Treatment") +
  theme_bw() +
  scale_y_continuous(breaks = 0:10) +
  scale_x_continuous(breaks = 0:14) +
  xlab("Time (days)") +
  ylab("Log10 of viral loads") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"))

pdf("Results/Observed_vl_NS_NR_F_by_weight.pdf", width = 10, height = 5)
G2
dev.off()
###########################################################################
subdata_median_Epoch <-  aggregate(log10_viral_load~Timepoint_ID+Trt+Epoch, 
                                        data = subdata, FUN = median)

G3 <- ggplot() +
  facet_wrap(.~Epoch, ncol = 2, nrow = 3)+
  geom_jitter(data = subdata2, aes(x = Timepoint_ID, y = log10_viral_load, 
                                   shape = censor, col = Trt), width = 0.2, alpha = 0.4, size = 1) +
  scale_shape_manual(values = c(2,1), name = "Censor") +
  geom_line(data = subset(subdata_median_Epoch, as.numeric(Timepoint_ID) <= 7), 
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 1, linetype = 1 ) + 
  geom_line(data =  subset(subdata_median_Epoch, as.numeric(Timepoint_ID) >= 7),
            aes(x = Timepoint_ID, y = log10_viral_load, col = Trt, group = Trt), linewidth = 0.75, linetype = 2 ) + 
  geom_point(data = subdata_median_Epoch, aes(x = Timepoint_ID, y = log10_viral_load, col = Trt), size = 2.5) + 
  scale_color_manual(values = c("#00235B", "#E21818", "#EBB02D"), name = "Treatment") +
  theme_bw() +
  scale_y_continuous(breaks = 0:10) +
  scale_x_continuous(breaks = 0:14) +
  xlab("Time (days)") +
  ylab("Log10 of viral loads") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(size = 10, face = "bold"))

pdf("Results/Observed_vl_NS_NR_F_by_epoch.pdf", width = 10, height = 6)
G3
dev.off()
###########################################################################
write.csv(subdata, "Analysis_Data/drug_test_NS_NR_F.csv", row.names = F)
###########################################################################