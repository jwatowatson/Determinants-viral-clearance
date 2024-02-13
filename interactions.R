traceplot(out, "trt_effect")
plot(out, pars = "trt_effect")

#str(out)
trt_effect <- rstan::extract(out, "trt_effect")$trt_effect
Trt_matrix <- stan_inputs[[model_settings$dataset[ind_res]]]$Trt_matrix
effect_names <- colnames(Trt_matrix)[-1]
#ind_nodata <- (apply(Trt_matrix[,-1], 2, sum) == 0)

ind_nodata <- as.data.frame(table(Baseline_data$VariantClass_new, Baseline_data$Trt))
colnames(ind_nodata) <- c("Variant", "Treatment", "Freq")
ind_nodata$Available <- ind_nodata$Freq > 0
ind_nodata$Treatment <- as.factor(ind_nodata$Treatment)
levels(ind_nodata$Treatment) <- gsub(' [A-z+]*',"",levels(ind_nodata$Treatment))
ind_nodata$Variant <- as.factor(ind_nodata$Variant)
ind_nodata$lab <- paste0("n = ", ind_nodata$Freq)

colnames(trt_effect) <- effect_names
ref_variant <- "BA.1"
Trt_list <- which(grepl("Trt", effect_names) & !grepl("VariantClass_new", effect_names))
Trts <- sub("Trt", "", effect_names[Trt_list])
Trts[grep("Nirmatrelvir", Trts)] <- "Nirmatrelvir"

trt_effect_summarize <- NULL
#nodata <- NULL
for(i in 1:length(Trts)){
  subdat <-  trt_effect[,grep(Trts[i], effect_names)]
  sum_dat <- subdat[,grepl("VariantClass_new", colnames(subdat))] +  subdat[,!grepl("VariantClass_new", colnames(subdat))]
  sum_dat <- cbind(subdat[,!grepl("VariantClass_new", colnames(subdat))], sum_dat)
  colnames(sum_dat)[1] <- paste0("Trt", Trts[i],":VariantClass_new",ref_variant)
  trt_effect_summarize <- cbind(trt_effect_summarize, sum_dat)
  
 # nodata <- c(nodata, ind_nodata[grep(Trts[i], names(ind_nodata))])
  
}
trt_effect_summarize

effect_names2 <- colnames(trt_effect_summarize)

trt_effect_summarize <- as.data.frame(t(apply(trt_effect_summarize,2,quantile, c(0.025, 0.5, 0.975))))
colnames(trt_effect_summarize) <- c("Low", "Med", "Up")
trt_effect_summarize$predictors <- effect_names2
rownames(trt_effect_summarize) <- NULL

trt_effect_summarize_int <- trt_effect_summarize[grep(":", trt_effect_summarize$predictors),]
trt_effect_summarize_int$Trt <-  sapply(strsplit(trt_effect_summarize_int$predictors, ":"), "[", 1)
trt_effect_summarize_int$Variant <-  sapply(strsplit(trt_effect_summarize_int$predictors, ":"), "[", 2)

trt_effect_summarize_int$Trt <- gsub("Trt", "", trt_effect_summarize_int$Trt)
trt_effect_summarize_int$Variant <- gsub("VariantClass_new", "", trt_effect_summarize_int$Variant)

trt_effect_summarize_int$Trt <- as.factor(trt_effect_summarize_int$Trt)
levels(trt_effect_summarize_int$Trt) <- gsub(' [A-z+]*',"",levels(trt_effect_summarize_int$Trt))

trt_effect_summarize_int <- merge(trt_effect_summarize_int, ind_nodata, by.x = c("Trt", "Variant"), by.y = c("Treatment", "Variant"))
levels(trt_effect_summarize_int$Trt)[5] <- "Casirivimab/imdevimab"
levels(trt_effect_summarize_int$Trt)[4] <- "Ritonavir-boosted nirmatrelvir"
 
trt_effect_summarize_int$Variant <- as.factor(trt_effect_summarize_int$Variant)
trt_effect_summarize_int$Variant <- factor(trt_effect_summarize_int$Variant, 
                                                 levels = rev(c("Delta", "BA.1",  "BA.2" ,"BA.2.75" , "BA.4"  , 
                                                          "BA.5"   , "XBB"  , "XBB.1.5-like",  "Other" )))
trt_effect_summarize_int$Med_percent <- formatter(exp(trt_effect_summarize_int$Med))
trt_effect_summarize_int$Low_percent <- formatter(exp(trt_effect_summarize_int$Low))
trt_effect_summarize_int$Up_percent <- formatter(exp(trt_effect_summarize_int$Up))

#trt_effect_summarize_int$Nodata <- nodata


G_variants <- ggplot(trt_effect_summarize_int[trt_effect_summarize_int$Available,],
        aes(y = Variant)) +
  geom_point(aes(x = Med_percent), size = 2.5, alpha = 0.75, col = "#1D2B53") +
  geom_errorbar(aes(xmin = Low_percent, xmax = Up_percent), width = 0, linewidth = 1, alpha = 0.75, col = "#1D2B53") +
  facet_wrap(.~ Trt) +
  theme_bw() +
  geom_vline(xintercept = 0, linetype = "dashed", col = "red") +
  xlab("Change in viral clearance rate (%)") +
  ylab("") +
  theme(strip.text = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        panel.spacing = unit(1, "lines")) +
  geom_text(aes(y = Variant, label = lab), x = -200, hjust  = 0, size = 3) +
  coord_cartesian(xlim=c(-200, 500))
G_variants

png("Plots/variant_effects.png", width = 10, height = 6, units = "in", res = 350)
G_variants
dev.off()
