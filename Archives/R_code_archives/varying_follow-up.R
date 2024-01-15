interventions = c('Remdesivir', 'Molnupiravir')


ref_arm = 'No study drug'
study_threshold = 1.2 # depending on comparison with no study drug or positive control

f_name = paste0('Analysis_Data/',intervention,'_analysis.csv')
platcov_dat = read.csv(f_name)
platcov_dat$Trt[platcov_dat$Trt=='Nirmatrelvir + Ritonavir']='Nirmatrelvir'

platcov_dat = platcov_dat %>% filter(!Trt %in% c('Ensitrelvir','Nitazoxanide'))