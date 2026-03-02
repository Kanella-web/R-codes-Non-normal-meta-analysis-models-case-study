## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

######### Frequentist models ###########
library(metafor)
####### normal(a) model ##########
set.seed(1508)

freq_norm1 = rma(measure = "SMD", m1i = pre_term_data_58$mean_EPT.VPT, m2i =  pre_term_data_58$mean_FT,
                 sd1i = pre_term_data_58$sd_EPT.VPT, sd2i =  pre_term_data_58$sd_FT,
                 n1i = pre_term_data_58$n_EPT.VPT, n2i = pre_term_data_58$n_FT, slab = pre_term_data_58$Study,
                 method = "REML")  # DL
forest(freq_norm1, cex = 0.4)
text(c(-6.3, 2.1), 68, c("Studies" , "Estimate[95% CI]"),   font=2, cex=1)

freq_norm <- summary(freq_norm1)

fr_norm = as.data.frame(confint(freq_norm1))

stud_eff = blup( freq_norm, level = 95 )

mu_fr_norm <- freq_norm$beta
LB_mu_fr_norm <- freq_norm$ci.lb
UB_mu_fr_norm <- freq_norm$ci.ub
prec_mu_fr_norm <- (UB_mu_fr_norm - LB_mu_fr_norm) 


fr_tau2 <- grepl("tau\\^2", row.names(fr_norm))
fr_tau2_n <- fr_norm[fr_tau2,]
tau2_fr_norm <- fr_tau2_n[1]
LB_tau2_fr_norm <- fr_tau2_n[2]
UB_tau2_fr_norm <- fr_tau2_n[3]
prec_tau2_fr_norm <- (UB_tau2_fr_norm - LB_tau2_fr_norm) 

pred = as.data.frame(predict(freq_norm))
theta_new = pred$pred
LB_theta_new = pred$pi.lb
UB_theta_new = pred$pi.ub
prec_theta_new_fr_norm = UB_theta_new - LB_theta_new

fr_normal <- cbind(mu_fr_norm,LB_mu_fr_norm, UB_mu_fr_norm, 
                   tau2_fr_norm, LB_tau2_fr_norm, UB_tau2_fr_norm, 
                   prec_mu_fr_norm, prec_tau2_fr_norm,
                   theta_new, LB_theta_new, UB_theta_new, prec_theta_new_fr_norm)

# write.csv(fr_normal, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\normal_model_58\\N_res_58SMDs.csv",row.names=FALSE )  
# write.csv(stud_eff, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\normal_model_58\\N_rel_eff_58SMDs.csv",row.names=FALSE )  

################ t-distribution(a) model ##############
library(metaplus)
set.seed(1508)

data_SMDs = escalc(measure = "SMD", m1i = pre_term_data_58$mean_EPT.VPT, m2i =  pre_term_data_58$mean_FT,
                   sd1i = pre_term_data_58$sd_EPT.VPT, sd2i =  pre_term_data_58$sd_FT,
                   n1i = pre_term_data_58$n_EPT.VPT, n2i = pre_term_data_58$n_FT, slab = pre_term_data_58$Study)

metaplust1 <- metaplus(yi =  data_SMDs$yi , sei = sqrt(data_SMDs$vi) , random = "t-dist")
metaplust <- summary(metaplust1)

forest(metaplust1, cex = 0.5 )

mu_metaplust <- metaplust$results[1]
LB_mu_metaplust <- metaplust$results[4]
UB_mu_metaplust <- metaplust$results[7]
prec_mu_metaplust <- UB_mu_metaplust - LB_mu_metaplust

tau2_metaplust <- metaplust$results[2]
####  if vinv = 0 => infinite df => normal distribution 
ind_normal <- metaplust$results[3]  

##### outliers test 
out_test = testOutliers(metaplust1)
out_test = out_test$pvalue


metaplus_t <- cbind(mu_metaplust ,LB_mu_metaplust, UB_mu_metaplust, tau2_metaplust, 
                    prec_mu_metaplust,out_test , ind_normal)

#write.csv(metaplus_t, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\t_model_58\\t_res_58SMDs.csv",row.names=FALSE )  

########## Mixture of two normal distributions with common-component mean(a) model #############
########## Using metaplus mixture to test for outliers #####
set.seed(1508)
library(metaplus)

data_SMDs = escalc(measure = "SMD", m1i = pre_term_data_58$mean_EPT.VPT, m2i =  pre_term_data_58$mean_FT,
                   sd1i = pre_term_data_58$sd_EPT.VPT, sd2i =  pre_term_data_58$sd_FT,
                   n1i = pre_term_data_58$n_EPT.VPT, n2i = pre_term_data_58$n_FT, slab = pre_term_data_58$Study)

metaplusmix1 <- metaplus(yi =  data_SMDs$yi , sei = sqrt(data_SMDs$vi) , random = "mixture")

forest(metaplusmix1, cex = 0.5)

metaplusmix <- summary(metaplusmix1)

#### PROBABILITY OF EACH STUDY BEING AN OUTLIER #####
set.seed(1508)
outlierProbs <- outlierProbs(metaplusmix1 )


mu_metaplusmix <- metaplusmix$results[1]
LB_mu_metaplusmix <- metaplusmix$results[5]
UB_mu_metaplusmix <- metaplusmix$results[9]
prec_mu_metaplusmix <- UB_mu_metaplusmix - LB_mu_metaplusmix 

tau2_metaplusmix <- metaplusmix$results[2]

tau2_out_metaplusmix <- metaplusmix$results[3]

prob_out <- metaplusmix$results[4]

out_test_mix1 = summary(testOutliers(metaplusmix1))
out_test_mix = out_test_mix1$pvalue

metaplus_mix <- cbind(mu_metaplusmix ,LB_mu_metaplusmix, UB_mu_metaplusmix, tau2_metaplusmix, tau2_out_metaplusmix,
                      prec_mu_metaplusmix, prob_out,out_test_mix  )
# 
# write.csv(metaplus_mix, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\common_mean_mix_model_58\\mix_res_58SMDs.csv",row.names=FALSE )  
# write.csv(outlierProbs$outlier.prob, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\common_mean_mix_model_58\\outlier_prob_res_58SMDs.csv",row.names=FALSE )  


outlier.prob <- rep(outlierProbs$outlier.prob, length.out = 58)

prob_oyt_study = cbind.data.frame(pre_term_data_58$Study, outlier.prob)
prob_out = round(prob_oyt_study$outlier.prob,3)
prob_out_study = cbind.data.frame(prob_oyt_study$`pre_term_data_58$Study`, prob_out)

#write.csv(prob_out_study, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\common_mean_mix_model_58\\outlier_prob_per_study_58SMDs.csv",row.names=FALSE )  


