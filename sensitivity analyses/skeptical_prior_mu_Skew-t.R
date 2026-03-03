## download the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

library(metafor)
# ################################## Skew-t model #############################
## prepare for the stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
### donload the stan model
url <-"https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/sk_mu_normal_ST_U_SMDs.stan"

filesk_ST_U <- "sk_mu_normal_ST_U_SMDs.stan"

download.file(url, filesk_ST_U)

modsk_ST_U <- cmdstan_model(filesk_ST_U)

fitsk_ST_U <-  modsk_ST_U$sample(
  data =list(ns = nrow(pre_term_data_58),
             y1 = pre_term_data_58$mean_FT,
             y2 = pre_term_data_58$mean_EPT.VPT,
             sd1 = pre_term_data_58$sd_FT,
             sd2 = pre_term_data_58$sd_EPT.VPT,
             n1 = pre_term_data_58$n_FT,
             n2 = pre_term_data_58$n_EPT.VPT),
  seed = 1508, 
  chains = 2, 
  parallel_chains = 2,
  refresh = 500,
  iter_warmup = 10000,
  iter_sampling = 50000,
  adapt_delta = 0.99
)

res_sk_ST_U <- as.data.frame(fitsk_ST_U$summary())
View(res_sk_ST_U)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitsk_ST_U$draws(c("mu", "pred")))


library(bayesplot)
draws <- fitsk_ST_U$draws()  

### traceplots 
bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", "nu",
                                      "xi", "skew","kurt"))
delta_names <- fitsk_ST_U$variable[grepl("delta",fitsk_ST_U$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout

### Pr of mu <- 0
Pr_mu = mean(fitsk_ST_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fitsk_ST_U$draws("pred") < 0 )


d_sk_ST_U <- grepl("delta",res_sk_ST_U$variable)

deltask_ST_U <- res_sk_ST_U[d_sk_ST_U,]

plot(density(deltask_ST_U$median))

dm_sk_ST_U <- grepl("mu",res_sk_ST_U$variable)

median_sk_ST_U <- res_sk_ST_U[dm_sk_ST_U,]$median
LBmu_sk_ST_U <- res_sk_ST_U[dm_sk_ST_U,]$q5
UBmu_sk_ST_U <- res_sk_ST_U[dm_sk_ST_U,]$q95
Rhat_musk_ST_U <-  res_sk_ST_U[dm_sk_ST_U,]$rhat
prec_mu_sk_ST_U <- UBmu_sk_ST_U - LBmu_sk_ST_U

dtau2_sk_ST_U <- grepl("tau_sqr",res_sk_ST_U$variable)

tau2_sk_ST_U <- res_sk_ST_U[dtau2_sk_ST_U,]$median
LBtau2_sk_ST_U <- res_sk_ST_U[dtau2_sk_ST_U,]$q5
UBtau2_sk_ST_U <- res_sk_ST_U[dtau2_sk_ST_U,]$q95
Rhat_tau2sk_ST_U <- res_sk_ST_U[dtau2_sk_ST_U,]$rhat
prec_tau2_sk_ST_U <- UBtau2_sk_ST_U - LBtau2_sk_ST_U

dx_sk_ST_U <- grepl("xi",res_sk_ST_U$variable)

xi_sk_ST_U <- res_sk_ST_U[dx_sk_ST_U, ]$median
LBxi_sk_ST_U <- res_sk_ST_U[dx_sk_ST_U,]$q5
UBxi_sk_ST_U <- res_sk_ST_U[dx_sk_ST_U,]$q95
Rhat_sk_ST_U <-  res_sk_ST_U[dx_sk_ST_U,]$rhat

dxn_sk_ST_U <- grepl("pred",res_sk_ST_U$variable)
pred_sk_ST_U <- res_sk_ST_U[dxn_sk_ST_U,]$median
LBpred_sk_ST_U <- res_sk_ST_U[dxn_sk_ST_U,]$q5
UBpred_sk_ST_U <- res_sk_ST_U[dxn_sk_ST_U,]$q95
Rhat_predsk_ST_U <-  res_sk_ST_U[dxn_sk_ST_U,]$rhat
prec_pred_sk_ST_U <- UBpred_sk_ST_U - LBpred_sk_ST_U


skew_sk_ST_U <- c()
dskew_sk_ST_U <- grepl("skew",res_sk_ST_U$variable)

skew_sk_ST_U <- res_sk_ST_U[dskew_sk_ST_U,]$median
LBskew_sk_ST_U <- res_sk_ST_U[dskew_sk_ST_U,]$q5
UBskew_sk_ST_U <- res_sk_ST_U[dskew_sk_ST_U,]$q95
Rhat_skewsk_ST_U <- res_sk_ST_U[dskew_sk_ST_U,]$rhat

kurt_sk_ST_U <- c()
dkurt_sk_ST_U <- grepl("kurt",res_sk_ST_U$variable)

kurt_sk_ST_U <- res_sk_ST_U[dkurt_sk_ST_U,]$median
LBkurt_sk_ST_U <- res_sk_ST_U[dkurt_sk_ST_U,]$q5
UBkurt_sk_ST_U <- res_sk_ST_U[dkurt_sk_ST_U,]$q95
Rhat_kurtsk_ST_U <- res_sk_ST_U[dkurt_sk_ST_U,]$rhat


dnu_sk_ST_U <- grepl("nu",res_sk_ST_U$variable)

nu_sk_ST_U <- res_sk_ST_U[dnu_sk_ST_U,]$median
LBnu_sk_ST_U <- res_sk_ST_U[dnu_sk_ST_U,]$q5
UBnu_sk_ST_U <- res_sk_ST_U[dnu_sk_ST_U,]$q95
Rhat_nusk_ST_U <- res_sk_ST_U[dnu_sk_ST_U,]$rhat

res1_sk_ST_U<-data.frame(median=median_sk_ST_U, 
                         lowerCI=LBmu_sk_ST_U,
                         upperCI=UBmu_sk_ST_U,
                         prec_mu_sk_ST_U = prec_mu_sk_ST_U,
                         tau2=tau2_sk_ST_U,
                         l_tau2 = LBtau2_sk_ST_U,
                         u_tau2 = UBtau2_sk_ST_U,
                         prec_tau2_sk_ST_U = prec_tau2_sk_ST_U,
                         skew = skew_sk_ST_U,
                         l_skew = LBskew_sk_ST_U,
                         u_skew = UBskew_sk_ST_U,
                         kurt = kurt_sk_ST_U,
                         l_kurt = LBkurt_sk_ST_U,
                         u_kurt = UBkurt_sk_ST_U,
                         nu = nu_sk_ST_U,
                         l_nu = LBnu_sk_ST_U,
                         u_nu = UBnu_sk_ST_U,
                         Rhat_muST = Rhat_musk_ST_U,
                         Rhat_tau2ST = Rhat_tau2sk_ST_U,
                         Rhat_skewST = Rhat_skewsk_ST_U,
                         Rhat_kurtST = Rhat_kurtsk_ST_U,
                         Rhat_nuST = Rhat_nusk_ST_U,
                         pred = pred_sk_ST_U,
                         LBpred = LBpred_sk_ST_U,
                         UBpred = UBpred_sk_ST_U,
                         prec_pred_sk_ST_U = prec_pred_sk_ST_U,
                         Rhat_predST = Rhat_predsk_ST_U,
                         Pr_mu = Pr_mu,
                         Pr_pred = Pr_pred
)



# 
# 
# write.csv(res1_sk_ST_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\ST_U_pred_58\\sk_res_normal_sk_ST_U_SMDs.csv",row.names=FALSE )
# write.csv(deltask_ST_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\ST_U_pred_58\\sk_rel_eff_normal_sk_ST_U_SMDs.csv",row.names=FALSE )

