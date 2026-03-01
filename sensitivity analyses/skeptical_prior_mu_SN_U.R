## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

library(dplyr)
library(metafor)
# ################################## normal-SN(U) model #############################
## prepare for the stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
#### download the stan model
url <-"https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/sk_mu_normal_SN_U_SMDs.stan"

filesk_SN_U <- "sk_mu_normal_SN_U_SMDs.stan"

download.file(url, filesk_SN_U)

modsk_SN_U <- cmdstan_model(filesk_SN_U)


fitsk_SN_U <-  modsk_SN_U$sample(
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

sk_sn.res_U <- as.data.frame(fitsk_SN_U$summary())
View(sk_sn.res_U)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitsk_SN_U$draws(c("mu", "pred")))

### traceplots
draws <- fitsk_SN_U$draws()

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

library(bayesplot)
draws <- fitsk_SN_U$draws() 

bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", 
                                      "xi", "skew","kurt"))
delta_names <- sk_sn.res_U$variable[grepl("delta",sk_sn.res_U$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout


### Pr of mu <- 0
Pr_mu = mean(fitsk_SN_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fitsk_SN_U$draws("pred") < 0 )


d_U <- grepl("delta",sk_sn.res_U$variable)

deltask_SN_U <- sk_sn.res_U[d_U,]

plot(density(deltask_SN_U$median))

dm_U <- grepl("mu",sk_sn.res_U$variable)

median_sk_SN_U <- sk_sn.res_U[dm_U,]$median
LBmu_sk_SN_U <- sk_sn.res_U[dm_U,]$q5
UBmu_sk_SN_U <- sk_sn.res_U[dm_U,]$q95
Rhat_musk_SN_U <-  sk_sn.res_U[dm_U,]$rhat
prec_mu_U <- UBmu_sk_SN_U - LBmu_sk_SN_U

dtau2_U <- grepl("tau_sqr",sk_sn.res_U$variable)

tau2_sk_SN_U <- sk_sn.res_U[dtau2_U,]$median
LBtau2_sk_SN_U <- sk_sn.res_U[dtau2_U,]$q5
UBtau2_sk_SN_U <- sk_sn.res_U[dtau2_U,]$q95
Rhat_tau2sk_SN_U <- sk_sn.res_U[dtau2_U,]$rhat
prec_tau2_U <- UBtau2_sk_SN_U - LBtau2_sk_SN_U

dx_U <- grepl("xi",sk_sn.res_U$variable)

xi_sk_SN_U <- sk_sn.res_U[dx_U, ]$median
LBxi_sk_SN_U <- sk_sn.res_U[dx_U,]$q5
UBxi_sk_SN_U <- sk_sn.res_U[dx_U,]$q95
Rhat_sk_SN_U <-  sk_sn.res_U[dx_U,]$rhat

dxn_sk_SN_U <- grepl("pred",sk_sn.res_U$variable)
pred_sk_SN_U <- sk_sn.res_U[dxn_sk_SN_U,]$median
LBpred_sk_SN_U <- sk_sn.res_U[dxn_sk_SN_U,]$q5
UBpred_sk_SN_U <- sk_sn.res_U[dxn_sk_SN_U,]$q95
Rhat_predsk_SN_U <-  sk_sn.res_U[dxn_sk_SN_U,]$rhat
prec_pred_U <- UBpred_sk_SN_U - LBpred_sk_SN_U


skew_sk_SN_U <- c()
dskew_sk_SN_U <- grepl("skew",sk_sn.res_U$variable)

skew_sk_SN_U <- sk_sn.res_U[dskew_sk_SN_U,]$median
LBskew_sk_SN_U <- sk_sn.res_U[dskew_sk_SN_U,]$q5
UBskew_sk_SN_U <- sk_sn.res_U[dskew_sk_SN_U,]$q95
Rhat_skewsk_SN_U <- sk_sn.res_U[dskew_sk_SN_U,]$rhat

kurt_sk_SN_U <- c()
dkurt_sk_SN_U <- grepl("kurt",sk_sn.res_U$variable)

kurt_sk_SN_U <- sk_sn.res_U[dkurt_sk_SN_U,]$median
LBkurt_sk_SN_U <- sk_sn.res_U[dkurt_sk_SN_U,]$q5
UBkurt_sk_SN_U <- sk_sn.res_U[dkurt_sk_SN_U,]$q95
Rhat_kurtsk_SN_U <- sk_sn.res_U[dkurt_sk_SN_U,]$rhat

res1_sk_SN_U <-data.frame(median=median_sk_SN_U, 
                          lowerCI=LBmu_sk_SN_U,
                          upperCI=UBmu_sk_SN_U,
                          prec_mu_U = prec_mu_U,
                          tau2=tau2_sk_SN_U,
                          l_tau2 = LBtau2_sk_SN_U,
                          u_tau2 = UBtau2_sk_SN_U,
                          prec_tau2_U,
                          skew = skew_sk_SN_U,
                          l_skew = LBskew_sk_SN_U,
                          u_skew = UBskew_sk_SN_U,
                          kurt = kurt_sk_SN_U,
                          l_kurt = LBkurt_sk_SN_U,
                          u_kurt = UBkurt_sk_SN_U,
                          Rhat_muSN = Rhat_musk_SN_U,
                          Rhat_tau2SN = Rhat_tau2sk_SN_U,
                          Rhat_skewSN = Rhat_skewsk_SN_U,
                          Rhat_kurtSN = Rhat_kurtsk_SN_U,
                          pred = pred_sk_SN_U,
                          LBpred = LBpred_sk_SN_U,
                          UBpred = UBpred_sk_SN_U,
                          prec_pred_U = prec_pred_U,
                          Rhat_predSN = Rhat_predsk_SN_U,
                          Pr_mu = Pr_mu,
                          Pr_pred = Pr_pred
)


# 
# write.csv(res1_sk_SN_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\SN_U_pred_58\\sk_res_normal_U_SMDs.csv",row.names=FALSE )  
# write.csv( deltask_SN_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\SN_U_pred_58\\sk_rel_eff_normal_U_SMDs.csv",row.names=FALSE )  

