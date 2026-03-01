## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

## analysis only with recent studies ###
# #keep only studies with publication year of >= 2016
# pre_term_data_58$numbers <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# pre_term_data_58$publication_year <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# 
# pre_term_data_58 <- pre_term_data_58[
#   pre_term_data_58$publication_year >= 2016,
# ]

library(dplyr)
library(metafor)
# ################################## normal-SN(U) model #############################
## prepare for the stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
#### download the stan model
url <-"https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/normal_SN_U_SMDs.stan"

fileSN_U <- "normal_t_U_MDs.stan"

download.file(url, fileSN_U)

modSN_U <- cmdstan_model(fileSN_U)

fitSN_U <-  modSN_U$sample(
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

sn.res_U <- as.data.frame(fitSN_U$summary())
View(sn.res_U)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitSN_U$draws(c("mu", "pred")))

### traceplots
draws <- fitSN_U$draws()

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

library(bayesplot)
draws <- fitSN_U$draws() 

bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", 
                                      "xi", "skew","kurt"))
delta_names <- sn.res_U$variable[grepl("delta",sn.res_U$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout


### Pr of mu <- 0
Pr_mu = mean(fitSN_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fitSN_U$draws("pred") < 0 )


d_U <- grepl("delta",sn.res_U$variable)

deltaSN_U <- sn.res_U[d_U,]

plot(density(deltaSN_U$median))

dm_U <- grepl("mu",sn.res_U$variable)

median_SN_U <- sn.res_U[dm_U,]$median
LBmu_SN_U <- sn.res_U[dm_U,]$q5
UBmu_SN_U <- sn.res_U[dm_U,]$q95
Rhat_muSN_U <-  sn.res_U[dm_U,]$rhat
prec_mu_U <- UBmu_SN_U - LBmu_SN_U

dtau2_U <- grepl("tau_sqr",sn.res_U$variable)

tau2_SN_U <- sn.res_U[dtau2_U,]$median
LBtau2_SN_U <- sn.res_U[dtau2_U,]$q5
UBtau2_SN_U <- sn.res_U[dtau2_U,]$q95
Rhat_tau2SN_U <- sn.res_U[dtau2_U,]$rhat
prec_tau2_U <- UBtau2_SN_U - LBtau2_SN_U

dx_U <- grepl("xi",sn.res_U$variable)

xi_SN_U <- sn.res_U[dx_U, ]$median
LBxi_SN_U <- sn.res_U[dx_U,]$q5
UBxi_SN_U <- sn.res_U[dx_U,]$q95
Rhat_SN_U <-  sn.res_U[dx_U,]$rhat

dxn_SN_U <- grepl("pred",sn.res_U$variable)
pred_SN_U <- sn.res_U[dxn_SN_U,]$median
LBpred_SN_U <- sn.res_U[dxn_SN_U,]$q5
UBpred_SN_U <- sn.res_U[dxn_SN_U,]$q95
Rhat_predSN_U <-  sn.res_U[dxn_SN_U,]$rhat
prec_pred_U <- UBpred_SN_U - LBpred_SN_U


skew_SN_U <- c()
dskew_SN_U <- grepl("skew",sn.res_U$variable)

skew_SN_U <- sn.res_U[dskew_SN_U,]$median
LBskew_SN_U <- sn.res_U[dskew_SN_U,]$q5
UBskew_SN_U <- sn.res_U[dskew_SN_U,]$q95
Rhat_skewSN_U <- sn.res_U[dskew_SN_U,]$rhat

kurt_SN_U <- c()
dkurt_SN_U <- grepl("kurt",sn.res_U$variable)

kurt_SN_U <- sn.res_U[dkurt_SN_U,]$median
LBkurt_SN_U <- sn.res_U[dkurt_SN_U,]$q5
UBkurt_SN_U <- sn.res_U[dkurt_SN_U,]$q95
Rhat_kurtSN_U <- sn.res_U[dkurt_SN_U,]$rhat

res1_SN_U <-data.frame(median=median_SN_U, 
                       lowerCI=LBmu_SN_U,
                       upperCI=UBmu_SN_U,
                       prec_mu_U = prec_mu_U,
                       tau2=tau2_SN_U,
                       l_tau2 = LBtau2_SN_U,
                       u_tau2 = UBtau2_SN_U,
                       prec_tau2_U,
                       skew = skew_SN_U,
                       l_skew = LBskew_SN_U,
                       u_skew = UBskew_SN_U,
                       kurt = kurt_SN_U,
                       l_kurt = LBkurt_SN_U,
                       u_kurt = UBkurt_SN_U,
                       Rhat_muSN = Rhat_muSN_U,
                       Rhat_tau2SN = Rhat_tau2SN_U,
                       Rhat_skewSN = Rhat_skewSN_U,
                       Rhat_kurtSN = Rhat_kurtSN_U,
                       pred = pred_SN_U,
                       LBpred = LBpred_SN_U,
                       UBpred = UBpred_SN_U,
                       prec_pred_U = prec_pred_U,
                       Rhat_predSN = Rhat_predSN_U,
                       Pr_mu = Pr_mu,
                       Pr_pred = Pr_pred
)



# write.csv(res1_SN_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\SN_U_pred_58\\res_normal_U_SMDs.csv",row.names=FALSE )  
# write.csv( deltaSN_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\SN_U_pred_58\\rel_eff_normal_U_SMDs.csv",row.names=FALSE )  

