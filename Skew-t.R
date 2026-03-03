## download the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

## analysis only with recent studies ###
# #keep only studies with publication year of >= 2016
# pre_term_data_58$numbers <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# pre_term_data_58$publication_year <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# 
# pre_term_data_58 <- pre_term_data_58[
#   pre_term_data_58$publication_year >= 2016, 
# ]

library(metafor)
# ################################## Skew-t model #############################
## prepare for the stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
### donload the stan model
url <-"https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/normal_ST_U_SMDs.stan"

fileST_U <- "normal_ST_U_SMDs.stan"

download.file(url, fileST_U)

modST_U <- cmdstan_model(fileST_U)

fitST_U <-  modST_U$sample(
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

res_ST_U <- as.data.frame(fitST_U$summary())
View(res_ST_U)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitST_U$draws(c("mu", "pred")))


library(bayesplot)
draws <- fitST_U$draws()  

### traceplots 
bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", "nu",
                                      "xi", "skew","kurt"))
delta_names <- fitST_U$variable[grepl("delta",fitST_U$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout

### Pr of mu <- 0
Pr_mu = mean(fitST_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fitST_U$draws("pred") < 0 )


d_ST_U <- grepl("delta",res_ST_U$variable)

deltaST_U <- res_ST_U[d_ST_U,]

plot(density(deltaST_U$median))

dm_ST_U <- grepl("mu",res_ST_U$variable)

median_ST_U <- res_ST_U[dm_ST_U,]$median
LBmu_ST_U <- res_ST_U[dm_ST_U,]$q5
UBmu_ST_U <- res_ST_U[dm_ST_U,]$q95
Rhat_muST_U <-  res_ST_U[dm_ST_U,]$rhat
prec_mu_ST_U <- UBmu_ST_U - LBmu_ST_U

dtau2_ST_U <- grepl("tau_sqr",res_ST_U$variable)

tau2_ST_U <- res_ST_U[dtau2_ST_U,]$median
LBtau2_ST_U <- res_ST_U[dtau2_ST_U,]$q5
UBtau2_ST_U <- res_ST_U[dtau2_ST_U,]$q95
Rhat_tau2ST_U <- res_ST_U[dtau2_ST_U,]$rhat
prec_tau2_ST_U <- UBtau2_ST_U - LBtau2_ST_U

dx_ST_U <- grepl("xi",res_ST_U$variable)

xi_ST_U <- res_ST_U[dx_ST_U, ]$median
LBxi_ST_U <- res_ST_U[dx_ST_U,]$q5
UBxi_ST_U <- res_ST_U[dx_ST_U,]$q95
Rhat_ST_U <-  res_ST_U[dx_ST_U,]$rhat

dxn_ST_U <- grepl("pred",res_ST_U$variable)
pred_ST_U <- res_ST_U[dxn_ST_U,]$median
LBpred_ST_U <- res_ST_U[dxn_ST_U,]$q5
UBpred_ST_U <- res_ST_U[dxn_ST_U,]$q95
Rhat_predST_U <-  res_ST_U[dxn_ST_U,]$rhat
prec_pred_ST_U <- UBpred_ST_U - LBpred_ST_U


skew_ST_U <- c()
dskew_ST_U <- grepl("skew",res_ST_U$variable)

skew_ST_U <- res_ST_U[dskew_ST_U,]$median
LBskew_ST_U <- res_ST_U[dskew_ST_U,]$q5
UBskew_ST_U <- res_ST_U[dskew_ST_U,]$q95
Rhat_skewST_U <- res_ST_U[dskew_ST_U,]$rhat

kurt_ST_U <- c()
dkurt_ST_U <- grepl("kurt",res_ST_U$variable)

kurt_ST_U <- res_ST_U[dkurt_ST_U,]$median
LBkurt_ST_U <- res_ST_U[dkurt_ST_U,]$q5
UBkurt_ST_U <- res_ST_U[dkurt_ST_U,]$q95
Rhat_kurtST_U <- res_ST_U[dkurt_ST_U,]$rhat


dnu_ST_U <- grepl("nu",res_ST_U$variable)

nu_ST_U <- res_ST_U[dnu_ST_U,]$median
LBnu_ST_U <- res_ST_U[dnu_ST_U,]$q5
UBnu_ST_U <- res_ST_U[dnu_ST_U,]$q95
Rhat_nuST_U <- res_ST_U[dnu_ST_U,]$rhat

res1_ST_U<-data.frame(median=median_ST_U, 
                      lowerCI=LBmu_ST_U,
                      upperCI=UBmu_ST_U,
                      prec_mu_ST_U = prec_mu_ST_U,
                      tau2=tau2_ST_U,
                      l_tau2 = LBtau2_ST_U,
                      u_tau2 = UBtau2_ST_U,
                      prec_tau2_ST_U = prec_tau2_ST_U,
                      skew = skew_ST_U,
                      l_skew = LBskew_ST_U,
                      u_skew = UBskew_ST_U,
                      kurt = kurt_ST_U,
                      l_kurt = LBkurt_ST_U,
                      u_kurt = UBkurt_ST_U,
                      nu = nu_ST_U,
                      l_nu = LBnu_ST_U,
                      u_nu = UBnu_ST_U,
                      Rhat_muST = Rhat_muST_U,
                      Rhat_tau2ST = Rhat_tau2ST_U,
                      Rhat_skewST = Rhat_skewST_U,
                      Rhat_kurtST = Rhat_kurtST_U,
                      Rhat_nuST = Rhat_nuST_U,
                      pred = pred_ST_U,
                      LBpred = LBpred_ST_U,
                      UBpred = UBpred_ST_U,
                      prec_pred_ST_U = prec_pred_ST_U,
                      Rhat_predST = Rhat_predST_U,
                      Pr_mu = Pr_mu,
                      Pr_pred = Pr_pred
)





# write.csv(res1_ST_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\ST_U_pred_58\\res_normal_ST_U_SMDs.csv",row.names=FALSE )
#write.csv(deltaST_U , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\ST_U_pred_58\\rel_eff_normal_ST_U_SMDs.csv",row.names=FALSE )  

