
## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

### set up for stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
## download the stan model
####  normal-t(U) model
url <- "https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/sk_mu_normal_t_U_SMDs.stan"
file_sk_t_distr_U  <- "sk_mu_normal_t_U_MDs.stan"

download.file(url,file_sk_t_distr_U  )

mod_sk_t_distr_U <- cmdstan_model(file_sk_t_distr_U )

fit_sk_t_distr_U <-  mod_sk_t_distr_U$sample(
  
  dati1  <-list(ns = nrow(pre_term_data_58),
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
  adapt_delta = 0.99,
  
)

sk_t_distr_U_res_SMDs <- as.data.frame(fit_sk_t_distr_U$summary())
View(sk_t_distr_U_res_SMDs)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fit_sk_t_distr_U$draws(c("mu", "pred")))

## traceplots ###
library(bayesplot)
library(posterior)

draws <- fit_sk_t_distr_U$draws()

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

library(bayesplot)
draws <- fit_sk_t_distr_U$draws()  # draws_array object

bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", "nu"))
delta_names <- sk_t_distr_U_res_SMDs$variable[grepl("delta", sk_t_distr_U_res_SMDs$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout

### Pr of mu <- 0
Pr_mu = mean(fit_sk_t_distr_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fit_sk_t_distr_U$draws("pred") < 0 )

dt_U <- grepl("delta",sk_t_distr_U_res_SMDs$variable)
deltatt_distr.res_SMDs <- sk_t_distr_U_res_SMDs[dt_U,]

plot(density(deltatt_distr.res_SMDs$median))

dmt_U <- grepl("mu",sk_t_distr_U_res_SMDs$variable)

median_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dmt_U,]$median
LBmu_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dmt_U,]$q5
UBmu_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dmt_U,]$q95
Rhat_mu_sk_t_distr_U <-  sk_t_distr_U_res_SMDs[dmt_U,]$rhat
prec_mu_sk_t_distr_U <- UBmu_sk_t_distr_U - LBmu_sk_t_distr_U


ddnt_U <- grepl("pred",sk_t_distr_U_res_SMDs$variable)

delta_new_sk_t_distr_U <- sk_t_distr_U_res_SMDs[ddnt_U,]$median
LBdelta_new_sk_t_distr_U <- sk_t_distr_U_res_SMDs[ddnt_U,]$q5
UBdelta_new_sk_t_distr_U <- sk_t_distr_U_res_SMDs[ddnt_U,]$q95
Rhat_delta_new_sk_t_distr_U <-  sk_t_distr_U_res_SMDs[ddnt_U,]$rhat
prec_delta_new_sk_t_distr_U <- UBdelta_new_sk_t_distr_U - LBdelta_new_sk_t_distr_U

dnut_U <- grepl("nu",sk_t_distr_U_res_SMDs$variable)
nu_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dnut_U,]$median
LBnu_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dnut_U,]$q5
UBnu_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dnut_U,]$q95
Rhat_nu_sk_t_distr_U <-  sk_t_distr_U_res_SMDs[dnut_U,]$rhat
prec_nu_sk_t_distr_U <- UBnu_sk_t_distr_U - LBnu_sk_t_distr_U

tau2t_U <- c()
dtau2t_U <- grepl("tau_sqr",sk_t_distr_U_res_SMDs$variable)

tau2_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dtau2t_U,]$median
LBtau2_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dtau2t_U,]$q5
UBtau2_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dtau2t_U,]$q95
Rhat_tau2_sk_t_distr_U <- sk_t_distr_U_res_SMDs[dtau2t_U,]$rhat
prec_tau2_sk_t_distr_U <- UBtau2_sk_t_distr_U - LBtau2_sk_t_distr_U

t_distr.res_SMDs1<-data.frame(median=median_sk_t_distr_U, 
                              lowerCI=LBmu_sk_t_distr_U,
                              upperCI=UBmu_sk_t_distr_U,
                              tau2=tau2_sk_t_distr_U,
                              l_tau2 = LBtau2_sk_t_distr_U,
                              u_tau2 = UBtau2_sk_t_distr_U,
                              nu_sk_t_distr_U = nu_sk_t_distr_U,
                              LBnu_sk_t_distr_U = LBnu_sk_t_distr_U,
                              UBnu_sk_t_distr_U = UBnu_sk_t_distr_U,
                              prec_mu_sk_t_distr_U = prec_mu_sk_t_distr_U,
                              prec_tau2_sk_t_distr_U = prec_tau2_sk_t_distr_U,
                              prec_delta_new_sk_t_distr_U = prec_delta_new_sk_t_distr_U,
                              Rhat_mut_distr = Rhat_mu_sk_t_distr_U,
                              Rhat_tau2t_distr = Rhat_tau2_sk_t_distr_U,
                              Rhat_tau2_sk_t_distr_U = Rhat_tau2_sk_t_distr_U,
                              delta_new = delta_new_sk_t_distr_U,
                              LBdelta_new = LBdelta_new_sk_t_distr_U,
                              UBdelta_new = UBdelta_new_sk_t_distr_U,
                              Rhat_delta_newt_distr = Rhat_delta_new_sk_t_distr_U,
                              Pr_mu = Pr_mu,
                              Pr_pred = Pr_pred
)

# 
# write.csv(t_distr.res_SMDs1 , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\t_U_pred_58\\sk_res_SMDs_normal_t_pred_58.csv",row.names=FALSE )  
# write.csv( deltatt_distr.res_SMDs , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\t_U_pred_58\\sk_rel_eff_normal_t_pred_58.csv",row.names=FALSE )  

