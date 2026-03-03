
## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

## analysis only with recent studies ###
# #keep onlt studies with publication year of >= 2016
# pre_term_data_58$numbers <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# pre_term_data_58$publication_year <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# 
# pre_term_data_58 <- pre_term_data_58[
#   pre_term_data_58$publication_year >= 2016,
# ]

### set up for stan model
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
check_cmdstan_toolchain()
## download the stan model
####  t-distribution model
url <- "https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/normal_t_U_SMDs.stan"
filet_distr_U <- "normal_t_U_SMDs.stan"

download.file(url, filet_distr_U)

modt_distr_U <- cmdstan_model(filet_distr_U)

fitt_distr_U <-  modt_distr_U$sample(
  
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

t_distr_U_res_SMDs <- as.data.frame(fitt_distr_U$summary())
View(t_distr_U_res_SMDs)
#### To visualize the predictive distribution
bayesplot::color_scheme_set("brightblue")
bayesplot::mcmc_dens(fitt_distr_U$draws(c("mu", "pred")))

## traceplots ###
library(bayesplot)
library(posterior)

draws <- fitt_distr_U$draws()

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

library(bayesplot)
draws <- fitt_distr_U$draws()  

bayesplot::mcmc_trace(draws, pars = c("mu","tau_sqr", "pred", "nu"))
delta_names <- t_distr_U_res_SMDs$variable[grepl("delta", t_distr_U_res_SMDs$variable)]

bayesplot::mcmc_trace(draws, pars = delta_names)

par(mfrow = c(1, 1))  # Reset layout

### Pr of mu <- 0
Pr_mu = mean(fitt_distr_U$draws("mu") < 0 )

### Pr of pred1 <- 0
Pr_pred = mean(fitt_distr_U$draws("pred") < 0 )

dt_U <- grepl("delta",t_distr_U_res_SMDs$variable)
deltatt_distr.res_SMDs <- t_distr_U_res_SMDs[dt_U,]

plot(density(deltatt_distr.res_SMDs$median))

dmt_U <- grepl("mu",t_distr_U_res_SMDs$variable)

median_t_distr_U <- t_distr_U_res_SMDs[dmt_U,]$median
LBmu_t_distr_U <- t_distr_U_res_SMDs[dmt_U,]$q5
UBmu_t_distr_U <- t_distr_U_res_SMDs[dmt_U,]$q95
Rhat_mut_distr_U <-  t_distr_U_res_SMDs[dmt_U,]$rhat
prec_mu_t_distr_U <- UBmu_t_distr_U - LBmu_t_distr_U


library(metafor)
t_distr_U_SMDs <- rma(yi = deltatt_distr.res_SMDs$median , vi = (deltatt_distr.res_SMDs$sd)^2 , 
                      method = "REML", slab = pre_term_data_58$Study)
forest(t_distr_U_SMDs, cex = 0.5)


ddnt_U <- grepl("pred",t_distr_U_res_SMDs$variable)

delta_new_t_distr_U <- t_distr_U_res_SMDs[ddnt_U,]$median
LBdelta_new_t_distr_U <- t_distr_U_res_SMDs[ddnt_U,]$q5
UBdelta_new_t_distr_U <- t_distr_U_res_SMDs[ddnt_U,]$q95
Rhat_delta_newt_distr_U <-  t_distr_U_res_SMDs[ddnt_U,]$rhat
prec_delta_new_t_distr_U <- UBdelta_new_t_distr_U - LBdelta_new_t_distr_U

dnut_U <- grepl("nu",t_distr_U_res_SMDs$variable)
nu_t_distr_U <- t_distr_U_res_SMDs[dnut_U,]$median
LBnu_t_distr_U <- t_distr_U_res_SMDs[dnut_U,]$q5
UBnu_t_distr_U <- t_distr_U_res_SMDs[dnut_U,]$q95
Rhat_nut_distr_U <-  t_distr_U_res_SMDs[dnut_U,]$rhat
prec_nu_t_distr_U <- UBnu_t_distr_U - LBnu_t_distr_U

tau2t_U <- c()
dtau2t_U <- grepl("tau_sqr",t_distr_U_res_SMDs$variable)

tau2_t_distr_U <- t_distr_U_res_SMDs[dtau2t_U,]$median
LBtau2_t_distr_U <- t_distr_U_res_SMDs[dtau2t_U,]$q5
UBtau2_t_distr_U <- t_distr_U_res_SMDs[dtau2t_U,]$q95
Rhat_tau2t_distr_U <- t_distr_U_res_SMDs[dtau2t_U,]$rhat
prec_tau2_t_distr_U <- UBtau2_t_distr_U - LBtau2_t_distr_U

t_distr.res_SMDs1<-data.frame(median=median_t_distr_U, 
                              lowerCI=LBmu_t_distr_U,
                              upperCI=UBmu_t_distr_U,
                              tau2=tau2_t_distr_U,
                              l_tau2 = LBtau2_t_distr_U,
                              u_tau2 = UBtau2_t_distr_U,
                              nu_t_distr_U = nu_t_distr_U,
                              LBnu_t_distr_U = LBnu_t_distr_U,
                              UBnu_t_distr_U = UBnu_t_distr_U,
                              prec_mu_t_distr_U = prec_mu_t_distr_U,
                              prec_tau2_t_distr_U = prec_tau2_t_distr_U,
                              prec_delta_new_t_distr_U = prec_delta_new_t_distr_U,
                              Rhat_mut_distr = Rhat_mut_distr_U,
                              Rhat_tau2t_distr = Rhat_tau2t_distr_U,
                              Rhat_tau2t_distr_U = Rhat_tau2t_distr_U,
                              delta_new = delta_new_t_distr_U,
                              LBdelta_new = LBdelta_new_t_distr_U,
                              UBdelta_new = UBdelta_new_t_distr_U,
                              Rhat_delta_newt_distr = Rhat_delta_newt_distr_U,
                              Pr_mu = Pr_mu,
                              Pr_pred = Pr_pred
)


# write.csv(t_distr.res_SMDs1 , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\t_U_pred_58\\res_SMDs_normal_t_pred_58.csv",row.names=FALSE )  
# write.csv( deltatt_distr.res_SMDs , "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\STAN_models_58\\t_U_pred_58\\rel_eff_normal_t_pred_58.csv",row.names=FALSE )  

