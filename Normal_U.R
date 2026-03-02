##load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

## analysis only with recent studies ###
# # #keep only studies with publication year of >= 2016
# pre_term_data_58$numbers <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# pre_term_data_58$publication_year <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# 
# pre_term_data_58 <- pre_term_data_58[
#   pre_term_data_58$publication_year >= 2016,
# ]

######### NORMAL MODEL FOR CONTINUOUS OUTCOMES ##############
#### Normal model ###
library(R2jags)
set.seed(1508)
writeLines("
  model{
    for(i in 1:ns){ 
    
      m[i] ~ dnorm(0,.0001)
      
      # Mean outcome for group 1
      y1[i] ~ dnorm(mu1[i], prec1[i])
      mu1[i] <- m[i]
      
      # Mean outcome for group 2
      y2[i] ~ dnorm(mu2[i], prec2[i])
      mu2[i] <- m[i] + delta12[i]*pooled_sd[i] ## delta12 = (m2-m1)/pooled_sd
      
      # Between-study effect
      delta12[i] ~ dnorm(mu, tau1)
      
      prec1[i] <- n1[i]/ pow(sd1[i], 2)
      prec2[i] <- n2[i]/ pow(sd2[i], 2)

      # pooled standard deviation
      pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))
    }
    
    # Priors
    mu ~ dnorm(0, .0001) 
    
    ## or skeptical prior on mu
    ### normal(b) model ###
    # mu~ dnorm(0, prec_mu)
    # prec_mu = 1 / (0.255)^2
    
    tau1 <- 1 / tau_sqr
    tau_sqr <- tau * tau
    tau ~ dunif(0,10)

    # Prediction for a new effect size
    delta_new ~ dnorm(mu, tau1)
  } 
", con = "Normal_U_SMDs.txt")

modfile = 'Normal_U_SMDs.txt'

dati1  <-list(ns = nrow(pre_term_data_58),
              y1 = pre_term_data_58$mean_FT,
              y2 = pre_term_data_58$mean_EPT.VPT,
              sd1 = pre_term_data_58$sd_FT,
              sd2 = pre_term_data_58$sd_EPT.VPT,
              n1 = pre_term_data_58$n_FT,
              n2 = pre_term_data_58$n_EPT.VPT)


run.modelNU_SMDs = jags(
  data = dati1,
  inits = NULL,
  parameters.to.save = c(
    "mu",
    "tau",
    "delta12",
    "delta_new"
  ),
  n.chains = 2,
  n.iter = 50000,
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
)

resultsNU1_SMDs = as.data.frame(run.modelNU_SMDs$BUGSoutput$summary) 
resultsNU_SMDs = round(resultsNU1_SMDs , 2)
View(resultsNU_SMDs)

##### to check the convergence ###########
library(coda)
mcmc_obj <- as.mcmc(run.modelNU_SMDs)

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^delta12\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters


par(mfrow = c(4, 3))  # Adjust margins if needed
traceplot(mcmc_obj[, c("mu", "tau","delta_new"
)])

par(mfrow = c(1, 1))  # Reset layout


#####STUDY-SPECIFIC EFFECTS
fd <- grepl("delta12", row.names(resultsNU1_SMDs))
md <- resultsNU1_SMDs[fd,]
rel_eff_NU_SMDs <- md$`50%`
LB_rel_eff_NU_SMDs <- md$`2.5%`
UB_rel_eff_NU_SMDs <- md$`97.5%`
sd_rel_eff_NU_SMDs <- md$sd
Rhat_deltaN_NU_SMDs <- md$Rhat

######## TO VISUALIZE THE DENSITY PLOT OF STUDY-SPECIFIC EFFECTS ##########
plot(density(rel_eff_NU_SMDs))

#### The density plot of pred ####
plot(density(run.modelNU_SMDs$BUGSoutput$sims.matrix[  ,"delta_new"]))

### Pr of mu <0 and Pr of pred <- 0 #########
Pr_mu_NU_SMDs = mean(run.modelNU_SMDs$BUGSoutput$sims.matrix[  ,"mu"] < 0 )
Pr_pred_NU_SMDs = mean(run.modelNU_SMDs$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )


f1NU <- grepl("mu", row.names(resultsNU1_SMDs))
m1NU <- resultsNU1_SMDs[f1NU,]
mu_NU_SMDs <- m1NU$`50%` 
LB_mu_NU_SMDs <- m1NU$`2.5%`
UB_mu_NU_SMDs <- m1NU$`97.5%`
Rhat_mu_NU_SMDs <- m1NU$Rhat
prec_mu_NU_SMDs <-  UB_mu_NU_SMDs - LB_mu_NU_SMDs


f2NU <- grepl("tau", row.names(resultsNU1_SMDs))
m2NU <- resultsNU1_SMDs[f2NU,]
tau_NU_SMDs <- m2NU$`50%`
tau2_NU_SMDs <- (m2NU$`50%`)^2
LB_tau2_NU_SMDs <- (m2NU$`2.5%`)^2
UB_tau2_NU_SMDs <- (m2NU$`97.5%`)^2
Rhat_tau_NU_SMDs <- m2NU$Rhat
tau2_NU_SMDs = round(tau2_NU_SMDs,2)
prec_tau2_NU_SMDs <-  UB_tau2_NU_SMDs - LB_tau2_NU_SMDs


fdNU <- grepl("delta12", row.names(resultsNU1_SMDs))
mdNU <- resultsNU1_SMDs[fdNU,]
rel_eff_NU_SMDs <- mdNU$`50%`
LB_rel_eff_NU_SMDs <- mdNU$`2.5%`
UB_rel_eff_NU_SMDs <- mdNU$`97.5%`
sd_rel_eff_NU_SMDs <- mdNU$sd
Rhat_delta_NU_SMDs <- mdNU$Rhat

fdnNU <- grepl("delta_new", row.names(resultsNU1_SMDs))
mdnNU <- resultsNU1_SMDs[fdnNU,]
pred_NU_SMDs <- mdnNU$`50%` 
LB_pred_NU_SMDs <- mdnNU$`2.5%`
UB_pred_NU_SMDs <- mdnNU$`97.5%`
Rhat_pred_NU_SMDs <- mdnNU$Rhat
prec_pred_NU_SMDs<-  UB_pred_NU_SMDs - LB_pred_NU_SMDs


lista_NU_SMDs <- cbind.data.frame(rel_eff_NU_SMDs,sd_rel_eff_NU_SMDs,LB_rel_eff_NU_SMDs, 
                                  UB_rel_eff_NU_SMDs, Rhat_delta_NU_SMDs)


res_NU_SMDs <- cbind.data.frame(mu_NU_SMDs, LB_mu_NU_SMDs, UB_mu_NU_SMDs, tau2_NU_SMDs,
                                LB_tau2_NU_SMDs, UB_tau2_NU_SMDs, pred_NU_SMDs, LB_pred_NU_SMDs,
                                UB_pred_NU_SMDs, Rhat_mu_NU_SMDs, Rhat_tau_NU_SMDs, Rhat_pred_NU_SMDs,
                                prec_mu_NU_SMDs,prec_tau2_NU_SMDs , prec_pred_NU_SMDs, 
                                Pr_mu_NU_SMDs, Pr_pred_NU_SMDs)

# write.csv(res_NU_SMDs, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\Nomal_U_58\\2N_U_res_58SMDs.csv",row.names=FALSE )  
# write.csv(lista_NU_SMDs, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\Nomal_U_58\\2rel_eff_N_U_58SMDs.csv",row.names=FALSE )  
# 


