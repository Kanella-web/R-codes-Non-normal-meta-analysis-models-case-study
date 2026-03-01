## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

### Binomial-DPp-n-(Unif/Gamma) ###
library(R2jags)
set.seed(1508)
cat("
model {
  for (i in 1:ns) { 
      m[i] ~ dnorm(0,.0001)
    
    # Mean outcome for group 1
      y1[i] ~ dnorm(mu1[i], prec1[i])
      mu1[i] <- m[i]
      
      # Mean outcome for group 2
      y2[i] ~ dnorm(mu2[i], prec2[i])
      mu2[i] <- m[i] + delta12[i]*pooled_sd[i]
      
    
      delta12[i] <- theta[Z[i]]
    
      Z[i] ~ dcat(p[]) #Z is an integer variable
      
      # Precision parameters for group 1 and 2 based on observed standard deviations
      prec1[i] <- n1[i]/ pow(sd1[i], 2)
      prec2[i] <- n2[i]/ pow(sd2[i], 2)


      # Calculate pooled standard deviation
      pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))
  }

  # Constructive DP
  # stick-breaking prior
  p[1] <- r[1]
  for (j in 2:(N-1)) {
    p[j] <- r[j] * (1 - r[j-1]) * p[j-1] / r[j-1]
  }
  for (k in 1:(N-1)) {
    r[k] ~ dbeta(1, alpha) T(0, 0.99)
  }
  # assumption to ensure sum p[] is 1 Ishwaran truncation
  ps <- sum(p[1:(N-1)])
  for (k in N:N) {
    p[k] <- 1 - ps
  }

  # Base distribution with normal prior
  for (k in 1:N) {
    
    theta[k] ~ dnorm(basemu, basetau1)  # normal base distribution
  }

   # Priors
   
   basetau1 <- 1 / basetau_sqr  
   basetau_sqr <- basetau * basetau
   basetau ~ dunif(0,10)

   basemu ~ dnorm(0, 0.0001)  
   
  # DPP parameter prior
  alpha ~ dgamma(1,1)
  
  # Random effects distribution mean
  for (i in 1:N) {
    meancl[i] <- p[i] * theta[i]
  }
  poptrue <- sum(meancl[])  # E[X]

  # Random effects distribution variance
  for (i in 1:N) {
    mom2[i] <- p[i] * theta[i] * theta[i]  # E[X2]
  }
  mom2.true <- sum(mom2[])
  var.true <- mom2.true - (poptrue * poptrue)  # E[X2] - E[X]^2

  # Programming for calculating summary statistics
  for (i in 1:ns) {
    for (j in 1:N) {
      SC[i, j] <- equals(j, Z[i])
    }
  }
  
   ## co-clustering probabilities
  for (i in 1:ns) {
    for (j in 1:ns) {
     equalsmatrix[i, j] <- equals(Z[i], Z[j])
     
    }
  equalsres[i] <- sum(equalsmatrix[i, ])
  
}


  # Total clusters K
  for (j in 1:N) {
    cl[j] <- step(sum(SC[, j]) - 1)
  }
  K <- sum(cl[])
  
  pred ~ dcat(p[1:N]) ##First randomly assigning the new study to one of the mixture components, according to the estimated weights p.
  delta_new <- theta[pred] ##Then drawing the new study's effect size from the random effects distribution.
  
}", file = "DPp_n_U_G_normal_base.txt")
modfile = 'DPp_n_U_G_normal_base.txt'

run.modelDPpn_U_G_normal_base = jags(
  dati1  <-list(ns = nrow(pre_term_data_58),
                y1 = pre_term_data_58$mean_FT,
                y2 = pre_term_data_58$mean_EPT.VPT,
                sd1 = pre_term_data_58$sd_FT,
                sd2 = pre_term_data_58$sd_EPT.VPT,
                n1 = pre_term_data_58$n_FT,
                n2 = pre_term_data_58$n_EPT.VPT,
                N = nrow(pre_term_data_58)
  )  ,  
  inits = NULL,
  parameters.to.save = c(
    "basemu", ## mu of the Normal base  distribution
    "basetau", ## tau of the Normal base distribution
    "poptrue",   ## overall mu
    "var.true", ## between-study heterogeneity
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "theta",
    "alpha",  ## concentration parameter
    "SC" ,    ## probability of each cluster assignment
    "pred",
    "delta_new",
    "equalsres",
    "equalsmatrix"
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

DPpresults_U_G_n_normal_base <- as.data.frame(run.modelDPpn_U_G_normal_base$BUGSoutput$summary) 

DPp_results_U_G_n_normal_base <- as.data.frame(run.modelDPpn_U_G_normal_base$BUGSoutput$summary) 

DPp_results_U_G_n_normal_base <-round(DPpresults_U_G_n_normal_base, digits=2) 

View(DPpresults_U_G_n_normal_base)

##### to check the convergence ###########
library(coda)
mcmc_obj <- as.mcmc(run.modelDPpn_U_G_normal_base)

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^theta\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^delta12\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters


par(mfrow = c(4, 3))  # Adjust margins if needed
traceplot(mcmc_obj[, c("basemu", "basetau", "alpha", "delta_new", 
                       "poptrue", 
                       "var.true")])

par(mfrow = c(1, 1))  # Reset layout

######## Pr(mu<0) #######
Pr_mu_DPp_UGn_normal_base = mean(run.modelDPpn_U_G_normal_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )
######## Pr(pred<0) #######
Pr_delta_new_DPp_UGn_normal_base = mean(run.modelDPpn_U_G_normal_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

### extract co-clustering probabilities
sims_equals <- run.modelDPpn_U_G_normal_base$BUGSoutput$sims.matrix
equals <- sims_equals[, grep("^equalsmatrix\\[", colnames(sims_equals))]

## create the 3dimentional array
ns <- nrow(pre_term_data_58)
niter <- nrow(equals)
equals_array <- array(
  equals,
  dim = c(niter, ns, ns)
)
dim(equals_array)

## create the 3dimentional matrix
distance.studies <- apply(equals_array, c(2, 3), mean)

##define the columns and the rows accordingly
pre_term_data_58$Study
study_names <- pre_term_data_58$Study
rownames(distance.studies) <- study_names
colnames(distance.studies) <- study_names

## heatmal
library(pheatmap)
library(RColorBrewer)
pheatmap(
  distance.studies,
  cellwidth = 15,
  color = brewer.pal(9, "Blues"),
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  main = "Co-clustering probabilities",
  treeheight_row = 0,
  treeheight_col = 0,
  clustering_distance_cols = as.dist(1 - distance.studies),
  clustering_distance_rows = as.dist(1 - distance.studies)
)



fdd <- grepl("delta12", row.names(DPpresults_U_G_n_normal_base))
mdd <- DPpresults_U_G_n_normal_base[fdd,]
rel_effDPp_UGn_normal_base <- mdd$`50%`
LB_rel_effDPp_UGn_normal_base <- mdd$`2.5%`
UB_rel_effDPp_UGn_normal_base <- mdd$`97.5%`
sd_rel_effDPp_UGn_normal_base <- mdd$sd
Rhat_deltaDPp_UGn_normal_base <- mdd$Rhat


#### PLOT THE DENSITY OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPp_UGn_normal_base))


########### SAVE THE REST OF THE RESULTS ###########

f11_UGn <- grepl("poptrue", row.names(DPpresults_U_G_n_normal_base))
m11_UGn <- DPpresults_U_G_n_normal_base[f11_UGn,]
mu_DPp_UGn_normal_base<- m11_UGn$`50%`
LB_mu_DPp_UGn_normal_base <- m11_UGn$`2.5%`
UB_mu_DPp_UGn_normal_base <- m11_UGn$`97.5%`
Rhat_muDPp_UGn_normal_base <- m11_UGn$Rhat
precDPp_mu_UGn_normal_base <- UB_mu_DPp_UGn_normal_base - LB_mu_DPp_UGn_normal_base


f17_UGn <- grepl("delta_new", row.names(DPpresults_U_G_n_normal_base))
m17_UGn <- DPpresults_U_G_n_normal_base[f17_UGn,]
delta_new_DPp_UGn_normal_base<- m17_UGn$`50%`
LB_delta_new_DPp_UGn_normal_base <- m17_UGn$`2.5%`
UB_delta_new_DPp_UGn_normal_base <- m17_UGn$`97.5%`
Rhat_delta_newDPp_UGn_normal_base <- m17_UGn$Rhat
precDPp_delta_new_UGn_normal_base <- UB_delta_new_DPp_UGn_normal_base - LB_delta_new_DPp_UGn_normal_base


f22_UGn <- grepl("var.true", row.names(DPpresults_U_G_n_normal_base))
m22_UGn <- DPpresults_U_G_n_normal_base[f22_UGn,]
tau2_DPp_UGn_normal_base <- m22_UGn$`50%`
LB_tau2_DPp_UGn_normal_base <- m22_UGn$`2.5%`
UB_tau2_DPp_UGn_normal_base <- m22_UGn$`97.5%`
Rhat_tauDPp_UGn_normal_base <- m22_UGn$Rhat
precDPp_tau2_UGn_normal_base <- UB_tau2_DPp_UGn_normal_base -  LB_tau2_DPp_UGn_normal_base


f33_UGn <- grepl("basemu", row.names(DPpresults_U_G_n_normal_base))
m33_UGn <- DPpresults_U_G_n_normal_base[f33_UGn,]
base_mu_DPp_UGn_normal_base <- m33_UGn$`50%`
LB_base_mu_DPp_UGn_normal_base <- m33_UGn$`2.5%`
UB_base_mu_DPp_UGn_normal_base <- m33_UGn$`97.5%`
Rhat_base_mu_DPp_UGn_normal_base <- m33_UGn$Rhat
prec_base_mu_DPp_UGn_normal_base <- UB_base_mu_DPp_UGn_normal_base - LB_base_mu_DPp_UGn_normal_base

f44_UGn <- grepl("basetau", row.names(DPpresults_U_G_n_normal_base))
m44_UGn <- DPpresults_U_G_n_normal_base[f44_UGn,]
base_tau_DPp_UGn_normal_base <- m44_UGn$`50%`
base_tau2_DPp_UGn_normal_base <- (m44_UGn$`50%`)^2
LB_base_tau2_DPp_UGn_normal_base <- (m44_UGn$`2.5%`)^2
UB_base_tau2_DPp_UGn_normal_base <- (m44_UGn$`97.5%`)^2
Rhat_base_tau_DPp_UGn_normal_base <- (m44_UGn$Rhat)^2
prec_base_tau2_DPp_UGn_normal_base <- UB_base_tau2_DPp_UGn_normal_base - LB_base_tau2_DPp_UGn_normal_base

f55_UGn <- grepl("alpha", row.names(DPpresults_U_G_n_normal_base))
m55_UGn <- DPpresults_U_G_n_normal_base[f55_UGn,]
alpha_DPp_UGn_normal_base <- m55_UGn$`50%`
LB_alpha_DPp_UGn_normal_base <- m55_UGn$`2.5%`
UB_alpha_DPp_UGn_normal_base <- m55_UGn$`97.5%`
Rhat_alpha_DPp_UGn_normal_base <- m55_UGn$Rhat
prec_alpha_DPp_UGn_normal_base <- UB_alpha_DPp_UGn_normal_base - LB_alpha_DPp_UGn_normal_base

fcl_UGn <- grepl("K", row.names(DPpresults_U_G_n_normal_base))   
mcl_UGn <- DPpresults_U_G_n_normal_base[fcl_UGn,]
median_K_DPp_UGn_normal_base <- mcl_UGn$`50%` 
LB_K_DPp_UGn_normal_base <- mcl_UGn$`2.5%`
UB_K_DPp_UGn_normal_base <- mcl_UGn$`97.5%`

fp_UGn <- grepl("p", row.names(DPpresults_U_G_n_normal_base))
mp_UGn <- DPpresults_U_G_n_normal_base[fp_UGn,]
mp_UGn <- mp_UGn[!grepl("pop", row.names(mp_UGn)),]
mp_UGn <- mp_UGn[!grepl("alpha", row.names(mp_UGn)),]
pi_DPp_UGn_normal_base <- mp_UGn$mean


listaDPp_UGn_normal_base <- cbind.data.frame(rel_effDPp_UGn_normal_base,sd_rel_effDPp_UGn_normal_base, LB_rel_effDPp_UGn_normal_base,
                                             UB_rel_effDPp_UGn_normal_base,Rhat_deltaDPp_UGn_normal_base)


numclus_DPp_UGn_normal_base <- unlist(median_K_DPp_UGn_normal_base)
LB_K_DPp_UGn_normal_base <- unlist(LB_K_DPp_UGn_normal_base)
UB_K_DPp_UGn_normal_base <- unlist(UB_K_DPp_UGn_normal_base)

resDPp_UGn_normal_base <- cbind.data.frame(base_mu_DPp_UGn_normal_base,LB_base_mu_DPp_UGn_normal_base,
                                           UB_base_mu_DPp_UGn_normal_base,base_tau_DPp_UGn_normal_base,base_tau2_DPp_UGn_normal_base,
                                           LB_base_tau2_DPp_UGn_normal_base,
                                           UB_base_tau2_DPp_UGn_normal_base,
                                           mu_DPp_UGn_normal_base,
                                           LB_mu_DPp_UGn_normal_base, 
                                           UB_mu_DPp_UGn_normal_base, 
                                           delta_new_DPp_UGn_normal_base,
                                           LB_delta_new_DPp_UGn_normal_base, 
                                           UB_delta_new_DPp_UGn_normal_base,
                                           tau2_DPp_UGn_normal_base, LB_tau2_DPp_UGn_normal_base, UB_tau2_DPp_UGn_normal_base,
                                           precDPp_mu_UGn_normal_base,precDPp_delta_new_UGn_normal_base,
                                           precDPp_tau2_UGn_normal_base, prec_base_mu_DPp_UGn_normal_base, prec_base_tau2_DPp_UGn_normal_base,
                                           prec_alpha_DPp_UGn_normal_base,  alpha_DPp_UGn_normal_base , LB_alpha_DPp_UGn_normal_base, UB_alpha_DPp_UGn_normal_base,
                                           Rhat_muDPp_UGn_normal_base,Rhat_delta_newDPp_UGn_normal_base,
                                           Rhat_tauDPp_UGn_normal_base, Rhat_base_mu_DPp_UGn_normal_base,
                                           Rhat_base_tau_DPp_UGn_normal_base, Rhat_alpha_DPp_UGn_normal_base, 
                                           numclus_DPp_UGn_normal_base, LB_K_DPp_UGn_normal_base, UB_K_DPp_UGn_normal_base,
                                           Pr_mu_DPp_UGn_normal_base, Pr_delta_new_DPp_UGn_normal_base)

extra_col = row.names(DPpresults_U_G_n_normal_base)
DPpresults_U_G_n_normal_base$extra_col = extra_col
prob_DPp_UGn_normal_base = round(max_values,2)

# 
# write.csv(resDPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\corrections\\2res_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
# write.csv(listaDPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\corrections\\2rel_eff_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
# write.csv(prob_DPp_UGn_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\corrections\\2max_prob_cluster_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
# write.csv(DPpresults_U_G_n_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_G_n_normal_base_58\\corrections\\2all_res_DPp_n_U_G_normal_base.csv",row.names=FALSE )  
# 
