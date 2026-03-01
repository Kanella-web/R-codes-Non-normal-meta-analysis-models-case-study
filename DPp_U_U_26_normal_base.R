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

### Binomial-DPp-26-N(U/Unif) ###
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

  # Constructive DP prior
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
    
    theta[k] ~ dnorm(basemu, basetau1)  
  }

   # Priors
   
   basetau1 <- 1 / basetau_sqr  
   basetau_sqr <- basetau * basetau
   basetau ~ dunif(0,10)  

   #basemu ~ dnorm(0, 0.0001)   
   ## or skeptical prior 
   basemu ~ dnorm(0, prec_mu)
   prec_mu <- 1 / (0.255)^2
   
  # concentration parameter prior
  alpha ~ dunif(0.3, 5)
  
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

  # Summary statistics
  for (i in 1:ns) {
    for (j in 1:N) {
      SC[i, j] <- equals(j, Z[i])
    }
  }

  # Total clusters K
  for (j in 1:N) {
    cl[j] <- step(sum(SC[, j]) - 1)
  }
  K <- sum(cl[])
  
  ## co-clustering probabilities
  for (i in 1:ns) {
    for (j in 1:ns) {
     equalsmatrix[i, j] <- equals(Z[i], Z[j])
     
    }
  equalsres[i] <- sum(equalsmatrix[i, ])
  
}

  pred ~ dcat(p[1:N]) ##First randomly assigning the new study to one of the mixture components, according to the estimated weights p.
  delta_new <- theta[pred] ##Then drawing the new study's effect size from the random effects distribution.
  
}", file = "DPp_26_U_U_normal_base2.txt")
modfile = 'DPp_26_U_U_normal_base2.txt'

run.modelDPp26_U_U_normal_base = jags(
  dati1  <-list(ns = nrow(pre_term_data_58),
                y1 = pre_term_data_58$mean_FT,
                y2 = pre_term_data_58$mean_EPT.VPT,
                sd1 = pre_term_data_58$sd_FT,
                sd2 = pre_term_data_58$sd_EPT.VPT,
                n1 = pre_term_data_58$n_FT,
                n2 = pre_term_data_58$n_EPT.VPT,
                N = 26
  ) ,  
  inits = NULL,
  parameters.to.save = c(
    "theta",
    "basemu", ## mu of the Normal base  distribution
    "basetau", ## tau of the Normal base distribution
    "poptrue",   ## overall SMD
    "var.true", ## between-study variance
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "SC",     ## probability of each cluster assignment
    "delta_new",
    "pred",
    "Z",
    "equalsres",
    "equalsmatrix" ## co-clustering probabilities
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)


DPpresults_U_U_26_normal_base <- as.data.frame(run.modelDPp26_U_U_normal_base$BUGSoutput$summary) 

DPp_results_U_U_26_normal_base <- as.data.frame(run.modelDPp26_U_U_normal_base$BUGSoutput$summary) 

DPp_results_U_U_26_normal_base <-round(DPpresults_U_U_26_normal_base, digits=2) 

View(DPp_results_U_U_26_normal_base)


### extract co-clustering probabilities
sims_equals <- run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix
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

#write.csv(distance.studies, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\corrections\\co-clust.csv",row.names=TRUE )

##define the columns and the rows accordingly
pre_term_data_58$Study
study_names <- pre_term_data_58$Study
rownames(distance.studies) <- study_names
colnames(distance.studies) <- study_names

## heatmap
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




######## Pr(mu<0) #######
Pr_mu_DPp_U_U26_normal_base = mean(run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )

######## Pr(mu_new<0) #######
Pr_delta_new_DPp_U_U26_normal_base = mean(run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

########## CLUSTER ASSIGNEMENT ##########
DPpresults_U_U_26_normal_base$ind <- row.names(DPpresults_U_U_26_normal_base)


################
fdd1 <- grepl("theta", row.names(DPpresults_U_U_26_normal_base))
mdd1 <- DPpresults_U_U_26_normal_base[fdd1,]
thetaDPp_U_U26_normal_base <- mdd1$`50%`

#### DENSITY PLOT OF THE RELATIVE EFFECTS #######
plot(density(thetaDPp_U_U26_normal_base))


fdd <- grepl("delta12", row.names(DPpresults_U_U_26_normal_base))
mdd <- DPpresults_U_U_26_normal_base[fdd,]
rel_effDPp_U_U26_normal_base <- mdd$`50%`
LB_rel_effDPp_U_U26_normal_base <- mdd$`2.5%`
UB_rel_effDPp_U_U26_normal_base <- mdd$`97.5%`
sd_rel_effDPp_U_U26_normal_base <- mdd$sd
Rhat_deltaDPp_U_U26_normal_base <- mdd$Rhat


#### DENSITY PLOT OF THE RELATIVE EFFECTS #######
plot(density(rel_effDPp_U_U26_normal_base))

########### SAVE THE REST OF THE RESULTS ###########

f11_U_U26 <- grepl("poptrue", row.names(DPpresults_U_U_26_normal_base))
m11_U_U26 <- DPpresults_U_U_26_normal_base[f11_U_U26,]
mu_DPp_U_U26_normal_base<- m11_U_U26$`50%`
LB_mu_DPp_U_U26_normal_base <- m11_U_U26$`2.5%`
UB_mu_DPp_U_U26_normal_base <- m11_U_U26$`97.5%`
Rhat_muDPp_U_U26_normal_base <- m11_U_U26$Rhat
precDPp_mu_U_U26_normal_base <- UB_mu_DPp_U_U26_normal_base - LB_mu_DPp_U_U26_normal_base

f17_U_U26 <- grepl("delta_new", row.names(DPpresults_U_U_26_normal_base))
m17_U_U26 <- DPpresults_U_U_26_normal_base[f17_U_U26,]
delta_new_DPp_U_U26_normal_base<- m17_U_U26$`50%`
LB_delta_new_DPp_U_U26_normal_base <- m17_U_U26$`2.5%`
UB_delta_new_DPp_U_U26_normal_base <- m17_U_U26$`97.5%`
Rhat_delta_newDPp_U_U26_normal_base <- m17_U_U26$Rhat
precDPp_delta_new_U_U26_normal_base <- UB_delta_new_DPp_U_U26_normal_base - LB_delta_new_DPp_U_U26_normal_base


f22_U_U26 <- grepl("var.true", row.names(DPpresults_U_U_26_normal_base))
m22_U_U26 <- DPpresults_U_U_26_normal_base[f22_U_U26,]
tau2_DPp_U_U26_normal_base <- m22_U_U26$`50%`
LB_tau2_DPp_U_U26_normal_base <- m22_U_U26$`2.5%`
UB_tau2_DPp_U_U26_normal_base <- m22_U_U26$`97.5%`
Rhat_tauDPp_U_U26_normal_base <- m22_U_U26$Rhat
precDPp_tau2_U_U26_normal_base <- UB_tau2_DPp_U_U26_normal_base -  LB_tau2_DPp_U_U26_normal_base


f33_U_U26 <- grepl("basemu", row.names(DPpresults_U_U_26_normal_base))
m33_U_U26 <- DPpresults_U_U_26_normal_base[f33_U_U26,]
base_mu_DPp_U_U26_normal_base <- m33_U_U26$`50%`
LB_base_mu_DPp_U_U26_normal_base <- m33_U_U26$`2.5%`
UB_base_mu_DPp_U_U26_normal_base <- m33_U_U26$`97.5%`
Rhat_base_mu_DPp_U_U26_normal_base <- m33_U_U26$Rhat
prec_base_mu_DPp_U_U26_normal_base <- UB_base_mu_DPp_U_U26_normal_base - LB_base_mu_DPp_U_U26_normal_base

f44_U_U26 <- grepl("basetau", row.names(DPpresults_U_U_26_normal_base))
m44_U_U26 <- DPpresults_U_U_26_normal_base[f44_U_U26,]
base_tau_DPp_U_U26_normal_base <- m44_U_U26$`50%`
base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`50%`)^2
LB_base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`2.5%`)^2
UB_base_tau2_DPp_U_U26_normal_base <- (m44_U_U26$`97.5%`)^2
Rhat_base_tau_DPp_U_U26_normal_base <- (m44_U_U26$Rhat)^2
prec_base_tau2_DPp_U_U26_normal_base <- UB_base_tau2_DPp_U_U26_normal_base - LB_base_tau2_DPp_U_U26_normal_base

f55_U_U26 <- grepl("alpha", row.names(DPpresults_U_U_26_normal_base))
m55_U_U26 <- DPpresults_U_U_26_normal_base[f55_U_U26,]
alpha_DPp_U_U26_normal_base <- m55_U_U26$`50%`
LB_alpha_DPp_U_U26_normal_base <- m55_U_U26$`2.5%`
UB_alpha_DPp_U_U26_normal_base <- m55_U_U26$`97.5%`
Rhat_alpha_DPp_U_U26_normal_base <- m55_U_U26$Rhat
prec_alpha_DPp_U_U26_normal_base <- UB_alpha_DPp_U_U26_normal_base - LB_alpha_DPp_U_U26_normal_base

fcl_U_U26 <- grepl("K", row.names(DPpresults_U_U_26_normal_base))   
mcl_U_U26 <- DPpresults_U_U_26_normal_base[fcl_U_U26,]
median_K_DPp_U_U26_normal_base <- mcl_U_U26$`50%` 
LB_K_DPp_U_U26_normal_base <- mcl_U_U26$`2.5%`
UB_K_DPp_U_U26_normal_base <- mcl_U_U26$`97.5%`

fp_U_U26 <- grepl("p", row.names(DPpresults_U_U_26_normal_base))
mp_U_U26 <- DPpresults_U_U_26_normal_base[fp_U_U26,]
mp_U_U26 <- mp_U_U26[!grepl("pop", row.names(mp_U_U26)),]
mp_U_U26 <- mp_U_U26[!grepl("alpha", row.names(mp_U_U26)),]
pi_DPp_U_U26_normal_base <- mp_U_U26$mean


listaDPp_U_U26_normal_base <- cbind.data.frame(rel_effDPp_U_U26_normal_base,sd_rel_effDPp_U_U26_normal_base, LB_rel_effDPp_U_U26_normal_base,
                                               UB_rel_effDPp_U_U26_normal_base,Rhat_deltaDPp_U_U26_normal_base)


numclus_DPp_U_U26_normal_base <- unlist(median_K_DPp_U_U26_normal_base)
LB_K_DPp_U_U26_normal_base <- unlist(LB_K_DPp_U_U26_normal_base)
UB_K_DPp_U_U26_normal_base <- unlist(UB_K_DPp_U_U26_normal_base)

resDPp_U_U26_normal_base <- cbind.data.frame(base_mu_DPp_U_U26_normal_base,LB_base_mu_DPp_U_U26_normal_base,
                                             UB_base_mu_DPp_U_U26_normal_base,base_tau_DPp_U_U26_normal_base,base_tau2_DPp_U_U26_normal_base,
                                             LB_base_tau2_DPp_U_U26_normal_base,
                                             UB_base_tau2_DPp_U_U26_normal_base,
                                             mu_DPp_U_U26_normal_base,
                                             LB_mu_DPp_U_U26_normal_base, 
                                             UB_mu_DPp_U_U26_normal_base, tau2_DPp_U_U26_normal_base, LB_tau2_DPp_U_U26_normal_base, UB_tau2_DPp_U_U26_normal_base,
                                             delta_new_DPp_U_U26_normal_base, LB_delta_new_DPp_U_U26_normal_base, UB_delta_new_DPp_U_U26_normal_base,
                                             
                                             precDPp_mu_U_U26_normal_base,precDPp_delta_new_U_U26_normal_base,
                                             precDPp_tau2_U_U26_normal_base, prec_base_mu_DPp_U_U26_normal_base, prec_base_tau2_DPp_U_U26_normal_base,
                                             prec_alpha_DPp_U_U26_normal_base,  alpha_DPp_U_U26_normal_base , LB_alpha_DPp_U_U26_normal_base, UB_alpha_DPp_U_U26_normal_base,
                                             Rhat_muDPp_U_U26_normal_base, Rhat_delta_newDPp_U_U26_normal_base,
                                             Rhat_tauDPp_U_U26_normal_base, Rhat_base_mu_DPp_U_U26_normal_base,
                                             Rhat_base_tau_DPp_U_U26_normal_base, Rhat_alpha_DPp_U_U26_normal_base, 
                                             numclus_DPp_U_U26_normal_base, LB_K_DPp_U_U26_normal_base, UB_K_DPp_U_U26_normal_base,
                                             Pr_mu_DPp_U_U26_normal_base, Pr_delta_new_DPp_U_U26_normal_base)

extra_col = row.names(DPpresults_U_U_26_normal_base)
DPpresults_U_U_26_normal_base$extra_col = extra_col
prob_DPp_U_U26_normal_base = round(max_values,2)

# write.csv(resDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\corrections\\2res_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
# write.csv(listaDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\corrections\\2rel_eff_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
# write.csv(prob_DPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\corrections\\2max_prob_cluster_DPp_26_U_U_normal_base.csv",row.names=FALSE )  
#write.csv(DPpresults_U_U_26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_58_pred\\corrections\\2all_res_DPp_26_U_U_normal_base.csv",row.names=TRUE )  

