## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

### Binomial-DPp-26-N(Unif/Unif) ###
if (requireNamespace("neojags", quietly = TRUE)){
  neojags::load.neojagsmodule()
} 
#> module neojags loaded
if (requireNamespace("neojags", quietly = TRUE)){
  library(rjags)
} 

## impliment Fernandez and Steel skew-t distribution for the base distribution 
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
    
    theta[k] ~ dfskew.t(xi, omega1, nu, shape)
  }

   # Priors
   
      xi ~ dnorm(0, .0001)
      omega1 <- 1 / omega_sqr
      omega_sqr <- omega * omega
      omega ~ dunif(0,10)
      shape ~ dnorm(0.1,0.04) ## 0.04 = 1/25 is the precision of normal
      nu ~ dexp(0.10) T(2.5, )
      ## shape ~ dgamma(0.5, 0.318) a prior they propose but it is similae to unif(0,5) we applied

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
  
}", file = "DPp_26_U_U_sk_t_base.txt")
modfile = 'DPp_26_U_U_sk_t_base.txt'



run.modelDPp26_U_U_sk_t_base = jags(
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
    "xi", ## mu of the Normal base  distribution
    "omega", ## tau of the Normal base distribution
    "shape",
    "nu",
    "poptrue",   ## overall mu
    "var.true", ## between-study variance
    "delta12", #### study-specific effects
    "K",       ### total number of clusters 
    "p",      ## weights of the process 
    "alpha",  ## concentration parameter
    "SC"   ,  ## probability of each cluster assignment
    "delta_new",
    "pred",
    "equalsres",
    "equalsmatrix"
  ),   
  
  n.chains = 2,
  n.iter = 50000,
  
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile
  
)

##### to check the convergence ###########
library(coda)
mcmc_obj <- as.mcmc(run.modelDPp26_U_U_sk_t_base)

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^theta\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^delta12\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters


par(mfrow = c(4, 3))  # Adjust margins if needed
traceplot(mcmc_obj[, c("xi", "omega", "alpha","nu", "delta_new", 
                       "poptrue", "shape",
                       "var.true")])

par(mfrow = c(1, 1))  # Reset layout

DPpresults_U_U_26_sk_t_base <- as.data.frame(run.modelDPp26_U_U_sk_t_base$BUGSoutput$summary) 

DPp_results_U_U_26_sk_t_base <- as.data.frame(run.modelDPp26_U_U_sk_t_base$BUGSoutput$summary) 

DPp_results_U_U_26_sk_t_base <-round(DPpresults_U_U_26_sk_t_base, digits=2) 

View(round(DPpresults_U_U_26_sk_t_base,2))
View(DPpresults_U_U_26_sk_t_base)


######## Pr(mu<0) #######
Pr_mu_DPp_U_U26_sk_t_base = mean(run.modelDPp26_U_U_sk_t_base$BUGSoutput$sims.matrix[  ,"poptrue"] < 0 )

######## Pr(mu_new<0) #######
Pr_delta_new_DPp_U_U26_sk_t_base = mean(run.modelDPp26_U_U_sk_t_base$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

### extract co-clustering probabilities
sims_equals <- run.modelDPp26_U_U_sk_t_base$BUGSoutput$sims.matrix
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



fdd <- grepl("delta12", row.names(DPpresults_U_U_26_sk_t_base))
mdd <- DPpresults_U_U_26_sk_t_base[fdd,]
rel_effDPp_U_U26_sk_t_base <- mdd$`50%`
LB_rel_effDPp_U_U26_sk_t_base <- mdd$`2.5%`
UB_rel_effDPp_U_U26_sk_t_base <- mdd$`97.5%`
sd_rel_effDPp_U_U26_sk_t_base <- mdd$sd

fdd1 <- grepl("theta", row.names(DPpresults_U_U_26_sk_t_base))
mdd1 <- DPpresults_U_U_26_sk_t_base[fdd1,]
theta_DPp_U_U26_sk_t_base <- mdd1$`50%`

plot(density(theta_DPp_U_U26_sk_t_base))

########### SAVE THE REST OF THE RESULTS ###########

f11_U_U26 <- grepl("poptrue", row.names(DPpresults_U_U_26_sk_t_base))
m11_U_U26 <- DPpresults_U_U_26_sk_t_base[f11_U_U26,]
mu_DPp_U_U26_sk_t_base<- m11_U_U26$`50%`
LB_mu_DPp_U_U26_sk_t_base <- m11_U_U26$`2.5%`
UB_mu_DPp_U_U26_sk_t_base <- m11_U_U26$`97.5%`
precDPp_mu_U_U26_sk_t_base <- UB_mu_DPp_U_U26_sk_t_base - LB_mu_DPp_U_U26_sk_t_base


f117_U_U26 <- grepl("delta_new", row.names(DPpresults_U_U_26_sk_t_base))
m117_U_U26 <- DPpresults_U_U_26_sk_t_base[f117_U_U26,]
delta_new_DPp_U_U26_sk_t_base<- m117_U_U26$`50%`
LB_delta_new_DPp_U_U26_sk_t_base <- m117_U_U26$`2.5%`
UB_delta_new_DPp_U_U26_sk_t_base <- m117_U_U26$`97.5%`
precDPp_delta_new_U_U26_sk_t_base <- UB_delta_new_DPp_U_U26_sk_t_base - LB_delta_new_DPp_U_U26_sk_t_base


f22_U_U26 <- grepl("var.true", row.names(DPpresults_U_U_26_sk_t_base))
m22_U_U26 <- DPpresults_U_U_26_sk_t_base[f22_U_U26,]
tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`50%`
LB_tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`2.5%`
UB_tau2_DPp_U_U26_sk_t_base <- m22_U_U26$`97.5%`
precDPp_tau2_U_U26_sk_t_base <- UB_tau2_DPp_U_U26_sk_t_base -  LB_tau2_DPp_U_U26_sk_t_base


f33_U_U26 <- grepl("xi", row.names(DPpresults_U_U_26_sk_t_base))
m33_U_U26 <- DPpresults_U_U_26_sk_t_base[f33_U_U26,]
xi_DPp_U_U26_sk_t_base <- m33_U_U26$`50%`
LB_xi_DPp_U_U26_sk_t_base <- m33_U_U26$`2.5%`
UB_xi_DPp_U_U26_sk_t_base <- m33_U_U26$`97.5%`
prec_xi_DPp_U_U26_sk_t_base <- UB_xi_DPp_U_U26_sk_t_base - LB_xi_DPp_U_U26_sk_t_base

f44_U_U26 <- grepl("omega", row.names(DPpresults_U_U_26_sk_t_base))
m44_U_U26 <- DPpresults_U_U_26_sk_t_base[f44_U_U26,]
omega_DPp_U_U26_sk_t_base <- m44_U_U26$`50%`
omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`50%`)^2
LB_omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`2.5%`)^2
UB_omega2_DPp_U_U26_sk_t_base <- (m44_U_U26$`97.5%`)^2
prec_omega2_DPp_U_U26_sk_t_base <- UB_omega2_DPp_U_U26_sk_t_base - LB_omega2_DPp_U_U26_sk_t_base

f55_U_U26 <- grepl("shape", row.names(DPpresults_U_U_26_sk_t_base))
m55_U_U26 <- DPpresults_U_U_26_sk_t_base[f55_U_U26,]
shape_DPp_U_U26_sk_t_base <- m55_U_U26$`50%`
LB_shape_DPp_U_U26_sk_t_base <- m55_U_U26$`2.5%`
UB_shape_DPp_U_U26_sk_t_base <- m55_U_U26$`97.5%`
prec_shape_DPp_U_U26_sk_t_base <- UB_shape_DPp_U_U26_sk_t_base - LB_shape_DPp_U_U26_sk_t_base

f551_U_U26 <- grepl("nu", row.names(DPpresults_U_U_26_sk_t_base))
m551_U_U26 <- DPpresults_U_U_26_sk_t_base[f551_U_U26,]
df_DPp_U_U26_sk_t_base <- m551_U_U26$`50%`
LB_df_DPp_U_U26_sk_t_base <- m551_U_U26$`2.5%`
UB_df_DPp_U_U26_sk_t_base <- m551_U_U26$`97.5%`
prec_df_DPp_U_U26_sk_t_base <- UB_df_DPp_U_U26_sk_t_base - LB_df_DPp_U_U26_sk_t_base

f55_U_U26 <- grepl("alpha", row.names(DPpresults_U_U_26_sk_t_base))
m55_U_U26 <- DPpresults_U_U_26_sk_t_base[f55_U_U26,]
alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`50%`
LB_alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`2.5%`
UB_alpha_DPp_U_U26_sk_t_base <- m55_U_U26$`97.5%`
prec_alpha_DPp_U_U26_sk_t_base <- UB_alpha_DPp_U_U26_sk_t_base - LB_alpha_DPp_U_U26_sk_t_base

fcl_U_U26 <- grepl("K", row.names(DPpresults_U_U_26_sk_t_base))   
mcl_U_U26 <- DPpresults_U_U_26_sk_t_base[fcl_U_U26,]
median_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`50%` 
LB_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`2.5%`
UB_K_DPp_U_U26_sk_t_base <- mcl_U_U26$`97.5%`

fp_U_U26 <- grepl("p", row.names(DPpresults_U_U_26_sk_t_base))
mp_U_U26 <- DPpresults_U_U_26_sk_t_base[fp_U_U26,]
mp_U_U26 <- mp_U_U26[!grepl("pop", row.names(mp_U_U26)),]
mp_U_U26 <- mp_U_U26[!grepl("alpha", row.names(mp_U_U26)),]
pi_DPp_U_U26_sk_t_base <- mp_U_U26$mean


listaDPp_U_U26_sk_t_base <- cbind.data.frame(rel_effDPp_U_U26_sk_t_base,sd_rel_effDPp_U_U26_sk_t_base, LB_rel_effDPp_U_U26_sk_t_base,
                                             UB_rel_effDPp_U_U26_sk_t_base)


numclus_DPp_U_U26_sk_t_base <- unlist(median_K_DPp_U_U26_sk_t_base)
LB_K_DPp_U_U26_sk_t_base <- unlist(LB_K_DPp_U_U26_sk_t_base)
UB_K_DPp_U_U26_sk_t_base <- unlist(UB_K_DPp_U_U26_sk_t_base)

resDPp_U_U26_sk_t_base <- cbind.data.frame(
  xi_DPp_U_U26_sk_t_base,LB_xi_DPp_U_U26_sk_t_base,UB_xi_DPp_U_U26_sk_t_base,
  omega_DPp_U_U26_sk_t_base,
  omega2_DPp_U_U26_sk_t_base,LB_omega2_DPp_U_U26_sk_t_base,UB_omega2_DPp_U_U26_sk_t_base,
  shape_DPp_U_U26_sk_t_base, LB_shape_DPp_U_U26_sk_t_base, UB_shape_DPp_U_U26_sk_t_base,
  df_DPp_U_U26_sk_t_base, LB_df_DPp_U_U26_sk_t_base, UB_df_DPp_U_U26_sk_t_base,
  mu_DPp_U_U26_sk_t_base,LB_mu_DPp_U_U26_sk_t_base, UB_mu_DPp_U_U26_sk_t_base,
  delta_new_DPp_U_U26_sk_t_base,LB_delta_new_DPp_U_U26_sk_t_base, UB_delta_new_DPp_U_U26_sk_t_base,
  tau2_DPp_U_U26_sk_t_base, LB_tau2_DPp_U_U26_sk_t_base, UB_tau2_DPp_U_U26_sk_t_base,
  prec_df_DPp_U_U26_sk_t_base,
  precDPp_mu_U_U26_sk_t_base, precDPp_delta_new_U_U26_sk_t_base,
  precDPp_tau2_U_U26_sk_t_base, prec_xi_DPp_U_U26_sk_t_base, prec_omega2_DPp_U_U26_sk_t_base,
  prec_alpha_DPp_U_U26_sk_t_base,  alpha_DPp_U_U26_sk_t_base , LB_alpha_DPp_U_U26_sk_t_base, UB_alpha_DPp_U_U26_sk_t_base,
  numclus_DPp_U_U26_sk_t_base, LB_K_DPp_U_U26_sk_t_base, UB_K_DPp_U_U26_sk_t_base,
  Pr_mu_DPp_U_U26_sk_t_base, Pr_delta_new_DPp_U_U26_sk_t_base)

extra_col = row.names(DPpresults_U_U_26_sk_t_base)
DPpresults_U_U_26_sk_t_base$extra_col = extra_col
prob_DPp_U_U26_sk_t_base = round(max_values,2)

# 
# write.csv(resDPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\corrections\\2res_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
# write.csv(listaDPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\corrections\\2rel_eff_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
# write.csv(prob_DPp_U_U26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\corrections\\2max_prob_cluster_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  
# write.csv(DPpresults_U_U_26_sk_t_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\DPp_U_U_26_st_base_pred_58\\corrections\\2all_res_DPp_26_U_U_sk_t_base_pred.csv",row.names=FALSE )  

