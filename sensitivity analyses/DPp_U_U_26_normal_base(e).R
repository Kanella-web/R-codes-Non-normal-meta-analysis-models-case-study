## load the data
pre_term_data_61 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data61.csv")
#pre_term_data_61 = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\pre_term_data61.csv")
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
#
run.modelDPp26_U_U_normal_base = jags(
  dati1  <-list(ns = nrow(pre_term_data_61),
                y1 = pre_term_data_61$mean_FT,
                y2 = pre_term_data_61$mean_EPT.VPT,
                sd1 = pre_term_data_61$sd_FT,
                sd2 = pre_term_data_61$sd_EPT.VPT,
                n1 = pre_term_data_61$n_FT,
                n2 = pre_term_data_61$n_EPT.VPT,
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
ns <- nrow(pre_term_data_61)
niter <- nrow(equals)
equals_array <- array(
  equals,
  dim = c(niter, ns, ns)
)
dim(equals_array)

## create the 3dimentional matrix
distance.studies <- apply(equals_array, c(2, 3), mean)

#write.csv(distance.studies, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data61\\DPp_U_U_61_pred\\corrections\\co-clust.csv",row.names=TRUE )

##define the columns and the rows accordingly
pre_term_data_61$Study
study_names <- pre_term_data_61$Study
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


write.csv(resDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\res_DPp_26_U_U_normal_base_61_studies.csv",row.names=TRUE )
write.csv(listaDPp_U_U26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\rel_eff_DPp_26_U_U_normal_base_61_studies.csv.csv",row.names=TRUE )
write.csv(DPpresults_U_U_26_normal_base, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\all_res_DPp_26_U_U_normal_base_61_studies.csv",row.names=TRUE )


## CLUSTER MEMBESHIPS
## extract indicator labels
sims <- run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix
Z_samps <- sims[, grep("^Z\\[", colnames(sims))]

## create the co-clustering matrix
ns <- ncol(Z_samps)

Pmat <- matrix(0, ns, ns)
for (i in 1:ns) {
  for (j in 1:ns) {
    Pmat[i, j] <- mean(Z_samps[, i] == Z_samps[, j])
  }
}

ns <- nrow(Pmat)

# find the maximum partner for each study (excluding itself)
max_partner <- integer(ns)
max_prob    <- numeric(ns)

for (i in 1:ns) {
  probs <- Pmat[i, ]
  probs[i] <- -Inf        # exclude self (diagonal = 1)
  
  max_partner[i] <- which.max(probs)
  max_prob[i]    <- max(probs)
}

pre_term_data_61$max_partner <- max_partner
pre_term_data_61$max_prob <- max_prob
###Ensures no study selected itself as partner
# should return FALSE
any(max_partner == seq_len(ns))


##cluster membership
library(igraph)

adj_max <- matrix(0, ns, ns)

for (i in 1:ns) {
  j <- max_partner[i]
  adj_max[i, j] <- 1
  adj_max[j, i] <- 1
}

g_max <- graph_from_adjacency_matrix(adj_max, mode = "undirected", diag = FALSE)
plot(g_max)

clusters_max <- components(g_max)$membership
plot(
  g_max,
  vertex.color = clusters_max,
  vertex.size  = 6,
  vertex.label = NA
)

clusters_max
pre_term_data_61$cluster_maxP <- clusters_max

### extract co-clustering probabilities from the model as before
sims_equals <- run.modelDPp26_U_U_normal_base$BUGSoutput$sims.matrix
equals <- sims_equals[, grep("^equalsmatrix\\[", colnames(sims_equals))]

## create the 3dimentional array
ns <- nrow(pre_term_data_61)
niter <- nrow(equals)
equals_array <- array(
  equals,
  dim = c(niter, ns, ns)
)
dim(equals_array)

## create the 3dimentional matrix
distance.studies <- apply(equals_array, c(2, 3), mean)

#write.csv(distance.studies, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data61\\DPp_U_U_61_pred\\corrections\\co-clust.csv",row.names=TRUE )

##define the columns and the rows accordingly
pre_term_data_61$Study
study_names <- pre_term_data_61$Study
rownames(distance.studies)= paste(pre_term_data_61$cluster_maxP, 
                                  pre_term_data_61$Study,
                                  sep = " | ")
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



################### FOREST PLTOS ###############
pre_term_data_61$IQ_type <- ifelse(
  pre_term_data_61$IQ_type == 1, "short", "long"
)


es = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\res_DPp_26_U_U_normal_base_61_studies.csv")
rel_eff = read.csv("C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\sens_61_studies\\rel_eff_DPp_26_U_U_normal_base_61_studies.csv.csv")

## ===============================
## COLORS & SYMBOLS

numbers <- clusters_max
clusters = numbers

## ---- explicit colors ----
color_vector <- rep("#737373", length(numbers))   # grey default
color_vector[numbers == 1] <-     "black"      
color_vector[numbers == 3] <-  "#b71c1c"           
color_vector[numbers == 2] <- "#1f77b4"      
color_vector[numbers == 4] <- "orange"      


## ---- symbols (unchanged logic) ----
unique_numbers <- sort(unique(numbers))
replacement_values <- c(17,16,15,13,12)

replacement_map <- setNames(replacement_values, unique_numbers)
new_numbers <- replacement_map[as.character(numbers)]

pr <- round(max_prob,3)

## ===============================
## ADDITIONAL STUDY INFO
## ===============================

PMB      <- pre_term_data_61$Birth_level
matched <- ifelse(pre_term_data_61$Matched == 1, "yes", "no")
mean_ga <- pre_term_data_61$Mean..ga
type_IQ <- pre_term_data_61$IQ_type

add_data <- cbind.data.frame(
  PMB,
  mean_ga,
  type_IQ,
  matched,
  clusters_max,
  pr
)


## ===============================
## COMBINE EVERYTHING
## ===============================

sorted_data <- cbind.data.frame(
  rel_eff,
  add_data,
  color_vector,
  new_numbers,
  Study = pre_term_data_61$Study
)


## ===============================
## 🔥 CORRECT SORTING 🔥
## ===============================

# Explicit color order (VERY IMPORTANT)
color_levels <- c(  "#737373","orange", "#1f77b4","black", "#b71c1c")

sorted_data$color_id <- match(sorted_data$color_vector, color_levels)

sorted_data$new_numbers <- as.numeric(sorted_data$new_numbers)

pmb_levels <- c("EPT", "VPT")
sorted_data$PMB_id <- match(sorted_data$PMB, pmb_levels)

sorted_data <- sorted_data[
  order(
    sorted_data$color_id,
    sorted_data$new_numbers,
    sorted_data$PMB_id,
    -as.numeric(sorted_data$mean_ga)
  ),
]



## ===============================
## EXTRACT SORTED OBJECTS
## ===============================

sorted_rel_eff <- sorted_data$rel_effDPp_U_U26_normal_base
sorted_LB      <- sorted_data$LB_rel_effDPp_U_U26_normal_base
sorted_UB      <- sorted_data$UB_rel_effDPp_U_U26_normal_base
sorted_authors <- sorted_data$Study
sorted_color   <- sorted_data$color_vector
sorted_pch     <- sorted_data$new_numbers

sorted_add_data <- sorted_data[, c(
  "PMB",
  "mean_ga",
  "type_IQ",
  "matched",
  "clusters_max",
  "pr"
)]


## ===============================
## FOREST PLOT
## ===============================

library(metafor)

forest(
  x       = sorted_rel_eff,
  ci.lb  = sorted_LB,
  ci.ub  = sorted_UB,
  slab   = sorted_authors,
  xlim   = c(-7, 3),
  ylim   = c(-4, 65),
  psize  = 2,
  cex    = 0.6,
  lwd    = 1.4,
  col    = sorted_color,
  pch    = sorted_pch,
  ilab   = sorted_add_data,
  ilab.xpos = c(-5.2, -4.8, -4.3, -3.5, 0.9, 1.5),
  refline = TRUE
)


clusters_max
cluster_elements = split(seq_along(clusters_max), clusters_max)


fdd <- grepl("delta12", row.names(DPpresults_U_U_26_normal_base))
mdd <- DPpresults_U_U_26_normal_base[fdd,]
rel_effDPp_U_U26_normal_base <- mdd$`50%`
LB_rel_effDPp_U_U26_normal_base <- mdd$`2.5%`
UB_rel_effDPp_U_U26_normal_base <- mdd$`97.5%`
sd_rel_effDPp_U_U26_normal_base <- mdd$sd
Rhat_deltaDPp_U_U26_normal_base <- mdd$Rhat

cl1 = rel_effDPp_U_U26_normal_base[cluster_elements[[1]]]
cl2 = rel_effDPp_U_U26_normal_base[cluster_elements[[2]]]
cl3 = rel_effDPp_U_U26_normal_base[cluster_elements[[3]]]

mean_cl1 = mean(cl1)
mean_cl2 = mean(cl2)
mean_cl3 = mean(cl3)

cluster_means = cbind.data.frame(mean_cl3, mean_cl2, mean_cl1)

abline(v=cluster_means , col = c( "#b71c1c","#1f77b4","black" ), lwd=1.5)
abline(v=0, col = "black",lty= 2, lwd=1.5)


addpoly(x= es$mu_DPp_U_U26_normal_base, ci.lb = es$LB_mu_DPp_U_U26_normal_base ,
        ci.ub = es$UB_mu_DPp_U_U26_normal_base , rows=-1)
arrows(x0 = es$LB_delta_new_DPp_U_U26_normal_base,
       x1 = es$UB_delta_new_DPp_U_U26_normal_base,
       y0 = -2.5, y1 = -2.5,
       angle = 90, code = 3, length = 0.05, lty = 1, lwd = 1.5, col =  "#D2691E")
# if you want to add the point estimate of the PI 
#points(es$delta_new_DPp_U_U26_normal_base, -1.2, pch = 15, col = "black")
abline(h=0.3, lwd=0.1, col="black", lty=1)


