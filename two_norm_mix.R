## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")


# ## analysis only with recent studies ###
# #keep only studies with publication year of >= 2016
# pre_term_data_58$numbers <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# pre_term_data_58$publication_year <- as.numeric(gsub("[^0-9]", "", pre_term_data_58$Study))
# 
# pre_term_data_58 <- pre_term_data_58[
#   pre_term_data_58$publication_year >= 2016,
# ]

############ Mixture of two normal distributions model ############
library(R2jags)
set.seed(1508)
####normal-2normal-mixture model###
writeLines("
  model {
    for (i in 1:ns) { 
      
      m[i] ~ dnorm(0,.0001)
      
      # Mean outcome for group 1
      y1[i] ~ dnorm(mu1[i], prec1[i])
      mu1[i] <- m[i]
      
      # Mean outcome for group 2
      y2[i] ~ dnorm(mu2[i], prec2[i])
      mu2[i] <- m[i] + delta12[i]*pooled_sd[i]
      
      # Between-study standardized effect
      delta12[i] ~ dnorm(mu[z[i]], tau1[z[i]])  ### component distribution
      
      z[i] ~ dcat(p[])  # Mixture component indicator
      
      # Precision parameters for group 1 and 2 based on observed standard deviations
      prec1[i] <- n1[i]/ pow(sd1[i], 2)
      prec2[i] <- n2[i]/ pow(sd2[i], 2)

      # Calculate pooled standard deviation
      pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))
    
    }

   # Priors for mixing proportions
    for (i in 1:N) {

        alpha[i] <- 1  

    }
    ###### Assigning a Dirichlet prior for p seasures that sum[p]=1 ###

    p ~ ddirch(alpha[])

#### Priors

    mu[1] ~ dnorm(m1, t1)  ##### sharing info across component means
    mu[2] ~ dnorm(m1, t1) T(mu[1], )

    for (k in 1:N) {
    
      tau1[k] ~ dgamma(3,b)
      
      tau_k_sqr[k] <- 1/ tau1[k]
      
      tau_comp[k] <- sqrt(tau_k_sqr[k])

    }
    
    
    m1 ~ dnorm(0,0.0001) 
    
    ## or skeptical prior 
    # m1 ~ dnorm(0, prec_mu)
    # prec_mu <- 1 / (0.255)^2


    t1 <- 1/ (t*t)
    t ~ dnorm(0, 1)T(0,)
    
    b~dgamma(0.03, 0.03)
    
    # Calculate mean of the mixture distribution
    for (i in 1:N) {
    
      mean_mu[i] <- p[i] * mu[i]
      
    }
    
    total_mean <- sum(mean_mu[])  # Total mean of the mixture
    
 
   # Random effects distribution variance #
  
  
  for(i in 1:N){
  
    var1[i]<-p[i]*tau_k_sqr[i]
  
    sq[i]<-mu[i]-total_mean
    
    var2[i]<-p[i]*sq[i]*sq[i]   #### E[X2]
  }
  
  total_var <-sum(var1[])+sum(var2[]) ### E[X2] - E[X]2
   
   
  # Summary statistics #
  for(i in 1:ns){
   for (j in 1:N){
    SC[i,j] <- equals(j,z[i])
   }
  }
  # total clusters K#
  for (j in 1:N){
   cl[j] <-step(sum(SC[,j])-1)
   }
  K<-sum(cl[])
  
   pred ~ dcat(p[1:N])
   delta_new ~ dnorm(mu[pred],tau1[pred])
  
}",
           con = "finite_mixture_norm_dir_weights_tau_comp_Gamma.txt")

modfile3 = 'finite_mixture_norm_dir_weights_tau_comp_Gamma.txt'

library(R2jags)

dati1  <-list(ns = nrow(pre_term_data_58),
              y1 = pre_term_data_58$mean_FT,
              y2 = pre_term_data_58$mean_EPT.VPT,
              sd1 = pre_term_data_58$sd_FT,
              sd2 = pre_term_data_58$sd_EPT.VPT,
              n1 = pre_term_data_58$n_FT,
              n2 = pre_term_data_58$n_EPT.VPT,
              N = 2)


run.model_fin_mix_norm_dir_weigths_tau_comp = jags(
  data = dati1,
  inits = NULL,
  parameters.to.save = c(
    "mu",    ## mu of each mixture component
    "tau_comp",  
    "delta12", ## study-specific effects
    "tau_k_sqr", ## tau2 of each mixture component
    "p",
    "total_mean",
    "total_var",
    "SC",
    "z",     ## cluster assignment 
    "K" ,    ## total number of components (known)
    "pred",
    "delta_new"
  ),
  n.chains = 2,
  n.iter = 50000,
  
  n.burnin = 10000,
  DIC = T,
  model.file = modfile3
)

results_fin_mix_norm_dir_weigths = as.data.frame(run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$summary)
results_fin_mix_norm_dir_weigths1 = round(results_fin_mix_norm_dir_weigths, 2)
View(results_fin_mix_norm_dir_weigths1)
#write.csv(results_fin_mix_norm_dir_weigths1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2gamma_comp_varall_res.csv",row.names=TRUE )  

##### to check the convergence ###########
library(coda)
mcmc_obj <- as.mcmc(run.model_fin_mix_norm_dir_weigths_tau_comp)

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^mu\\[", varnames(mcmc_obj[[1]]))][,1:2]) # 
traceplot(mcmc_obj[, grep("^tau1\\[", varnames(mcmc_obj[[1]]))][,1:2]) #
traceplot(mcmc_obj[, grep("^p\\[", varnames(mcmc_obj[[1]]))][,1:2]) #

par(mfrow = c(4, 3) , mar = c(4, 4, 3, 1))  # Adjust margins if needed

traceplot(mcmc_obj[, grep("^delta12\\[", varnames(mcmc_obj[[1]]))][,1:12]) # first 6 clusters


par(mfrow = c(4, 3))  # Adjust margins if needed
traceplot(mcmc_obj[, c("delta_new", 
                       "total_mean", 
                       "total_var")])

par(mfrow = c(1, 1))  # Reset layout


## diagnostic plot ##
library(ggExtra)
library(MASS)
library(tidyverse)
library(jarbes)
library(R2jags)
library(ggstance)
require(gridExtra)
library(meta)
library(extraDistr)

plot.delta.pi <- function(pi.bias, mu, tau, 
                          level = c(0.5, 0.75, 0.95),
                          limits.x = c(0, 1),
                          limits.y = c(0, 10),
                          label.x = "Mixture probability",
                          label.y = "Difference between component means",
                          kde2d.n = 25,
                          color.line = "black",
                          title = paste("Posterior Contours (50%, 75% and 95%)"),
                          #marginals = TRUE,
                          bin.hist = 30,
                          color.hist = "white",
                          S = 2000,
                          color.data.points = "black",
                          alpha.data.points = 0.1,
                          ...)

{
  
  ## --- Posterior quantities (as in the paper)
  p.bias.1 = pi.bias 
  delta <- mu[,2] - mu[,1]

  ## Component-specific SDs
  min_sd <- min(tau[,1], tau[ ,2])
  
  cut.point <- 2 * mean(min_sd)
  
  dat.post <- data.frame(
    x = p.bias.1,
    y = delta
  )
  
  dat.post = dat.post[sample(1:S), ]
  
  ## Base plot (unchanged style)
  baseplot <- ggplot(dat.post, aes(x = x, y = y)) + 
    geom_point(size = 0.01, alpha = alpha.data.points,
               aes(color = "MCMC Samples"), color = color.data.points) +
    scale_x_continuous(name = label.x, limits = limits.x) +
    scale_y_continuous(name = label.y, limits = limits.y) +
    geom_hline(yintercept =  cut.point, linetype = "dashed") +
    theme_bw()
  
  ## KDE contours (same as original)
  
  dens <- kde2d(dat.post[ , 1], dat.post[ ,2], n = 50)
  dx <- diff(dens$x[1:2])
  dy <- diff(dens$y[1:2])
  sz <- sort(dens$z)
  c1 <- cumsum(sz) * dx * dy
  Levels.nonpar <- approx(c1, sz, xout = 1 - level)$y
  
  densdf <- data.frame(
    expand.grid(pi.bias = dens$x, bias = dens$y),
    z = as.vector(dens$z)
  )
  
  finalplot <- baseplot +
    geom_contour(data = densdf,
                 aes(pi.bias, bias, z = z),
                 breaks = Levels.nonpar)
  
  ggExtra::ggMarginal(finalplot,
                      type = "histogram",
                      fill = color.hist,
                      bins = bin.hist)
}


p.bias.1 <- run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$sims.list$p[,2]
mu     <- run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$sims.list$mu
tau   <- run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$sims.list$tau_comp


### some warnings are related to limits.x and limits.y
plot.delta.pi(
  pi.bias = p.bias.1,
  mu      = mu,
  tau   = tau,
  limits.x = c(0, 1),
  limits.y = c(0, 5.5)
)


### posterior probabilities that the treatment effect is below zero and that the treatment effect of a new stuy is below zero
Pr_mu_fin_mix_2norm_dir_weigths = mean(run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$sims.matrix[  ,"total_mean"] < 0 )
Pr_delta_new_fin_mix_2norm_dir_weigths = mean(run.model_fin_mix_norm_dir_weigths_tau_comp$BUGSoutput$sims.matrix[  ,"delta_new"] < 0 )

######## cluster assignments ######
fdd <- grepl("z", row.names(results_fin_mix_norm_dir_weigths))
mdd <- results_fin_mix_norm_dir_weigths[fdd,]

##### cluster assignment #########
clust <- round(mdd$`50%`, 0)

cluster_mix_norm_1 = which(clust == 1)
cluster_mix_norm_2 = which(clust == 2)

# 
# write.csv(cluster_mix_norm_1, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2cluster1_assign_fin_mix_norm_dw.csv",row.names=FALSE )  
# write.csv(cluster_mix_norm_2, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2cluster2_assign_fin_mix_norm_dw.csv",row.names=FALSE )  

### study-specific effect estimates 
fdd <- grepl("delta12", row.names(results_fin_mix_norm_dir_weigths))
mdd <- results_fin_mix_norm_dir_weigths[fdd,]
rel_eff_fin_mix_norm_dw <- mdd$`50%`
LB_rel_eff_fin_mix_norm_dw <- mdd$`2.5%`
UB_rel_eff_fin_mix_norm_dw <- mdd$`97.5%`
sd_rel_eff_fin_mix_norm_dw <- mdd$sd
Rhat_rel_eff_fin_mix_norm_dw <- mdd$Rhat

## cluster means 
mean_cl1 = mean(rel_eff_fin_mix_norm_dw[cluster_mix_norm_1])
mean_cl2 = mean(rel_eff_fin_mix_norm_dw[cluster_mix_norm_2])
mean_cl = cbind.data.frame(mean_cl1 , mean_cl2)

#write.csv(mean_cl, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2mean_cl.csv",row.names=FALSE )  

########## CLUSTER ASSIGNEMENT ##########
#### BASED ON  THE max(SC[i,j]) ####
results_fin_mix_norm_dir_weigths$ind <- row.names(results_fin_mix_norm_dir_weigths)

############# extraction of the parameters of interest ###########
f55sc_norm_mix <- grepl("SC", row.names(results_fin_mix_norm_dir_weigths))
m55sc_norm_mix <- results_fin_mix_norm_dir_weigths[f55sc_norm_mix,]
SC_norm_mix <- m55sc_norm_mix$mean

##### Distribution of data points to clusters ###

extra_col = rownames(m55sc_norm_mix)
m55sc_norm_mix$names = extra_col

prob_list <- data.frame(Column10 = m55sc_norm_mix[, 10], Column1 = m55sc_norm_mix[, 1])

split_dataframe <- function(df, chunk_size) {
  split(df, rep(1:ceiling(nrow(df) / chunk_size), each = chunk_size, length.out = nrow(df)))
}

split_prob_list <- split_dataframe(prob_list, nrow(pre_term_data_58))

split_prob_list

# Function to find the maximum probabilities for data points to be distributed into clusters ####
find_max_values_and_indices <- function(split_list) {
  num_rows <- nrow(split_list[[1]])
  num_sublists <- length(split_list)
  
  max_values <- numeric(num_rows)
  max_indices <- integer(num_rows)
  
  for (i in 1:num_rows) {
    values <- sapply(split_list, function(x) x[i, 2])
    max_value <- max(values, na.rm = TRUE)
    max_index <- which(values == max_value)[1]  
    max_values[i] <- max_value
    max_indices[i] <- max_index
  }
  
  return(list(max_values = max_values, max_indices = max_indices))
}

result <- find_max_values_and_indices(split_prob_list)
max_values <- result$max_values #### probabilities of cluster assignment
max_indices <- result$max_indices #### cluster assignments 
# 
# write.csv(max_values, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2cluster_prob_pred.csv",row.names=FALSE )  
# write.csv(max_indices, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2cluster_assignment.csv",row.names=FALSE )  

hist(rel_eff_fin_mix_norm_dw)
plot(density(rel_eff_fin_mix_norm_dw))


# Calculate densities for both clusters 
density2 <- density(rel_eff_fin_mix_norm_dw[cluster_mix_norm_2])
density1 <- density(rel_eff_fin_mix_norm_dw[cluster_mix_norm_1])

# Plot the first density
plot(density1, col = "red", lwd = 2, main = "Density Plot of Two Clusters",
     xlab = "Value", ylab = "Density", xlim= c(-3,0), ylim = c(0,5))

lines(density2, col = "blue", lwd = 2)


## Mean of each component ######
f33_fin_mix_norm_dir_weigths <- grepl("mu", row.names(results_fin_mix_norm_dir_weigths))
m33_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[f33_fin_mix_norm_dir_weigths,]
cl_mu_fin_mix_norm_dw <- m33_fin_mix_norm_dir_weigths$`50%`
LB_cl_mu_fin_mix_norm_dw <- m33_fin_mix_norm_dir_weigths$`2.5%`
UB_cl_mu_fin_mix_norm_dw <- m33_fin_mix_norm_dir_weigths$`97.5%`
Rhat_cl_mu_fin_mix_norm_dw <- m33_fin_mix_norm_dir_weigths$Rhat
prec_cl_mu_fin_mix_norm_dw <- UB_cl_mu_fin_mix_norm_dw - LB_cl_mu_fin_mix_norm_dw


#### Variance of each component
f44_fin_mix_norm_dir_weigths <- grepl("tau_k_sqr", row.names(results_fin_mix_norm_dir_weigths))
m44_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[f44_fin_mix_norm_dir_weigths,]
m44_fin_mix_norm_dir_weigths <- m44_fin_mix_norm_dir_weigths[!grepl("tau1", row.names(m44_fin_mix_norm_dir_weigths)),]
cl_tau2_fin_mix_norm_dw <- m44_fin_mix_norm_dir_weigths$`50%`
LB_cl_tau2_fin_mix_norm_dw <- m44_fin_mix_norm_dir_weigths$`2.5%`
UB_cl_tau2_fin_mix_norm_dw <- m44_fin_mix_norm_dir_weigths$`97.5%`
Rhat_cl_tau_fin_mix_norm_dw <- m44_fin_mix_norm_dir_weigths$Rhat
prec_cl_tau2_fin_mix_norm_dw <- UB_cl_tau2_fin_mix_norm_dw - LB_cl_tau2_fin_mix_norm_dw


########### SAVE THE REST OF THE RESULTS ###########

f11_fin_mix_norm_dir_weigths <- grepl("total_mean", row.names(results_fin_mix_norm_dir_weigths))
m11_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[f11_fin_mix_norm_dir_weigths,]
mu_fin_mix_norm_dw<- m11_fin_mix_norm_dir_weigths$`50%`
LB_mu_fin_mix_norm_dw <- m11_fin_mix_norm_dir_weigths$`2.5%`
UB_mu_fin_mix_norm_dw <- m11_fin_mix_norm_dir_weigths$`97.5%`
Rhat_mu_fin_mix_norm_dw <- m11_fin_mix_norm_dir_weigths$Rhat
prec_mu_fin_mix_norm_dw <- UB_mu_fin_mix_norm_dw - LB_mu_fin_mix_norm_dw

########### SAVE THE REST OF THE RESULTS ###########

f17_fin_mix_norm_dir_weigths <- grepl("delta_new", row.names(results_fin_mix_norm_dir_weigths))
m17_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[f17_fin_mix_norm_dir_weigths,]
delta_new_fin_mix_norm_dw<- m17_fin_mix_norm_dir_weigths$`50%`
LB_delta_new_fin_mix_norm_dw <- m17_fin_mix_norm_dir_weigths$`2.5%`
UB_delta_new_fin_mix_norm_dw <- m17_fin_mix_norm_dir_weigths$`97.5%`
Rhat_delta_new_fin_mix_norm_dw <- m17_fin_mix_norm_dir_weigths$Rhat
prec_delta_new_fin_mix_norm_dw <- UB_delta_new_fin_mix_norm_dw - LB_delta_new_fin_mix_norm_dw

f22_fin_mix_norm_dir_weigths <- grepl("total_var", row.names(results_fin_mix_norm_dir_weigths))
m22_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[f22_fin_mix_norm_dir_weigths,]
tau2_fin_mix_norm_dw <- m22_fin_mix_norm_dir_weigths$`50%`
LB_tau2_fin_mix_norm_dw <- m22_fin_mix_norm_dir_weigths$`2.5%`
UB_tau2_fin_mix_norm_dw <- m22_fin_mix_norm_dir_weigths$`97.5%`
Rhat_tau_fin_mix_norm_dw <- m22_fin_mix_norm_dir_weigths$Rhat
prec_tau2_fin_mix_norm_dw <- UB_tau2_fin_mix_norm_dw -  LB_tau2_fin_mix_norm_dw


fp_fin_mix_norm_dir_weigths <- grepl("p", row.names(results_fin_mix_norm_dir_weigths))
mp_fin_mix_norm_dir_weigths <- results_fin_mix_norm_dir_weigths[fp_fin_mix_norm_dir_weigths,]
mp_fin_mix_norm_dir_weigths <- mp_fin_mix_norm_dir_weigths[!grepl("pop", row.names(mp_fin_mix_norm_dir_weigths)),]
mp_fin_mix_norm_dir_weigths <- mp_fin_mix_norm_dir_weigths[!grepl("pred", row.names(mp_fin_mix_norm_dir_weigths)),]
pi_fin_mix_norm_dir_weigths <- mp_fin_mix_norm_dir_weigths$mean


lista_fin_mix_norm_dir_weigths <- cbind.data.frame(rel_eff_fin_mix_norm_dw,sd_rel_eff_fin_mix_norm_dw, LB_rel_eff_fin_mix_norm_dw,
                                                   UB_rel_eff_fin_mix_norm_dw,Rhat_rel_eff_fin_mix_norm_dw)

############ Overall results ##########
res_fin_mix_norm_dw <- cbind.data.frame(mu_fin_mix_norm_dw, LB_mu_fin_mix_norm_dw,UB_mu_fin_mix_norm_dw, 
                                        delta_new_fin_mix_norm_dw, LB_delta_new_fin_mix_norm_dw,UB_delta_new_fin_mix_norm_dw, 
                                        tau2_fin_mix_norm_dw, LB_tau2_fin_mix_norm_dw, UB_tau2_fin_mix_norm_dw,
                                        prec_mu_fin_mix_norm_dw, prec_delta_new_fin_mix_norm_dw,
                                        prec_tau2_fin_mix_norm_dw, 
                                        Rhat_mu_fin_mix_norm_dw, Rhat_delta_new_fin_mix_norm_dw,
                                        Rhat_tau_fin_mix_norm_dw, 
                                        Pr_mu_fin_mix_2norm_dir_weigths, Pr_delta_new_fin_mix_2norm_dir_weigths )

##### Clusters results ##########
clusters_fin_mix_norm_dw = cbind.data.frame(cl_mu_fin_mix_norm_dw, LB_cl_mu_fin_mix_norm_dw, UB_cl_mu_fin_mix_norm_dw,
                                            cl_tau2_fin_mix_norm_dw, LB_cl_tau2_fin_mix_norm_dw, UB_cl_tau2_fin_mix_norm_dw,
                                            prec_cl_mu_fin_mix_norm_dw, prec_cl_tau2_fin_mix_norm_dw,
                                            pi_fin_mix_norm_dir_weigths) 

extra_col_cl = row.names(clusters_fin_mix_norm_dw)
clusters_fin_mix_norm_dw$extra_col = extra_col_cl


# write.csv(res_fin_mix_norm_dw, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2gamma_comp_varres.csv",row.names=FALSE )  
# write.csv(clusters_fin_mix_norm_dw, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2gamma_comp_varclusters.csv",row.names=FALSE )  
# write.csv(lista_fin_mix_norm_dir_weigths, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\two_normal_mix_58\\dir_weights\\gamma_comp_var\\tauinv2_gamma\\2gamma_comp_varrel_eff.csv",row.names=FALSE )  


