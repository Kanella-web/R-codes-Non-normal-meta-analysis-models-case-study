##load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

######### Scaled-NORMAL MODEL FOR CONTIHMROUS OUTCOMES ##############
library(R2jags)
set.seed(1508)
writeLines("
model{

  ############################
  # STUDY LEVEL
  ############################
  for(i in 1:ns){
    mu1[i] <- theta.1[i]
    
    ## Likelihood
    y1[i] ~ dnorm(mu1[i], prec1[i])
    y2[i] ~ dnorm(mu2[i], prec2[i])
    
    mu2[i] <- mu1[i] + theta.2[i]*pooled_sd[i]

    prec1[i] <- n1[i] / pow(sd1[i], 2)
    prec2[i] <- n2[i] / pow(sd2[i], 2)

   
    pooled_sd[i] <- sqrt(((n1[i] - 1)*pow(sd1[i], 2) + (n2[i] - 1)*pow(sd2[i], 2)) / (n1[i] + n2[i] - 2))

    ## Scaled normal random effects (conditional formulation)
    theta.1[i] ~ dnorm(mu.1, pre.theta.1[i])
    theta.2[i] ~ dnorm(mu.2.1[i], pre.theta.2.1[i])

    #Conditional mean
    mu.2.1[i] <- mu.2 + rho * sigma.2 / sigma.1 * (theta.1[i] - mu.1)

    # Conditional precision
    pre.theta.2.1[i] <- pre.theta.2[i] / (1 - rho * rho)

    lambda[i] ~ dchisqr(df)
    w[i] <- df/lambda[i]

    pre.theta.1[i] <- pre.1 / w[i]
    pre.theta.2[i] <- pre.2 / w[i]

    # Probability of outlier
    p.w[i] <- step(w[i]-0.99)
  }

    ############################
    # HYPERPRIORS
    ############################

    ## Means (logit scale, logistic prior)
    mu.1 ~ dlogis(0, 2)
    mu.2 ~ dlogis(0, 2)

   # Dispersion parameters
    sigma.1 ~ dunif(0, 5)
    sigma.2 ~ dunif(0, 5)
    pre.1 <- 1 /(sigma.1*sigma.1)
    pre.2 <- 1 /(sigma.2*sigma.2)

  ## Correlation (Fisher transform)
  
   z ~ dnorm(0, pre.Fisher.rho)
   rho <- 2*exp(z)/(1+exp(z)) - 1

  pre.Fisher.rho <- 1/(sd.Fisher.rho * sd.Fisher.rho)
  ############################
  
   # Predictions ...
    mu.new[1] <- mu.1
    mu.new[2] <- mu.2

    Sigma.new[1, 1] <- pow(sigma.1, 2)
    Sigma.new[2, 2] <- pow(sigma.2, 2)
    Sigma.new[1, 2] <- rho * sigma.1 * sigma.2
    Sigma.new[2, 1] <- Sigma.new[1, 2]

    lambda.new ~ dchisqr(df)
    w.new <- df/lambda.new
    Omega.new[1:2,1:2] <- inverse(Sigma.new[1:2, 1:2]) / w.new

    # DEGREES OF FREEDOM
    ############################
    a <- 1/df.upper
    b <- 1/df.lower
    df <- 1/U
    U ~ dunif(a, b)

    theta.new[1:2] ~ dmnorm(mu.new[1:2], Omega.new[1:2 ,1:2])

    # Functional parameters
    alpha.0 <- (mu.2 - alpha.1 * mu.1)
    alpha.1 <- rho * sigma.2/sigma.1

}
", con = "HMR_scaled_normal_SMDs.txt")


dati <- list(
  ns  = nrow(pre_term_data_58),
  y1  = pre_term_data_58$mean_FT,
  y2  = pre_term_data_58$mean_EPT.VPT,
  sd1 = pre_term_data_58$sd_FT,
  sd2 = pre_term_data_58$sd_EPT.VPT,
  n1  = pre_term_data_58$n_FT,
  n2  = pre_term_data_58$n_EPT.VPT,
  df.lower = 3,            # Lower bound of the df's prior
  df.upper = 10 ,
  sd.Fisher.rho = 1.25
)


params <- c(
  # Means
  "mu.1", "mu.2",
  
  "alpha.0", "alpha.1",
  "rho",
  
  # Heterogeneity (latent-scale)
  "sigma.1",
  'sigma.2',
  
  # weights
  "Omega.new", "w.new",
  
  # Degrees of freedom
  "df",
  
  # Outliers
  "w", "p.w",
  
  # Predictions
  "theta.new",
  "theta.2"
)


fit <- jags(
  data = dati,
  inits = NULL,
  parameters.to.save = params,
  model.file = "HMR_scaled_normal_SMDs.txt",
  n.chains = 2,
  n.iter = 50000,
  n.burnin = 10000,
  DIC = TRUE
)

resultsHMR_SMDs = as.data.frame(fit$BUGSoutput$summary)
resultsHMR_SMDs = round(resultsHMR_SMDs , 4)
View(resultsHMR_SMDs)

#write.csv(resultsHMR_SMDs, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\HMR_58\\HMR_58_SMDs.csv",row.names=TRUE )  

### Pr of mu <0 and Pr of pred <- 0 #########
Pr_mu_HMR_SMDs = mean(fit$BUGSoutput$sims.matrix[  ,"mu.2"] < 0 )
Pr_pred_HMR_SMDs = mean(fit$BUGSoutput$sims.matrix[  ,"theta.new[2]"] < 0 )

library(dplyr)


post <- fit$BUGSoutput$sims.list
w <- post$w   # iterations × ns

# posterior summaries
w.med <- apply(w, 2, median)
w.u <- apply(w, 2, quantile, prob = 0.75)
w.l <- apply(w, 2, quantile, prob = 0.25)

studies <- pre_term_data_58$Study

dat.weights <- data.frame(
  x = studies,
  y = w.med,
  ylo = w.l,
  yhi = w.u
)

# outliers
all_studies_med_weight <- median(w.med)

## outliers are presented as studies which have 2times the median weight of the "homogeneous" studies
outliers <- w.med > 2 *all_studies_med_weight 
true_positions <- which(outliers)
true_positions

w.col <- dat.weights$x[1:58] %in% studies[true_positions]

w.col.plot <- ifelse(w.col, "red", "blue")

# Function to generate the forest plot in ggplot...
w.plot <- function(d){
  # d is a data frame with 4 columns
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi) )+ 
    #theme_bw()+
    geom_pointrange( colour=w.col.plot,
                     lwd =0.8)+
    coord_flip() + 
    #geom_hline(yintercept = 0, colour = "black", size = 0.5, lty =2)
    geom_hline(yintercept = 1, lty=2)+ 
    xlab("Study") +
    ylab("Scale weights") + theme_bw()
  return(p)
}                          

diagnostic.plot2 <- w.plot(dat.weights)
diagnostic.plot2



#####STUDY-SPECIFIC EFFECTS
fd <- grepl("theta.2", row.names(resultsHMR_SMDs))
md <- resultsHMR_SMDs[fd,]
rel_eff_HMR_SMDs <- md$`50%`
LB_rel_eff_HMR_SMDs <- md$`2.5%`
UB_rel_eff_HMR_SMDs <- md$`97.5%`
sd_rel_eff_HMR_SMDs <- md$sd
Rhat_HMR_SMDs <- md$Rhat

######## TO VISUALIZE THE DENSITY PLOT OF STUDY-SPECIFIC EFFECTS ##########
plot(density(rel_eff_HMR_SMDs))

fdn <- grepl("theta.new", row.names(resultsHMR_SMDs))
mdn <- resultsHMR_SMDs[fdn,]
pred_HMR_SMDs <- mdn$`50%`
pred_HMR_SMDs_LB <- mdn$`2.5%`
pred_HMR_SMDs_UB <- mdn$`97.5%`
pred_HMR_SMDs_sd <- mdn$sd
pred_HMR_SMDs_Rhat <- mdn$Rhat

### Pr of mu <0 and Pr of pred <- 0 #########
Pr_mu_HMR_SMDs = mean(fit$BUGSoutput$sims.matrix[  ,"theta.2[2]"] < 0 )
Pr_pred_HMR_SMDs = mean(fit$BUGSoutput$sims.matrix[  ,"theta.new[2]"] < 0 )


f1HMR <- grepl("mu.2", row.names(resultsHMR_SMDs))
m1HMR <- resultsHMR_SMDs[f1HMR,]
mu_HMR_SMDs <- m1HMR$`50%` 
LB_mu_HMR_SMDs <- m1HMR$`2.5%`
UB_mu_HMR_SMDs <- m1HMR$`97.5%`
Rhat_mu_HMR_SMDs <- m1HMR$Rhat
prec_mu_HMR_SMDs <-  UB_mu_HMR_SMDs - LB_mu_HMR_SMDs


f2HMR <- grepl("sigma.2", row.names(resultsHMR_SMDs))
m2HMR <- resultsHMR_SMDs[f2HMR,]
tau_HMR_SMDs <- m2HMR$`50%`
tau2_HMR_SMDs <- (m2HMR$`50%`)^2
LB_tau2_HMR_SMDs <- (m2HMR$`2.5%`)^2
UB_tau2_HMR_SMDs <- (m2HMR$`97.5%`)^2
Rhat_tau_HMR_SMDs <- m2HMR$Rhat
tau2_HMR_SMDs = round(tau2_HMR_SMDs,2)
prec_tau2_HMR_SMDs <-  UB_tau2_HMR_SMDs - LB_tau2_HMR_SMDs


lista_HMR_SMDs <- cbind.data.frame(rel_eff_HMR_SMDs,sd_rel_eff_HMR_SMDs,LB_rel_eff_HMR_SMDs, 
                                  UB_rel_eff_HMR_SMDs, Rhat_HMR_SMDs)


res_HMR_SMDs <- cbind.data.frame(mu_HMR_SMDs, LB_mu_HMR_SMDs, UB_mu_HMR_SMDs, tau2_HMR_SMDs,
                                LB_tau2_HMR_SMDs, UB_tau2_HMR_SMDs, pred_HMR_SMDs[2], pred_HMR_SMDs_LB[2],
                                pred_HMR_SMDs_UB[2], Rhat_mu_HMR_SMDs, Rhat_tau_HMR_SMDs, 
                                pred_HMR_SMDs_Rhat[2],
                                prec_mu_HMR_SMDs,prec_tau2_HMR_SMDs , pred_HMR_SMDs[2], 
                                Pr_mu_HMR_SMDs, Pr_pred_HMR_SMDs)


#write.csv(res_HMR_SMDs, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\HMR_58\\res_HMR.csv",row.names=FALSE )  

