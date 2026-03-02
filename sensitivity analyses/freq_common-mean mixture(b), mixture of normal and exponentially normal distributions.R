
######## Frequentist models #########

library(metafor)
library(maxLik)
inference=function (effect, sigmasq, start_v=-1, start_left=1, start_right=1) 
{
  REML=rma.uni(yi=effect, vi=sigmasq, method="REML") ## DL
  start3=c(REML$b, log((REML$tau2)^0.5+0.05), start_v) ## DL
  start4=c(REML$b, log((REML$tau2)^0.5+0.05), start_left, start_right) ## DL
  stage3_1=optim(start3, LL_3parameter, effect=effect, sigmasq=sigmasq, control=list(fnscale=-1), method="Nelder-Mead")$par
  stage3_1=c(mu=stage3_1[1], logtau=stage3_1[2], logv=stage3_1[3])
  mle3=maxLik(LL_3parameter, start=stage3_1, effect=effect, sigmasq=sigmasq, fixed=c(FALSE,FALSE,FALSE))
  stage4_1=optim(start4, LL_4parameter, effect=effect, sigmasq=sigmasq, control=list(fnscale=-1), method="Nelder-Mead")$par
  stage4_1=c(mu=stage4_1[1], logtau=stage4_1[2], logright=stage4_1[3],  logleft=stage4_1[4] )
  mle4=maxLik(LL_4parameter, start=stage4_1, effect=effect, sigmasq=sigmasq, fixed=c(FALSE,FALSE,FALSE,FALSE))
  print("3 parameter symmetric model inference")
  print(summary(mle3))
  print("4 parameter skew model inference")
  print(summary(mle4))
}

LL_3parameter=function (x, effect, sigmasq) 
{
  # mixture of 2 normal symmetric
  mu = x[1]
  tau = exp(x[2])
  tausq = tau^2
  v = exp(x[3])
  vsq=v^2
  const=.5*log(2.*pi)
  
  objf=effect
  # to set up a vector
  for (k in 1:length(effect))
  {
    sigsq = sigmasq[k] + tausq
    uv=sigsq+vsq
    sig = sqrt(sigsq)
    p = sigsq/(sigsq+vsq)
    y = effect[k] - mu 
    term=.5*y^2/uv
    xlog=(1.-p)*exp(-.5*y^2*vsq/((sigsq+vsq)*sigsq))/sig+p/sqrt(uv)
    objf[k]= -term+log(xlog)-const
  }
  # for Newton-raphson (NR) return scalar objf, ie return(sum(objf))
  return(sum(objf))
}

LL_4parameter =function (x, effect, sigmasq) 
{
  # mixture of normal and lagged normal
  mu = x[1]
  tau = exp(x[2])
  #if(tau < 0)return(NA)
  tausq = tau^2
  expright = exp(x[3])
  expleft = exp(x[4])
  #if(expright == Inf)return(NA)
  #if(expleft == Inf)return(NA)
  abfact=1./(expright+expleft)
  
  objf=effect
  # to set up a vector
  for (k in 1:length(effect))
  {
    sigsq = sigmasq[k] + tausq
    sig = sqrt(sigsq)
    p = sigsq/(sigsq+expright^2+expleft^2)
    y = effect[k] - mu 
    yl=y-expright+expleft
    factor = -0.5*log(sigsq) - 0.5*y^2/sigsq - 0.5*log(2.*pi)
    f3 = exp(factor)
    #          arg1 = yl/sig - sig/expright
    #          arg2 = -yl/sig - sig/expleft
    arg1num=yl*expright-sigsq
    arg1den=sig*expright
    arg2num=-yl*expleft-sigsq
    arg2den=sig*expleft
    arg1=arg1num/arg1den
    arg2=arg2num/arg2den
    
    if(!is.finite(arg1)){
      if(arg1num > 0){
        arg1=Inf
      }else{
        arg1=-Inf
      }
    }
    if(!is.finite(arg2)){
      if(arg2num > 0){
        arg2=Inf
      }else{
        arg2=-Inf
      }
    }
    if (arg1 < -8.) {
      f1 = abfact*exp(-0.5*yl^2/sigsq)/(sqrt(2.*pi)*abs(arg1))
    }else{
      bigphi = pnorm(arg1)
      f1 = abfact*exp(0.5*sigsq/expright^2-yl/expright)*bigphi
    }
    if (arg2 < -8.) {
      f2 = abfact*exp(-0.5*yl^2/sigsq)/(sqrt(2.*pi)*abs(arg2))
    }else{
      bigphi = pnorm(arg2)
      f2 = abfact*exp(0.5*sigsq/expleft^2+yl/expleft)*bigphi
    }
    f = (1.-p)*f3 + p*(f1+f2)
    objf[k]= log(f)
  }
  # for Newton-raphson (NR) return scalar objf, ie return(sum(objf))
  return(sum(objf))
}

## load the data
pre_term_data_58 = read.csv("https://github.com/Kanella-web/R-codes-Non-normal-meta-analysis-models-case-study/raw/refs/heads/main/preterm%20data/pre_term_data58.csv")

data_SMDs = escalc(measure = "SMD", m1i = pre_term_data_58$mean_EPT.VPT, m2i =  pre_term_data_58$mean_FT,
                   sd1i = pre_term_data_58$sd_EPT.VPT, sd2i =  pre_term_data_58$sd_FT,
                   n1i = pre_term_data_58$n_EPT.VPT, n2i = pre_term_data_58$n_FT)

set.seed(1508)
inference(data_SMDs$yi, data_SMDs$vi)

#### mixture of two normal distributions(b) model
#### mu_3P
mu_3P = -0.87165
LB_mu_3P = mu_3P - 1.96*0.05159
UB_mu_3P = mu_3P + 1.96*0.05159
prec_mu_3P = UB_mu_3P - LB_mu_3P
#### u_3P
tau_3P = exp(-1.39621)
LB_tau_3P = tau_3P - 1.96*0.18569
UB_tau_3P = tau_3P + 1.96*0.18569

#### u2_3P = sigma2 + tau2
tau2_3P = (tau_3P)^2
LB_tau2_3P = (LB_tau_3P)^2
UB_tau2_3P = (UB_tau_3P)^2
prec_tau2_3P = UB_tau2_3P - LB_tau2_3P

#### additional variance of outlying studies
v = exp(-1.02217 )
LB_v = v - 1.96*0.80362
UB_v = v + 1.96*0.80362

v2 = (v)^2
LB_v2 = (LB_v)^2
UB_v2 = (UB_v)^2
prec_v2 = UB_v2 - LB_v2


res_3P = data.frame(mu_3P, LB_mu_3P, UB_mu_3P, tau_3P, LB_tau_3P, UB_tau_3P,
                    tau2_3P, LB_tau2_3P, UB_tau2_3P, v2,  LB_v2,UB_v2,
                    prec_mu_3P, prec_u2_3P, prec_v2)

# write.csv(res_3P, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\3_parameter_common_mean_mix_58\\res_3P_common_mean_mix.csv",row.names=FALSE )  

#### mixture of normal and exponentially modified normal distributions 
#### mu_4P 
mu_4P = -0.85362 
LB_mu_4P = mu_4P - 1.96*0.04871
UB_mu_4P = mu_4P + 1.96*0.04871
prec_mu_4P = UB_mu_4P - LB_mu_4P

#### tau_4P
tau_4P = exp(-1.22709)
LB_tau_4P = tau_4P - 1.96*0.14138
UB_tau_4P = tau_4P + 1.96*0.14138

#### tau2_4P
tau2_4P = (tau_4P)^2
LB_tau2_4P = (LB_tau_4P)^2
UB_tau2_4P = (UB_tau_4P)^2
prec_tau2_4P = UB_tau2_4P - LB_tau2_4P

######### for skewness and kurtosis parameters 
a_1 = exp(-0.06535 ) ### 1/a
a = 1/ a_1          ### a
b_1 = exp(0.50509)  #### 1/b
b = 1 / b_1         ### b
### alpha / r = a
### alpha*r = b (1)
### alpha = r*a
### (1)=> r*a*r = b
### r2 = b/a
r2 = b/a
r = sqrt(r2)
skew = r  ### skewness 
### alpha = r*a
alpha = a*r
### kurt = 1/alpha
kurt = 1 / alpha


res_4P = data.frame(mu_4P, LB_mu_4P, UB_mu_4P, tau_4P, LB_tau_4P, UB_tau_4P,
                    tau2_4P, LB_tau2_4P, UB_tau2_4P, skew, kurt,
                    prec_mu_4P, prec_tau2_4P)

#write.csv(res_4P, "C:\\Users\\Lela Panag\\Desktop\\2nd PhD article\\data\\pre_term_data58\\frequentist models_58\\4_parameter_skewed_model_58\\res_4P_skewed_model.csv",row.names=FALSE )  

