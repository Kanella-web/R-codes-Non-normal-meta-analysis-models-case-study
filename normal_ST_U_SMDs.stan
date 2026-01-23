functions {
  real skew_t_lpdf(real y, real xi, real omega, real nu, real alpha) {
    return log(2) - log(omega) +
      student_t_lpdf((y - xi)/omega | nu, 0, 1) +
      student_t_lcdf((alpha*(y - xi) / omega)*sqrt(((nu + 1)*omega*omega) /
                                                     (nu*omega*omega + (y - xi)*(y - xi))) | (nu + 1), 0, 1);
  }
}

data{
  
  int<lower=1> ns;                     // Number of observations
  array[ns] real y1;                    // Observed outcomes for group 1
  array[ns] real y2;                    // Observed outcomes for group 2
  array[ns] int n1;                    // number of observations for group 1
  array[ns] int n2;                    // number of observations for group 2
  array[ns] real<lower=0> sd1;         // Standard deviations for group 1
  array[ns] real<lower=0> sd2;         // Standard deviations for group 2
  
}

parameters {
  array[ns] real m;
  real xi;
  real<lower=0> omega;
  real alpha;
  real<lower = 4.5, upper=1000> nu;
  array[ns] real delta12;
  real pred;
}


transformed parameters {
  
  array[ns] real pooled_sd;
  real mu;
  real tau_sqr;
  real skew ;
  real kurt;
  real beta;
  real delt ;
  real g1;
  real g2;
  real g3;
  real g4;
  real g5;
  //real nu; 
  
  //nu = (ns - 2) ;
  
  beta = ( sqrt(nu/pi())*((tgamma(0.5*(nu-1))) / (tgamma(0.5*(nu)))) );
  
  delt = ( alpha / sqrt(1+ (alpha^2)) ) ;
  
  g1 = (beta*delt) ;
  
  g2 = ( (nu / (nu-2)) - ((beta^2)*(delt^2)) ) ;
  
  g3 = beta * delt * ( (3*nu /((nu-3)*(nu-2)) )  - (delt^2) * ( (nu/(nu-3)) - 2*(beta^2) ) ) ;

  g4 = ( (3*(nu^2) / ((nu-2)*(nu-4))) - 6*(delt^2)*(beta^2)* ( (nu*(nu-1)) / ((nu-2)*(nu-3)) ) + (delt^4) *
           (beta^2) * ( ((4*nu) / (nu-3)) - 3*(beta^2) ) ) ;

  mu = xi + omega*(g1) ;
  
  tau_sqr = (omega^2)*(g2) ;
  
  skew = ( (g3) / ((g2)^(1.5)) ) ;

  g5 =  ((g4) / ((g2)^2))  ;

  kurt = (g5 - 3) ;

  // if (nu > 4) {
  //   g4 = (3*nu^2 / ((nu-2)*(nu-4))) - 6*(delt^2)*(beta^2)*(nu*(nu-1)/((nu-2)*(nu-3))) + delt^4 * beta^2 * ((4*nu/(nu-3)) - 3*beta^2);
  //   g5 = g4 / (g2^2);
  //   kurt = g5 - 3;
  // } else {
  //   kurt = 0;  // or set to `nan`
  // }
  // 
  // if (nu > 3) {
  //   g3 = beta * delt * ( (3*nu / ((nu-3)*(nu-2))) - delt^2 * ((nu/(nu-3)) - 2*beta^2) );
  //   skew = g3 / pow(g2, 1.5);
  // } else {
  //   skew = 0;  // or set to `nan` to indicate undefined
  // }
  // 
  for (i in 1:ns) {
    pooled_sd[i] = sqrt(((n1[i] - 1) * pow(sd1[i], 2) + (n2[i] - 1) * pow(sd2[i], 2)) / (n1[i] + n2[i] - 2));
  }
  
}

model{
  for(i in 1:ns) {
    
    m[i] ~ normal(0, 100);
    
    y1[i] ~ normal(m[i], sd1[i]/ sqrt(n1[i]));
    
    y2[i] ~ normal(m[i] + delta12[i]*pooled_sd[i], sd2[i]/ sqrt(n2[i]));

    delta12[i] ~ skew_t(xi, omega, nu, alpha);
  }
  
  xi ~ normal(0, 100);
  omega ~ uniform(0, 10);
  nu ~ exponential(0.10);
  alpha ~ normal(0, 5);
  pred ~ skew_t(xi, omega, nu, alpha);
}

