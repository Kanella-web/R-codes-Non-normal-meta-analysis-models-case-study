
data {
  
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

  array[ns] real delta12;
}

transformed parameters {
  real tau_sqr;
  real mu;
  real skew;
  real kurt;
  real beta;
  real delt;
  array[ns] real pooled_sd;
  
  beta = sqrt(2/pi()) ; 
  
  delt = ( alpha / (sqrt(1 + (alpha^2))) ) ; 
  
  mu = (xi + (beta * omega * delt) ) ; 
  
  tau_sqr = ( (omega^2) * ( 1 - ((beta^2) * (delt^2)) ) ) ;
  
  kurt = ((2*(pi() - 3)) * (((delt*beta)^4) / ((1 - (beta*delt)^2)^2) )) ;
  
  skew = ( (4 - pi()) / 2 ) * ( (( beta*delt)^3) / (( 1 - ((beta*delt)^2) )^(1.5)) ) ; 
  
  for (i in 1:ns) {
    pooled_sd[i] = sqrt(((n1[i] - 1) * pow(sd1[i], 2) + (n2[i] - 1) * pow(sd2[i], 2)) / (n1[i] + n2[i] - 2));
  }
}

model {
  for (i in 1:ns) {
    
    m[i] ~ normal(0, 100);
    
    y1[i] ~ normal(m[i], sd1[i]/ sqrt(n1[i]));
    
    y2[i] ~ normal(m[i] + delta12[i]*pooled_sd[i], sd2[i]/ sqrt(n2[i]));

    delta12[i] ~ skew_normal(xi, omega, alpha);
  }
  
  xi ~ normal(0, 100);
  omega ~ uniform(0, 10);
  alpha ~ normal(0, 5);
}

generated quantities{
  real pred;
  pred = skew_normal_rng(xi, omega, alpha);
  
}

