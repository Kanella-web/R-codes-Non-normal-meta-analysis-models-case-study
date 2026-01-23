
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
  
  array[ns] real m;                    // baseline effect
  real xi;                             // Mean of the between-study effect
  real<lower=0> omega;                 // Between-study standard deviation
  array[ns] real delta12;   
  real<lower=2.5, upper=1000> nu;
  
}

transformed parameters {
  array[ns] real pooled_sd;
  real mu;
  real tau_sqr;
 // real kurt;

  mu = xi;
  tau_sqr = ((nu / (nu-2)) * (omega^2));
  
  //kurt = (6 / (nu-4)) ;, for nu>4
  // or
  // if (nu > 4) {
  //     kurt = (6 / (nu-4)) ;
  //   } else {
  //     kurt = 0;  // or set to `nan`
  //    }

  for (i in 1:ns) {
    pooled_sd[i] = sqrt(((n1[i] - 1) * pow(sd1[i], 2) + (n2[i] - 1) * pow(sd2[i], 2)) / (n1[i] + n2[i] - 2));
  }
}


model{
  for(i in 1:ns) {
    
    m[i] ~ normal(0, 100);
    
    y1[i] ~ normal(m[i], sd1[i]/ sqrt(n1[i]));
    
    y2[i] ~ normal(m[i] + delta12[i]*pooled_sd[i], sd2[i]/ sqrt(n2[i]));
    
    delta12[i] ~ student_t(nu, xi, omega);
    
  }
 
  xi ~ normal(0, 100);
  omega ~ uniform(0, 10);
  nu ~ exponential(0.10);
}


generated quantities{
  real pred;
  pred = student_t_rng(nu, xi, omega);
}
