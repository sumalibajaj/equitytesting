
data{
  int N;
  int n[N]; // number of LFD positives
  int n_tp[N]; // number of PCR positives 
  
  // prevalence priors
  vector<lower=0, upper=1>[N] prevalence_mean;
  vector<lower=0, upper=1>[N] prevalence_sd;
}

parameters {
  vector<lower=0, upper=1>[N] prevalence;
  real<lower=0, upper=1> a;
  real<lower=0, upper=1> b;
}

transformed parameters{
  vector<lower=0, upper=1>[N] theta;

  for(i in 1:N){
    theta[i] = a*prevalence[i]/(b + (a-b)*prevalence[i]);
  }
}


model{
  n_tp ~ binomial(n, theta);
  
  prevalence ~ normal(prevalence_mean, prevalence_sd);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = binomial_lpmf(n_tp[i]| n[i], theta);
  }
}
