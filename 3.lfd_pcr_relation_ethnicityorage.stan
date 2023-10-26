
data{
  int N;
  int eth_N; // number of ethnic/age groups
  int n[N]; // number of LFD positives
  int n_tp[N]; // number of PCR positives 
  
  int<upper = eth_N> eth[N]; // indicator vector for ethnicity/age group
  
  // prevalence priors
  vector<lower=0, upper=1>[N] prevalence_mean;
  vector<lower=0, upper=1>[N] prevalence_sd;
}

parameters {
  vector<lower=0, upper=1>[N] prevalence;
  real<lower=0, upper=1> a[eth_N];
  real<lower=0, upper=1> b[eth_N];
}

transformed parameters{
  vector<lower=0, upper=1>[N] theta;

  for(i in 1:N){
    theta[i] = a[eth[i]]*prevalence[i]/(b[eth[i]] + (a[eth[i]]-b[eth[i]])*prevalence[i]);
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
