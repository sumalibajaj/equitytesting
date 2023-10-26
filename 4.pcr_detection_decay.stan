
data {
  int n_obs;
  real days[n_obs]; // number of days since LFD+
  int n_chunks; // number of different parameter chunks
  int chunks[n_obs]; // time chunk index associated with each obs (i.e. 1, 2, 3, ..., n_chunks)
  int n_pcr_positive[n_obs];
  int n_lfd_positive[n_obs];
}

parameters {
  real<lower=0, upper=1> theta[n_chunks];
  real<lower=0> lambda[n_chunks];
}

transformed parameters {
  real<lower=0, upper=1> prob[n_obs];
  
  for(i in 1:n_obs){
    int chunk = chunks[i];
    prob[i] = theta[chunk] * exp(-lambda[chunk] * days[i]);
  }
}

model {
  for(i in 1:n_obs) {
    n_pcr_positive[i] ~ binomial(n_lfd_positive[i], prob[i]);
  }
  
  lambda ~ cauchy(0, 10);
}
