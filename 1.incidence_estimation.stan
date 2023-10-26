
data {
  int N;
  int N_effective;
  int Tstart; // pass in time after which we want to estimate incidence
  vector<lower=0>[N] cases_reported;
  matrix[N_effective, N] Omega;
  int nreporting_period;
  int reporting_period[N];
  
  // priors on reporting fraction
  real<lower=0> a;
  real<lower=0> b;
  
  // prevalence priors
  vector<lower=0>[N - Tstart + 1] prevalence_mean;
  vector<lower=0>[N - Tstart + 1] prevalence_sd;
}

parameters {
  real<lower=0> sigma;
  vector<lower=0, upper=1>[nreporting_period] reporting_fraction;
  vector<lower=0>[N_effective] prevalence;
}

transformed parameters {
  vector<lower=0>[N] incidence;
  for(t in 1:N) {
    incidence[t] = cases_reported[t]/reporting_fraction[reporting_period[t]];
  }
}

model {
  for(t in 1:N_effective) {
    prevalence[t] ~ normal(Omega[t] * incidence, sigma);
  }
  sigma ~ normal(0, 1);
  reporting_fraction ~ beta(a, b);
  prevalence ~ normal(prevalence_mean, prevalence_sd);
}

generated quantities {
  vector<lower=0>[N_effective] prevalence_sim;
  for(t in 1:N_effective)
    prevalence_sim[t] = Omega[t] * incidence;
}
