data {
    int<lower=1> Ntotal; // number of observations
    int<lower=1> Ncol; // total number of columns in model matrix
    matrix[Ntotal, Ncol] X; // model matrix
    int y[Ntotal]; // response variable binomial distributed
}
// transformed data {
// }
parameters {
  // parameters to estimate in the model
    vector[Ncol] betas; // regression parameters
}
transformed parameters {
  vector[Ntotal] mu; // fitted values from linear predictor
  mu = X * betas; 
  mu = inv_logit(mu);
}
model {
  betas[1] ~ normal(0, 1.5); //prior for the betas
  betas[2:Ncol] ~ normal(0, 0.5);
  // likelihood function
  y ~ bernoulli(mu);
}
generated quantities {
  int y_new[Ntotal];
  vector[Ntotal] log_lik;
  y_new = bernoulli_rng(mu);
  for (i in 1:Ntotal) log_lik[i] = bernoulli_lpmf(y[i] | mu[i]);
}
