data {
  int N;         // for 2 x 2 contingency table, N is 4
  int K;         // number of coefficients in design matrix (4 for saturated model)
  int C;         // number cancers
  int y[C, N];   // contingency table is cancer x 4 
  matrix[N, K] x; // design matrix
}
transformed data {
  vector[C] total_count;
  for(c in 1:C) {
    // cancer-specific totals
    total_count[c] = sum(y[c]); 
  }
}
parameters {
  vector[K] beta[C];
  vector[K-1] mu_beta;
  real<lower=0> sigma;
  //vector[N] gamma;
}
transformed parameters {
  vector<lower=0>[N] mu[C];
  for(c in 1:C) {
    mu[c] = total_count[c]*exp(x * beta[c]);
  }
}
model {
  // No hierarchical model on intercept
  // why?  intercept is the grand mean
  beta[, 1] ~ normal(0, 5);
  // Hierarchical model on marginal and interaction effects
  for(c in 1:C) {
    beta[c, 2:K] ~ normal(mu_beta, sigma);
  }
  mu_beta ~ normal(0, 5);
  sigma ~ cauchy(0, 2.5);
  for(c in 1:C) {
    y[c] ~ poisson(mu[c]);
  }
}
generated quantities {
  vector[2] marginal[C];
  vector[N] pi[C];               // C x 4 matrix
  for(c in 1:C){
    // assigns length-4 vector to each element of pi
    pi[c] = mu[c] / sum(mu[c]);  
  }
  // for saturated model of 2 x 2 table, helpful to obtain marginal
  // prevalence of RAS and partner gene
  for(c in 1:C){
    marginal[c, 1] = pi[c, 2] + pi[c, 4];  // partner gene
    marginal[c, 2] = pi[c, 3] + pi[c, 4];  // RAS
  }
}

