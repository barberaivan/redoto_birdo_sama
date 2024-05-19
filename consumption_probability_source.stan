data {
  int<lower=0> N;        // observations
  int<lower=0> K;        // groups (foraging behaviour)
  int<lower=0> Sb;       // number of bird species
  int<lower=0> Sp;       // number of plant species
  int<lower=0> Ns;       // number of source articles

  int group[N];          // 1 = gulper, 2 = masher
  int plant[N];          // plant species ids
  int bird[N];           // bird species ids
  int source[N];
  int y[N];
  vector[N] sizediff;
  matrix[Sb, K] group_mat;
  // binary, to perform matrix multiplication with sigma_bird

  // define parameters for prior distributions (all normal)
  real prior_a_sd;
  real prior_b_sd;
  real lower_u;
  real prior_sigma_sd;
}

parameters {
  // all parameters vary with foraging behaviour
  vector[K] a_raw;  // intercept at logit scale
  vector[K] b_raw;  // slope at logit scale
  vector<lower = lower_u, upper = 1>[K] u; // upper limit for logistic curve

  // random effect for plant and bird species
  real<lower = 0> sigma_plant_raw;
  vector<lower = 0>[K] sigma_bird_raw;  // bird species sd varies among groups
  real<lower = 0> sigma_source_raw;
  vector[Sp] e_plant_raw;
  vector[Sb] e_bird_raw;
  vector[Ns] e_source_raw;
}

transformed parameters {
  vector[K] a = a_raw * prior_a_sd;  // intercept at logit scale
  vector[K] b = b_raw * prior_b_sd;  // slope at logit scale

  real<lower = 0> sigma_plant = sigma_plant_raw * prior_sigma_sd;
  vector<lower = 0>[K] sigma_bird = sigma_bird_raw * prior_sigma_sd;
  real<lower = 0> sigma_source = sigma_source_raw * prior_sigma_sd;
  vector[Sp] e_plant = e_plant_raw * sigma_plant;
  vector[Sb] sigma_bird_rep = group_mat * sigma_bird;
  vector[Sb] e_bird = e_bird_raw .* sigma_bird_rep;
  vector[Ns] e_source = e_source_raw * sigma_source;
}

model {
  // priors
  a_raw ~ std_normal();
  b_raw ~ std_normal();
  // (not specifying a prior for u assigns a flat prior over its domain (lower_u, 1))

  e_plant_raw ~ std_normal();
  e_bird_raw ~ std_normal();
  e_source_raw ~ std_normal();
  sigma_plant_raw ~ std_normal();
  sigma_bird_raw ~ std_normal();
  sigma_source_raw ~ std_normal();

  // likelihood (loop over observations)
  for(n in 1:N) {
    // compute mu
    real lp = a[group[n]] +
              b[group[n]] * sizediff[n] +
              e_plant[plant[n]] +
              e_bird[bird[n]] +
              e_source[source[n]];
    real p = u[group[n]] * inv_logit(lp);

    // likelihood
    y[n] ~ bernoulli(p);
  }
}