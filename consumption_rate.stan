data {
  int<lower=0> N;        // observations
  int<lower=0> K;        // groups (foraging behaviour)
  int<lower=0> Sb;       // number of bird species
  int<lower=0> Sp;       // number of plant species
  int<lower=0> Ns;       // number of source articles

  int group[N];          // 1 = gulper, 2 = masher
  matrix[Sb, K] group_mat; // matrix to identify group of each bird species
  int plant[N];          // plant species ids
  int bird[N];           // bird species ids
  int source[N];
  vector[N] y;
  vector[N] sizediff;

  // define parameters for prior distributions (all normal)
  real prior_a_sd;
  real prior_b_sd;
  real prior_phi_sd;
  real lower_u;
  real upper_u; // u has uniform prior over (lower_u, upper_u)
  real prior_sigma_sd;
}

parameters {
  // all parameters vary with foraging behaviour
  vector[K] a_raw;  // intercept at logit scale
  vector[K] b_raw;  // slope at logit scale
  vector<lower = lower_u, upper = upper_u>[K] u; // upper limit for logistic curve
  vector<lower=0>[K] phi_raw; // dispersion parameter for gamma

  // random effect for plant and bird species
  real<lower = 0> sigma_plant_raw;
  vector<lower = 0>[K] sigma_bird_raw;
  real<lower = 0> sigma_source_raw;
  vector[Sp] e_plant_raw;
  vector[Sb] e_bird_raw;
  vector[Ns] e_source_raw;
}

transformed parameters {
  vector[K] a = a_raw * prior_a_sd;  // intercept at logit scale
  vector[K] b = b_raw * prior_b_sd;  // slope at logit scale
  vector<lower=0>[K] phi = phi_raw * prior_phi_sd; // dispersion parameter for gamma

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
  phi_raw ~ std_normal();
  // (not specifying a prior for u assigns a flat prior over its domain (0, upper_u))

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
    real mu = u[group[n]] * inv_logit(lp);

    // get dispersion
    real phi_n = phi[group[n]];

    // likelihood
    y[n] ~ gamma(1 / phi_n, 1 / (mu * phi_n));
  }
}