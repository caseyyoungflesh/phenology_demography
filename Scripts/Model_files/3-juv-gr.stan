// juv ~ gr (sens)

data {
int<lower=0> N;                      // number of data points
int<lower=0> Nst;                    // number of species/stations
int<lower=0> Nsp;
vector[N] juv_obs; 
vector<lower=0>[N] juv_sd;
int<lower=0> st_id[N];
int<lower=0> sp_id[N];
int<lower=0> sp_id2[Nst];
vector[N] gr;
}

parameters {
vector[Nst] alpha_raw;
vector[Nst] beta_raw;
vector[Nsp] mu_alpha;
real<lower=0> sigma_alpha;
vector[Nsp] mu_beta;
real<lower=0> sigma_beta;
real mma;
real<lower=0> sma;
real mmb;
real<lower=0> smb;
real<lower=0> sigma;
vector[N] juv_raw;
}

transformed parameters {
vector[N] mu;
vector[Nst] alpha;
vector[Nst] beta;
vector[N] juv;

alpha = alpha_raw * sigma_alpha + mu_alpha[sp_id2];
beta = beta_raw * sigma_beta + mu_beta[sp_id2];

// intercept and slope for each species/site
mu = alpha[st_id] + beta[st_id] .* gr;

// implies juv ~ normal(mu, sigma);
juv = juv_raw * sigma + mu;
}

model {
// priors
sigma ~ normal(0, 20);
sigma_alpha ~ normal(0, 10);
sigma_beta ~ std_normal();
mma ~ normal(200, 20);
sma ~ normal(0, 30);
mmb ~ std_normal();
smb ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();
juv_raw ~ std_normal();

//centered
mu_alpha ~ normal(mma, sma);
mu_beta ~ normal(mmb, smb);

// obs model
juv_obs ~ normal(juv, juv_sd);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(juv, juv_sd);
}
