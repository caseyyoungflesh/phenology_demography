// juv ~ gr (sens)

data {
int<lower=0> N;                      // number of data points
int<lower=0> Nst;                    // number of species/stations
int<lower=0> Nsp;
vector[N] juv; 
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
vector<lower=0>[Nsp] sigma;
real<lower=0> nu;
real<lower=0> lambda;
real<lower=0> kappa;
}

transformed parameters {
vector[N] mu;
vector[Nst] alpha;
vector[Nst] beta;

alpha = alpha_raw * sigma_alpha + mu_alpha[sp_id2];
beta = beta_raw * sigma_beta + mu_beta[sp_id2];

// intercept and slope for each species/site
mu = alpha[st_id] + beta[st_id] .* gr;
}

model {
// priors
sigma_alpha ~ normal(0, 10);
sigma_beta ~ std_normal();
mma ~ normal(200, 20);
sma ~ normal(0, 30);
mmb ~ std_normal();
smb ~ std_normal();
alpha_raw ~ std_normal();
beta_raw ~ std_normal();
nu ~ gamma(2, 0.1);             // degrees of freedom parameter
lambda ~ normal(0, 10);
kappa ~ normal(0, 10);

//centered
mu_alpha ~ normal(mma, sma);
mu_beta ~ normal(mmb, smb);

sigma ~ normal(lambda, kappa);

juv ~ student_t(nu, mu, sigma[sp_id]);
}

generated quantities {
real y_rep[N];

y_rep = student_t_rng(nu, mu, sigma[sp_id]);
}
