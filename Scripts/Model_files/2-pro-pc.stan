// pro ~ PC

data {
int<lower=0> N;                        // number of data points
int<lower=0> Nsp;                      // number of species
int<lower=0> Nst;                      // number of species/stations
int<lower=0> y[N]; 
int<lower=0> ntot[N];
int<lower=0> sp_id[N];
int<lower=0> sp_id2[Nst];             // species id for each site
int<lower=0> st_id[N];
vector[N] PC1_1;
vector[N] PC1_2;
vector[N] PC2_1;
vector[N] PC2_2;
vector[N] EF;
}

parameters {
real<lower=0> nu;
real mma;
real<lower=0> sma;
vector[Nsp] mu_alpha;
real<lower=0> sigma_alpha;
vector[Nst] alpha_raw;
real mu_beta1;
real mu_beta2;
real mu_gamma1;
real mu_gamma2;
real mu_kappa;
real<lower=0> sigma_kappa;
vector<offset = mu_kappa, multiplier = sigma_kappa>[Nsp] kappa;
vector<lower=0>[4] sigma_bg;
cholesky_factor_corr[4] L_Rho_bg;             // cholesky factor of corr matrix
matrix[4, Nsp] z_bg;                          // z-scores
vector<lower=0>[Nsp] sigma_eps;
vector[N] epsilon;
real mu_se;
real<lower=0> sd_se;
}

transformed parameters {
vector[N] logit_p;
vector[Nst] alpha;
matrix[Nsp, 4] bg;
matrix[4, 4] Rho_bg;               // corr matrix
vector[Nsp] beta1;
vector[Nsp] beta2;
vector[Nsp] gamma1;
vector[Nsp] gamma2;

// implies alpha ~ normal(mu_alpha[sp_id2], sigma_alpha);
alpha = alpha_raw * sigma_alpha + mu_alpha[sp_id2];

// cholesky factor of covariance matrix multiplied by z score
// implies gt ~ MVN(0, sigma)
bg = (diag_pre_multiply(sigma_bg, L_Rho_bg) * z_bg)';
// implies Rho = L_Rho * L_Rho';
Rho_bg = multiply_lower_tri_self_transpose(L_Rho_bg);

beta1 = mu_beta1 + bg[,1];
beta2 = mu_beta2 + bg[,2];
gamma1 = mu_gamma1 + bg[,3];
gamma2 = mu_gamma2 + bg[,4];

// intercept and slope for each species/site
logit_p = alpha[st_id] + beta1[sp_id] .* PC1_1 + beta2[sp_id] .* PC1_2 + gamma1[sp_id] .* PC2_1 + gamma2[sp_id] .* PC2_2 + kappa[sp_id] .* EF + epsilon;
}

model {
// priors
sigma_eps ~ normal(mu_se, sd_se);
mu_se ~ normal(0, 3);
sd_se ~ normal(0, 3);
nu ~ gamma(2, 0.1);
mma ~ std_normal();
sma ~ normal(0, 3);
mu_beta1 ~ normal(0, 6);
mu_beta2 ~ normal(0, 6);
mu_gamma1 ~ normal(0, 6);
mu_gamma2 ~ normal(0, 6);
sigma_alpha ~ normal(0, 3);
alpha_raw ~ std_normal();
sigma_bg ~ normal(0, 5);
mu_kappa ~ normal(0, 3);
sigma_kappa ~ normal(0, 3);

// centered
epsilon ~ student_t(nu, 0, sigma_eps[sp_id]);     // overdispersion parameter
mu_alpha ~ normal(mma, sma);                      // species-level intercept drawn from distr

// non-centered
kappa ~ normal(mu_kappa, sigma_kappa);

to_vector(z_bg) ~ std_normal();
L_Rho_bg ~ lkj_corr_cholesky(2);

y ~ binomial_logit(ntot, logit_p);
}

generated quantities {
vector[N] p;
vector[N] lpnre;

p = inv_logit(logit_p);
lpnre = alpha[st_id] + beta1[sp_id] .* PC1_1 + beta2[sp_id] .* PC1_2 + gamma1[sp_id] .* PC2_1 + gamma2[sp_id] .* PC2_2;

}


