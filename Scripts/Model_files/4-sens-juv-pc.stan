// juv_sens ~ pc + phylo

data {
int<lower=0> N;                   // number of obs
real juv_sens[N];
real<lower=0> juv_sens_sd[N];
vector[N] cov;                    // distance
matrix[N, N] R;                   // scaled phylo distance matrix (from phylogeny)
matrix[N, N] I;                   // identity matrix
}

parameters {
real<lower=0> sigma;
real alpha;
real beta;
real<lower=0, upper=1> lambda;
vector[N] theta;
}

transformed parameters {
vector[N] mu;
vector[N] mu_js;
matrix[N, N] R_lambda;
matrix[N, N] S;
matrix[N, N] L;

// Pagel's lambda
R_lambda = lambda * R + (1 - lambda) * I;
S = R_lambda * sigma;

mu = alpha + beta * cov;

// non-centered multi-normal
L = cholesky_decompose(S);
mu_js = mu + L * theta;
}

model {
alpha ~ normal(0, 5);
beta ~ normal(0, 5);
sigma ~ normal(0, 5);
lambda ~ uniform(0, 1);
theta ~ std_normal();

// response
juv_sens ~ normal(mu_js, juv_sens_sd);
}

generated quantities {
real y_rep[N];

y_rep = normal_rng(mu_js, juv_sens_sd);
}
