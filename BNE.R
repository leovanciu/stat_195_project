library(data.table)
library(purrr)
library(rstan); options(mc.cores = 4); rstan_options(auto_write = TRUE)
library(posterior)


logit_clip <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}


dist_mat <- function(x) as.matrix(dist(x))

obj   <- readRDS("iwrf_baseline_and_pcs.rds")
base  <- obj$baseline
pcs   <- obj$pcs

# Stan code
bne_code <- "
functions {
  matrix m32(int N, matrix D, real rho) {
    matrix[N,N] K;
    for (i in 1:N) for (j in i:N) {
      real r   = sqrt(3)*D[i,j]/rho;
      real v   = (1 + r)*exp(-r);
      K[i,j]=v; if (i!=j) K[j,i]=v;
    }
    return K;
  }
  vector mono(vector a) {
    int M = num_elements(a);
    vector[M] out;
    vector[M] inc = log1p_exp(a);
    out[1] = inc[1];
    for (m in 2:M) out[m] = out[m-1] + inc[m];
    return out / out[M];
  }
}
data {
  int<lower=1> N;               // num obs
  int<lower=0,upper=1> y[N];    // labels
  int<lower=1>     Kb;          // # base learners = 2
  matrix[N,Kb]     F;           // logits from RF + intercept
  matrix[N,N]      D;           // dist in PC space
  int<lower=5>     M;           // calibrator grid size
  vector[M]        ug;          // [0,1] grid
  matrix[M,M]      Dg;          // grid dist
}
parameters {
  simplex[Kb]  w;      // ensemble weights
  real<lower=0> s1;    real<lower=0> r1;  vector[N]  z1;
  real<lower=0> s2;    real<lower=0> r2;  vector[M]  z2;
  vector[M]      f0;                // for monotone
}
transformed parameters {
  matrix[N,N] K1 = s1^2 * m32(N,D,r1) + 1e-6*diag_matrix(rep_vector(1,N));
  vector[N]   d1 = cholesky_decompose(K1) * z1;
  vector[N]   mu = F * w + d1;
  vector[N]   p1 = inv_logit(mu);

  matrix[M,M] K2 = s2^2 * m32(M,Dg,r2) + 1e-6*diag_matrix(rep_vector(1,M));
  vector[M]   g  = mono(cholesky_decompose(K2) * z2 + f0);

  vector[N] p2;
  for (i in 1:N) {
    real u = p1[i];
    int j = 1;
    while (j < M && u > ug[j+1]) j++;
    if (j == M) p2[i] = g[M];
    else {
      real t = (u - ug[j]) / (ug[j+1] - ug[j]);
      p2[i] = (1-t)*g[j] + t*g[j+1];
    }
  }
}
model {
  // priors
  s1 ~ normal(0,1);  r1 ~ lognormal(0,1);
  s2 ~ normal(0,1);  r2 ~ lognormal(0,1);
  z1 ~ normal(0,1);  z2 ~ normal(0,1);
  f0 ~ normal(0,1);

  y ~ bernoulli(p2);
}
generated quantities {
  vector[N] p_hat = p2;
}
"

sm <- stan_model(model_code = bne_code)
M   <- 41
ug  <- seq(0.01, 0.99, length.out = M)
Dg  <- dist_mat(matrix(ug, ncol=1))

fit_one <- function(name) {
  cat("-> fitting BNE for", name, "...\n")
  p0   <- base[[name]]$p
  y0   <- base[[name]]$truth
  N    <- length(p0)
  Fmat <- cbind(1, logit_clip(p0))
  Dmat <- dist_mat(pcs[[name]])
  
  dat <- list(
    N  = N, y = y0,
    Kb = ncol(Fmat), F = Fmat,
    D  = Dmat, M = M, ug = ug, Dg = Dg
  )
  sampling(sm, data = dat,
           chains = 4, cores = 4,
           iter = 1000, warmup = 1000,
           control = list(adapt_delta=0.9, max_treedepth=12))
}

fits <- map(names(base), fit_one)

# collect summaries
out <- map2_dfr(fits, names(base), function(f, nm) {
  draws <- as_draws_df(f) |> select(starts_with("p_hat")) |> as.matrix()
  n     <- ncol(draws)
  data.table(
    ds   = nm,
    id   = 1:n,
    est  = colMeans(draws),
    lo   = apply(draws,2,quantile,0.025),
    hi   = apply(draws,2,quantile,0.975),
    truth= base[[nm]]$truth
  )
})

saveRDS(list(fits=fits, cal=out), "bne_results.rds")