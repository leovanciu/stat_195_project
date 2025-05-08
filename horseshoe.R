
library(data.table)
library(cmdstanr)
library(posterior)
library(rsample)
library(pROC)
library(progress)
library(purrr)

logit_clip <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

scale_tr <- function(dt) {
  nums <- names(dt)[sapply(dt, is.numeric) & names(dt) != "target"]
  sd0  <- sapply(dt[, ..nums], sd)
  keep <- nums[sd0 > 0]
  mns  <- sapply(dt[, ..keep], mean)
  sds  <- sd0[sd0 > 0]
  
  dt2 <- copy(dt)
  dt2[, (keep) := Map(function(col, m, s) (col - m) / s,
                      .SD, mns, sds), .SDcols = keep]
  list(data = dt2, cols = keep, m = mns, s = sds)
}

scale_te <- function(dt, cols, m, s) {
  dt2 <- copy(dt)
  dt2[, (cols) := Map(function(col, mm, ss) (col - mm) / ss,
                      .SD, m, s), .SDcols = cols]
  dt2
}

one_metrics <- function(y, p, thr = 0.5) {
  pr <- as.integer(p >= thr)
  tp <- sum(pr == 1 & y == 1)
  tn <- sum(pr == 0 & y == 0)
  fp <- sum(pr == 1 & y == 0)
  fn <- sum(pr == 0 & y == 1)
  
  acc  <- (tp + tn) / length(y)
  sens <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
  spec <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
  f1   <- if (tp + fp + fn > 0) 2*tp / (2*tp + fp + fn) else NA_real_
  
  mccd <- (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)
  mcc  <- if (!is.na(mccd) && mccd>0) (tp*tn - fp*fn)/sqrt(mccd) else NA_real_
  
  auc  <- as.numeric(roc(y, p, quiet=TRUE)$auc)
  data.table(acc, sens, spec, f1, mcc, auc)
}

# Datasets
stat <- fread("statlog.csv")[, target := as.integer(target)]
hf   <- fread("S1Data.csv")[, target := as.integer(Event)][, Event := NULL]
hf   <- setcolorder(hf, c(setdiff(names(hf), "target"), "target"))
all_ds <- list(statlog = stat, hfail = hf)

# Stan model
stan_code <- "
data { int<lower=1> N, K; matrix[N,K] X; int<lower=0,upper=1> y[N];
       real mu0; real<lower=0> t0; }
parameters { real a; vector[K] z; vector<lower=0>[K] l; real<lower=0> t; }
transformed parameters { vector[K] b = t * l .* z; }
model {
  a ~ normal(mu0,2);
  z ~ normal(0,1);
  l ~ cauchy(0,1);
  t ~ cauchy(0,t0);
  y ~ bernoulli_logit(a + X * b);
}
generated quantities {
  vector[N] th;
  for (i in 1:N) th[i] = inv_logit(a + X[i] * b);
}
"
mod <- cmdstan_model(write_stan_file(stan_code),
                     cpp_options = list(stan_threads = TRUE))

# loop over datasets
for (nm in names(all_ds)) {
  cat("\n---", nm, "---\n")
  dt   <- all_ds[[nm]]
  N    <- nrow(dt)
  folds <- vfold_cv(dt, v = 10, strata = "target")$splits
  
  o_mean <- numeric(N)
  o_lo   <- numeric(N)
  o_hi   <- numeric(N)
  mets   <- list()
  pb     <- progress_bar$new(total = length(folds),
                             format = " Fold :current/:total [:bar] :eta")
  
  for (sp in folds) {
    pb$tick()
    ti <- sp$in_id
    te <- setdiff(seq_len(N), ti)
    
    dtr <- dt[ti]; dte <- dt[te]
    sc  <- scale_tr(dtr)
    dtr <- sc$data
    dte <- scale_te(dte, sc$cols, sc$m, sc$s)
    
    Xtr <- as.matrix(dtr[, sc$cols, with=FALSE])
    Xte <- as.matrix(dte[, sc$cols, with=FALSE])
    ytr <- dtr$target; yte <- dte$target
    
    mu0 <- qlogis((sum(ytr)+0.5)/(length(ytr)+1))
    sd0 <- 1
    sd_dat <- list(N = nrow(Xtr), K = ncol(Xtr),
                   X = Xtr, y = ytr,
                   mu0 = mu0, t0 = sd0)
    
    fit <- mod$sample(data = sd_dat,
                      chains = 4, parallel_chains = 2,
                      threads_per_chain = 2,
                      iter_warmup = 1000, iter_sampling = 1000,
                      adapt_delta = 0.9, refresh = 0)
    
    draws <- as_draws_matrix(fit$draws(c("a", "b")))
    a0     <- draws[, "a"]
    bm     <- draws[, grep("^b\\[", colnames(draws)), drop=FALSE]
    lin    <- Xte %*% t(bm) + matrix(a0, nrow(Xte), nrow(bm), byrow=TRUE)
    th_samp <- plogis(lin)
    
    o_mean[te] <- rowMeans(th_samp)
    o_lo  [te] <- apply(th_samp, 1, quantile, 0.025)
    o_hi  [te] <- apply(th_samp, 1, quantile, 0.975)
    
    mets[[length(mets)+1]] <-
      map_dfr(seq_len(ncol(th_samp)), function(j) {
        one_metrics(yte, th_samp[,j])
      })
  }
  
  # save predictions + CIs
  fwrite(data.table(
    id    = 1:N,
    truth = dt$target,
    mean  = o_mean,
    lo95  = o_lo,
    hi95  = o_hi
  ), paste0(nm, "_oof.csv"))
  
  allm <- rbindlist(mets)
  summ <- allm[, lapply(.SD, function(x) list(
    avg = mean(x), lo = quantile(x, 0.025), hi = quantile(x, 0.975)
  ))]
  fwrite(data.table(t(unlist(summ))), paste0(nm, "_metrics.csv"))
}

