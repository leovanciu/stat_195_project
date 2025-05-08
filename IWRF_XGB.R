library(data.table)
library(FSelectorRcpp)
library(ranger)
library(rsample)
library(rBayesianOptimization)
library(yardstick)
library(dplyr)
library(purrr)
library(ggplot2)
library(xgboost)

# settings
SUB_N    <- NULL
FOLDS    <- 10
IFOLD    <- 5
KFEAT    <- 9
TS       <- 800
TF       <- 500
B0       <- 10
BI       <- 30

grid <- expand.grid(a=seq(0,1,0.1), b=seq(0,1,0.1)) %>%
  filter(a + b <= 1) %>% mutate(c = 1 - a - b)

# load
d1 <- fread("statlog.csv")[, target := as.factor(target)]
d2 <- fread("S1Data.csv")[, .(
  TIME, Gender, Smoking, Diabetes, BP, Anaemia, Age,
  Ejection.Fraction, Sodium, Creatinine, Pletelets, CPK,
  target = as.factor(Event)
)]

if (!is.null(SUB_N)) {
  d1 <- d1[1:SUB_N]
  d2 <- d2[1:SUB_N]
}
sets <- list(stat=d1, hf=d2)

# helper
mm <- function(x) (x - min(x)) / (max(x)-min(x)+1e-12)
rank_fs <- function(df, y, w, k) {
  X <- df[, setdiff(names(df), y), with=FALSE]
  yv <- as.numeric(df[[y]])
  m0 <- colMeans(X[yv==1,])
  m1 <- colMeans(X[yv==2,])
  v0 <- apply(X[yv==1,],2,var); v1 <- apply(X[yv==2,],2,var)
  h  <- (m0-m1)^2 / (v0+v1+1e-6)
  mi <- information_gain(reformulate(".",y), df)$importance
  sd <- apply(X,2,sd)
  s  <- w[1]*mm(h) + w[2]*mm(mi) + w[3]*mm(sd)
  names(sort(s, decreasing=TRUE))[1:k]
}

# get global features
get_feats <- function(dat) {
  out <- vfold_cv(dat, v=FOLDS, strata=target)
  allf <- map(out$splits, function(sp) {
    tr <- analysis(sp)
    score1 <- function(a,b,c) {
      f <- rank_fs(tr, "target", c(a,b,c), KFEAT)
      inner <- vfold_cv(tr, v=IFOLD, strata=target)
      map_dbl(inner$splits, function(sp2) {
        tt <- analysis(sp2); vv <- assessment(sp2)
        m <- ranger(target~., tt[,c(f,"target"),with=FALSE],
                    probability=TRUE, num.trees=TS)
        p <- predict(m, vv[,f,with=FALSE])$predictions[,"1"]
        t0 <- factor(as.integer(vv$target)-1, 0:1)
        p0 <- factor(as.integer(p>0.5),      0:1)
        f_meas_vec(t0, p0, pos="1")
      }) %>% mean()
    }
    sc <- pmap_dbl(grid, score1)
    best <- grid[which.max(sc), ]
    rank_fs(tr, "target", c(best$a, best$b, best$c), KFEAT)
  })
  Reduce(intersect, allf)
}

# main runner
runit <- function(dat, nm, feats) {
  n   <- nrow(dat)
  oof <- numeric(n)
  out <- vfold_cv(dat, v=FOLDS, strata=target)
  
  for (i in seq_along(out$splits)) {
    tr <- analysis(out$splits[[i]])
    va <- assessment(out$splits[[i]])
    
    # BO for IWRF
    fun2 <- function(am, an, pr, mr, mn) {
      yv <- as.numeric(tr$target)-1
      pos <- which(yv==1); neg <- which(yv==0)
      nm1 <- floor(pr*nrow(tr))
      idx <- c(sample(pos,nm1,TRUE), sample(neg,nrow(tr)-nm1,TRUE))
      b0 <- tr[idx]
      inn <- vfold_cv(b0, v=IFOLD, strata=target)
      f1s <- map_dbl(inn$splits, function(sp2) {
        t1 <- analysis(sp2); v1 <- assessment(sp2)
        y1 <- as.numeric(t1$target)-1
        w0 <- am*nrow(t1)/(2*sum(y1==0))
        w1 <- an*nrow(t1)/(2*sum(y1==1))
        w  <- ifelse(t1$target==0, w0, w1)
        m  <- ranger(target~., t1[,c(feats,"target"),with=FALSE],
                     probability=TRUE, num.trees=TS,
                     mtry=floor(mr), min.node.size=floor(mn),
                     case.weights=w)
        p1 <- predict(m, v1[,feats,with=FALSE])$predictions[,"1"]
        t0 <- factor(as.integer(v1$target)-1, 0:1)
        p0 <- factor(as.integer(p1>0.5),     0:1)
        f_meas_vec(t0, p0, pos="1")
      })
      list(Score=mean(f1s), Pred=0)
    }
    
    bo <- BayesianOptimization(
      FUN        = fun2,
      bounds     = list(am=c(0.5,1.5), an=c(1,3),
                        pr=c(0.3,0.7),
                        mr=c(3,length(feats)),
                        mn=c(1,15)),
      init_points=B0, n_iter=BI,
      acq        ="ucb", kappa=2.576, verbose=FALSE
    )
    best <- bo$Best_Par
    
    # final ensemble
    yv  <- as.numeric(tr$target)-1
    pos <- which(yv==1); neg <- which(yv==0)
    nm1 <- floor(best[["pr"]]*nrow(tr))
    idx <- c(sample(pos,nm1,TRUE), sample(neg,nrow(tr)-nm1,TRUE))
    b0  <- tr[idx]
    w0  <- best[["am"]]*nrow(b0)/(2*sum(b0$target==0))
    w1  <- best[["an"]]*nrow(b0)/(2*sum(b0$target==1))
    cw  <- ifelse(b0$target==0, w0, w1)
    
    # RF
    rf <- ranger(target~., b0[,c(feats,"target"),with=FALSE],
                 probability=TRUE, num.trees=TF,
                 mtry=floor(best[["mr"]]),
                 min.node.size=floor(best[["mn"]]),
                 case.weights=cw)
    prf <- predict(rf, va[,feats,with=FALSE])$predictions[,"1"]
    
    # XGB
    db <- xgb.DMatrix(data=as.matrix(b0[,feats,with=FALSE]),
                      label=as.numeric(b0$target)-1)
    xg <- xgb.train(data=db, nrounds=200,
                    params=list(
                      objective   ="binary:logistic",
                      eval_metric ="logloss",
                      eta         =0.05, max_depth=3,
                      subsample   =0.8, colsample_bytree=0.8
                    ), verbose=0)
    pg <- predict(xg, as.matrix(va[,feats,with=FALSE]))
    
    p0 <- (prf + pg)/2
    idx <- setdiff(seq_len(n), out$splits[[i]]$in_id)
    oof[idx] <- p0
    
    # report
    t0 <- factor(as.integer(va$target)-1, 0:1)
    p0f <- factor(as.integer(p0>0.5),     0:1)
    m <- tibble(
      fold= i,
      acc = accuracy_vec(t0,p0f),
      pre = precision_vec(t0,p0f, pos="1"),
      rec = recall_vec(t0,p0f, pos="1"),
      f1  = f_meas_vec(t0,p0f, pos="1"),
      auc = roc_auc_vec(t0, p0)
    )
    print(m)
  }
  oof
}

# run
res <- imap(sets, ~{
  f <- get_feats(.x)
  runit(.x, .y, f)
})