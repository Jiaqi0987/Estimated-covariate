library(nprobust)
library(splines)
library(splines2)
library(gbm)
library(caret)
library(glmnet)
library(mvtnorm)
bc.estimator <- function(y, rhat, a, muhat, eval, hseq, deg = 1, g_type = "null"){
  n <- length(y)
  get.QR <- function(y, rhat, a, muhat, t.val, h){
    u <- (rhat - t.val) / h
    K <- dnorm(u)
    K.der  <- -u * K / h
    U <- sapply(0:deg, function(k) u^k / factorial(k))
    U.der <- cbind(0, U[, (1 <= deg) * (1:deg)]) / (h)
    aug.term <- (K.der * U + K * U.der) * (a - rhat)
    eps1_mat <- K * U * y + aug.term * muhat   # n x (deg+1), per-obs eps1 contributions
    Rhat <- colMeans(eps1_mat)
    t11 <- mean(K.der * (a - rhat) + K)
    t12 <- t21 <- mean((K.der * u + K / h) * (a - rhat) + K * u)
    t22 <- mean((K.der * u^2 + K * 2*u / h) * (a - rhat) + K * u^2)
    Qhat <- matrix(c(t11, t12,
                     t21, t22), ncol = 2, nrow = 2, byrow = TRUE)
    return(list(Qhat=Qhat, Rhat=Rhat, eps1_mat=eps1_mat,
                K=K, K.der=K.der, u=u, U=U, U.der=U.der, aug.term=aug.term))
  }
  #  use cv to choose best bandwidth h for BC
    folds <- createFolds(y, k = 5, list = TRUE)
    cv.error = rep(0,length(hseq))
    for (i in 1:length(hseq)){
    h_tmp = hseq[i]
    fold_errors <- sapply(folds, function(test.idx) {
      train.idx <- setdiff(seq_len(length(y)), test.idx)

      rhat_te = rhat[test.idx]
      y_te = y[test.idx]

      est <- rep(NA, length(eval))
        for(j in 1:length(eval)) {
          t.val <- eval[j]
          tmp <- get.QR(y[train.idx], rhat[train.idx], a[train.idx], muhat[train.idx], t.val, h_tmp)
          Qhat <- tmp$Qhat
          Rhat <- tmp$Rhat
          est[j] <- (solve(Qhat) %*% Rhat)[1, 1]
        }

      mean((approx(eval, est, xout = rhat_te)$y- y_te)^2,na.rm=TRUE)
    })
    cv.error[i] = mean(fold_errors)
  }
  idx_first <- which.min(cv.error)
  best.h <- hseq[idx_first]
  est <- rep(NA, length(eval))
  se  <- rep(NA, length(eval))
  for(j in 1:length(eval)) {
    t.val <- eval[j]
    tmp   <- get.QR(y, rhat, a, muhat, t.val, best.h)
    Rhat  <- tmp$Rhat
    K     <- tmp$K;  K.der <- tmp$K.der
    u_    <- tmp$u;  U_ <- tmp$U;  U.der_ <- tmp$U.der

    # bhat: exact BC estimator solution to (P_n eps2) bhat = P_n eps1
    bhat  <- solve(tmp$Qhat) %*% Rhat
    est[j] <- bhat[1, 1]

    # Q_gram = kernel Gram matrix E[K_{ht}(r) g g^T], used only for the SE sandwich
    Qhat_gram <- matrix(c(mean(K),      mean(K * u_),
                          mean(K * u_), mean(K * u_^2)),
                        nrow = 2, byrow = TRUE)

    # psi_i does not use the known simulation g_type. Always include the
    # feasible rhohat term, since in applications we do not know whether rho=0.
    # where b1_i = m̂_t(rhat_i), b2_i = m̂'_t, rhohat_i = muhat_i - b1_i
    b1      <- as.vector(U_ %*% bhat)       # m̂_t(rhat_i)
    b2      <- as.vector(U.der_ %*% bhat)   # m̂'_t (same for all i at fixed t)
    rhohat  <- muhat - b1                   # feasible proxy for varpi_i
    psi     <- K * U_ * (y - b1 - b2 * (a - rhat))
    psi     <- psi + tmp$aug.term * rhohat

    Qinv  <- solve(Qhat_gram)
    Omega <- crossprod(psi) / n
    q1    <- Qinv[1, , drop = FALSE]
    se[j] <- sqrt(as.numeric(q1 %*% Omega %*% t(q1)) / n)
  }
  return(list(est=est, df.or.h = best.h, se = se))
}

basis_fun = function(x, deg, basis = "polynomial", derivative = FALSE){
  
  if(basis == "fourier"){
    if(derivative ==FALSE){
      mat <- matrix(0, nrow = length(x), ncol = deg)
      for (k in 1:deg) {
        mat[, k] <- cos(pi * k * x)
      }
      return(mat)
    } else{
      mat <- matrix(0, nrow = length(x), ncol = deg)
      for (k in 1:deg) {
        mat[, k] <- -pi * k * sin(pi * k * x)
      }
      return(mat)
    }
  } else if (basis == "bspline"){
    library(splines)
    library(splines2)
    if(derivative==FALSE){
      B <- bs(x, df = deg, degree = min(deg, 3), Boundary.knots = c(0, 1))
      return(B)
    } else{
      return(splines2::dbs(x, derivs=1, df = deg, degree = min(deg, 3),
                           Boundary.knots = c(0, 1)))
    }
  }else if (basis == "polynomial"){
    if (!derivative) {
      mat = poly(x, degree = deg, raw = TRUE)
    } else {
      if(deg > 1) {
        der.coef <- matrix(2:deg, nrow=length(x), ncol=deg-1, byrow=TRUE)
        mat <- cbind(1, poly(x, degree=deg-1, raw=TRUE) * der.coef)
      } else {
        mat <- matrix(1, ncol=1, nrow=length(x))
      }
      
    }
    return(mat)
  }
}

cpi.spline <- function(y, a, rhat, der.est, eval, df.seq, 
                       basis = "polynomial", do.gam = FALSE) {
  
  if(do.gam) {
    # dat.tmp <- data.frame(y=y-der.est*(a-rhat), rhat=rhat)
    # fit <- mgcv::gam(y ~ s(rhat, k=20), data=dat.tmp)
    dat.tmp <- data.frame(y = y, rhat = rhat)
    fit.tmp <- mgcv::gam(y ~ s(rhat, k=20), data=dat.tmp)
    X_full <- predict(fit.tmp, type = "lpmatrix")
    fit.star <- glm(a ~ -1 + offset(logit(rhat)) + ., 
                    data = cbind(data.frame(a = a), der.est * X_full),
                    family = binomial())
    rhat_star <- predict(fit.star, type = "response")
    dat.tmp2 <- data.frame(y = y, rhat = rhat_star)
    fit <- mgcv::gam(y ~ s(rhat, k=20), data=dat.tmp2)
    
    est <- as.vector(predict(fit, newdata = data.frame(rhat=eval)))
    return(list(est=est, fit=fit))
  } else {

    get_est <- function(y, r_curr, a, eval, df, basis, knots = NA){
      d.m <- cbind(1, basis_fun(r_curr, df, basis = basis))
      mdl = lm(y~-1+d.m)
      bhat = coef(mdl)
      sigma_hat = summary(mdl)$sigma
      if(basis =="bspline"){
        eval.basis <- cbind(1, bs(eval, df = df, knots = knots, degree=min(df, 3),
                                  Boundary.knots = c(0,1)))
      }else{
        eval.basis = cbind(1, basis_fun(eval, df, basis = basis))
      }
      est = eval.basis %*% bhat
    
      return(list(est=est))
    }
    #plot(r,rhat)
    #abline(a = 0, b =1  , col = "red")
    #plot(r,r_curr)
    #abline(a = 0, b= 1,col = "red")
    
    n <- length(y)
    fold_id <- rep(1:5, length.out = n)
    cv_err <- rep(NA, length(df.seq))
    for ( i in 1:length(df.seq)){
      df = df.seq[i]
      fold_err = rep(NA,5)
      X = bs(rhat, df = df, Boundary.knots = c(0,1), degree = min(3,df))
      knots   <- attr(X, "knots")
      for(k in 1:5){
        idx_val <- which(fold_id == k)
        idx_tr  <- setdiff(seq_len(n), idx_val)
        y_tr <- y[idx_tr]; y_val <- y[idx_val]
        rhat_tr <- rhat[idx_tr]; rhat_val <- rhat[idx_val]
        a_tr <- a[idx_tr]
        pred = get_est(y = y_tr,r_curr = rhat_tr,a = a_tr, eval = rhat_val, df =df,basis  = basis, knots = knots )$est
        fold_err[k] = mean((y_val - pred)^2)
      }
      cv_err[i] <- mean(fold_err)
    }
    
    #df <- min(df.seq[cv_err/min(cv_err) < 1.01])
    df <-min(df.seq[which.min(cv_err)])
    if(basis =="bspline"){
      X = bs(rhat, df = df, Boundary.knots = c(0,1), degree = min(3,df))
      knots   <- attr(X, "knots")
    }
    tmp = get_est(y,rhat,a, eval,df,basis, knots = knots)
    est = tmp$est
    
    # compute adjusted variance
    tmp = get_est(y = y, r_curr = rhat, a = a,eval = rhat, df = df, basis = basis, knots = knots)
    m_rhat = tmp$est
    
    if(basis =="bspline"){
      Phi_eval <- cbind(1, bs(eval, df = df, knots = knots, degree=min(df, 3),
                                Boundary.knots = c(0,1)))
    }else{
      Phi_eval <- cbind(1, basis_fun(eval, df, basis = basis))
    }
    if(basis == "bspline"){
      Phi_rhat <- cbind(1, bs(rhat, df = df, knots = knots, degree = min(df, 3),
                              Boundary.knots = c(0, 1)))
    } else {
      Phi_rhat <- cbind(1, basis_fun(rhat, df, basis = basis))
    }
    Qhat <- crossprod(Phi_rhat) / n
    Qinv <- tryCatch(solve(Qhat),
                     error = function(e) solve(Qhat + 1e-8 * diag(nrow(Qhat))))
    # y is already the corrected outcome (y_orig - der.est*(a-rhat)), so residual is just y - m_rhat
    u <- as.numeric(y - m_rhat)
    # robust Omega = Pn[ phi phi^T u^2 ]  (implemented as (Phi * u)'(Phi * u) / n
    Omega_hat <- crossprod(Phi_rhat * u) / length(u)
    Vbeta_hat <- Qinv %*% Omega_hat %*% Qinv
    se.adjust <- sqrt(pmax(diag(Phi_eval %*% Vbeta_hat %*% t(Phi_eval)), 0) / n)


    return(list(est = est, df.or.h = df,se.adj = se.adjust))
  }
}

## BC for splines

bc.spline= function (y, rhat, a, muhat, eval, df.seq, basis, g_type = "null"){
  n = length(y)
  inv_Qhat <- function(Qhat, ridge = 1e-8) {
    out <- try(solve(Qhat), silent = TRUE)

    if (inherits(out, "try-error")) {
      # fallback: ridge-stabilized inverse
      out <- solve(Qhat + ridge * diag(ncol(Qhat)))
    }

    out
  }
  get_est = function(y,rhat,a,muhat,eval,df,basis,knots = NA){
    d.m = cbind(1,basis_fun(rhat, df, basis = basis))
    d.m.der = cbind(0,basis_fun(rhat, df, basis = basis, derivative = TRUE))
    Qhat <- t(d.m)%*%d.m
    Rhat <- t(d.m)%*%y
    What <- diag(a - rhat)
    dotQ <- t(d.m)%*%What%*%d.m.der + t(d.m.der)%*%What%*%d.m # IF(Q)
    Qinv <- inv_Qhat(Qhat)
    Qder <- -Qinv %*%dotQ%*%Qinv
    Rtilde <- t(d.m)%*%y + t(d.m.der)%*%((a-rhat)*muhat)
    bhat <- Qinv %*%Rtilde + Qder%*%Rhat # bias-corrected betahat

    if(basis =="bspline"){
      eval.basis <- cbind(1,bs(eval, df = df, knots = knots, Boundary.knots = c(0,1), degree = min(3,df)))
    }else{
      eval.basis =cbind(1, basis_fun(eval, df, basis = basis))
    }

    est = eval.basis %*% bhat
    return(list(est=est, bhat=bhat, Qhat=Qhat, dotQ=dotQ,
                d.m=d.m, d.m.der=d.m.der, eval.basis=eval.basis))
  }
  if(length(df.seq)>1){


    fold_id <- sample(rep(1:5, length.out = length(y)))
    cv_err <- rep(NA, length(df.seq))
    for(i in 1:length(df.seq)){
      df = df.seq[i]
      fold_err = rep(NA,5)
      X = bs(rhat, df = df, Boundary.knots = c(0,1), degree = min(3,df))
      knots   <- attr(X, "knots")
      for(k in 1:5){
        idx_val <- which(fold_id == k)
        idx_tr  <- setdiff(seq_len(n), idx_val)
        y_tr <- y[idx_tr]; y_val <- y[idx_val]
        rhat_tr <- rhat[idx_tr]; rhat_val <- rhat[idx_val]
        a_tr <- a[idx_tr]
        muhat_tr <- muhat[idx_tr]
        pred = get_est(y = y_tr,rhat = rhat_tr,a = a_tr,muhat = muhat_tr,
                       eval = rhat_val, df =df,basis  = basis,knots = knots )
        fold_err[k] = mean((y_val - pred$est)^2)
      }
      cv_err[i] <- mean(fold_err)
    }
    df <-min(df.seq[which.min(cv_err)])
    #df = min(df.seq[cv_err/min(cv_err) < 1.01])
  }else{
    df = df.seq
  }
  if(basis == "bspline"){
    X = bs(rhat, df = df, Boundary.knots = c(0,1), degree = min(3,df))
    knots   <- attr(X, "knots")
  }
  res = get_est(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval, df = df, basis = basis, knots = knots)
  est        <- res$est
  bhat       <- res$bhat
  d.m        <- res$d.m
  d.m.der    <- res$d.m.der
  eval.basis <- res$eval.basis

  # Sandwich SE based on Proposition 3 linear representation.
  # Q = E[phi_k phi_k'] estimated by Gram matrix (matches proposition, not P_n eps2).
  # psi_i does not use the known simulation g_type. Always include the
  # feasible rhohat term, since in applications we do not know whether rho=0.
  mhat_rhat  <- as.vector(d.m     %*% bhat)   # m̂(r̂_i)
  dmhat_rhat <- as.vector(d.m.der %*% bhat)   # m̂'(r̂_i)
  rhohat     <- muhat - mhat_rhat             # feasible proxy for varpi_i

  psi <- d.m * (y - mhat_rhat - dmhat_rhat * (a - rhat))
  psi <- psi + d.m.der * ((a - rhat) * rhohat)

  Qhat_gram <- crossprod(d.m) / n            # Q_{phi,k}(r̂) = E[phi_k phi_k']
  Q_inv     <- inv_Qhat(Qhat_gram)
  Omega_hat <- crossprod(psi) / n
  V_beta    <- Q_inv %*% Omega_hat %*% Q_inv
  se <- sqrt(pmax(rowSums((eval.basis %*% V_beta) * eval.basis), 0) / n)

  return(list(est = est, df.or.h = df, se = se))
}
# train_rhat <- function(a, x, new.x, method = "lasso"){
#   if(method == "lasso") {
#     fit_rhat_cv <- cv.glmnet(
#       x = x,
#       y = a,
#       family = "binomial",
#       alpha = 1,
#       nfolds = 5,
#       standardize = TRUE,
#       intercept = TRUE
#     )
#     rhat <- as.numeric(predict(fit_rhat_cv, newx = new.x, s = "lambda.1se", type = "response"))
#     rhat_train <- as.numeric(predict(fit_rhat_cv, newx = x, s = "lambda.1se", type = "response"))
#   } else if(method == "gbm") {
#       dat.tmp <- cbind(data.frame(a = a), x)
#       colnames(dat.tmp) <- c("a", paste0("X", 1:ncol(x)))
#       new.x <- as.data.frame(new.x); colnames(new.x) <- paste0("X", 1:ncol(x))
#       # fit <- glm(a ~ ., data = dat.tmp, family=binomial())
#       fit <- gbm::gbm(a ~ ., data = dat.tmp)
#       rhat_train <- predict(fit, type = "response", newdata = dat.tmp)
#       rhat <- predict(fit, newdata = as.data.frame(new.x), type = "response")
#   } else if(method == "glm") {
#     dat.tmp <- cbind(data.frame(a = a), x)
#     colnames(dat.tmp) <- c("a", paste0("X", 1:ncol(x)))
#     new.x <- as.data.frame(new.x); colnames(new.x) <- paste0("X", 1:ncol(x))
#     fit <- glm(a ~ ., data = dat.tmp, family=binomial())
#     rhat_train <- predict(fit, type = "response", newdata = dat.tmp)
#     rhat <- predict(fit, newdata = as.data.frame(new.x), type = "response")
#   } else if(method == "rf") {
#     dat.tmp <- cbind(data.frame(a = a), x)
#     colnames(dat.tmp) <- c("a", paste0("X", 1:ncol(x)))
#     new.x <- as.data.frame(new.x); colnames(new.x) <- paste0("X", 1:ncol(x))
#     fit <- ranger::ranger(a ~ ., data = dat.tmp, probability = TRUE)
#     rhat_train <- predict(fit, type = "response", data = dat.tmp)$predictions[, 2]
#     rhat <- predict(fit, data = as.data.frame(new.x), type = "response")$predictions[, 2]
#   }
#   return(list(rhat=rhat, rhat_train=rhat_train))
# }
# 
# train_rhat <- function(a, x, new.x){
#   dat.tmp <- cbind(data.frame(a = a), x)
#   colnames(dat.tmp) <- c("a", paste0("X", 1:ncol(x)))
#   new.x <- as.data.frame(new.x); colnames(new.x) <- paste0("X", 1:ncol(x))
#   # fit <- glm(a ~ ., data = dat.tmp, family=binomial())
#   fit <- gbm::gbm(a ~ ., data = dat.tmp)
#   rhat_train <- predict(fit, type = "response", newdata = dat.tmp)
#   rhat <- predict(fit, newdata = as.data.frame(new.x), type = "response")
#   # rhat <- as.numeric(predict(fit_rhat_cv, newx = new.x, s = "lambda.1se", type = "response"))
#   # rhat_train <- as.numeric(predict(fit_rhat_cv, newx = x, s = "lambda.1se", type = "response"))
#   return(list(rhat=rhat, rhat_train=rhat_train))
# }

train_rhat <- function(a, x, new.x){
  dat.tmp <- cbind(data.frame(a = a), x)
  colnames(dat.tmp) <- c("a", paste0("X", 1:ncol(x)))
  #new.x <- as.data.frame(new.x); colnames(new.x) <- paste0("X", 1:ncol(x))
  fit_rhat_cv <- cv.glmnet(
    x = x,
    y = a,
    family = "binomial",
    alpha = 1,
    nfolds = 5,
    standardize = TRUE,
    intercept = TRUE
  )
  rhat <- as.numeric(predict(fit_rhat_cv, newx = new.x, s = "lambda.1se", type = "response"))
  rhat_train <- as.numeric(predict(fit_rhat_cv, newx = x, s = "lambda.1se", type = "response"))
  return(list(rhat=rhat, rhat_train=rhat_train))
}

train_muhat <- function(y, x, new.x){
#  fit_muhat_cv <- cv.glmnet(
#    x = x,
#    y = y,
#    family = "gaussian",
#    alpha = 1,          
#    nfolds = 5,
#    standardize = TRUE, 
#    intercept = TRUE    
#  )
#  muhat <- as.numeric(predict(fit_muhat_cv, newx = new.x, s = "lambda.1se", type = "response"))
#  muhat_train <- as.numeric(predict(fit_muhat_cv, newx = x, s = "lambda.1se", type = "response"))
#  return(list(muhat=muhat, muhat_train=muhat_train))
  library(ranger)
    fit <- ranger(
      y ~ .,
      data = data.frame(y = y, x),
      num.trees = 500,
      mtry = floor(ncol(x)/3),
      min.node.size = 5
    )
    
    muhat <- predict(fit, data = data.frame(new.x))$predictions
    muhat_train <- predict(fit, data = data.frame(x))$predictions
    
    return(list(muhat = muhat, muhat_train = muhat_train))
}
 
