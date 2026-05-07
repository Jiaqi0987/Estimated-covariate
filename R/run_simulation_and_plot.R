#!/usr/bin/env Rscript
# Run simulations and/or generate plots for the estimated-covariates study.
# Output filename prefix: sim_semior_*

parse_args <- function(x) {
  out <- list(); i <- 1
  while (i <= length(x)) { out[[sub("^--", "", x[i])]] <- x[i + 1]; i <- i + 2 }
  out
}

`%||%` <- function(x, y) if (!is.null(x) && !is.na(x) && nzchar(x)) x else y

args <- parse_args(commandArgs(trailingOnly = TRUE))

mode <- args$mode %||% "run"
if (!(mode %in% c("run", "plot", "both"))) {
  stop("--mode must be one of: run, plot, both")
}

setting <- if (mode %in% c("run", "both")) as.integer(args$setting %||% stop("Need --setting")) else NA_integer_
nsim    <- as.integer(args$nsim    %||% "500")
d       <- as.integer(args$d       %||% "200")
outdir  <- args$outdir %||% "results"
rate    <- as.numeric(args$rate    %||% "0.3")
n       <- if (mode %in% c("run", "both")) as.integer(args$n %||% stop("Need --n")) else NA_integer_
n.tr    <- as.integer(args$n_tr %||% args$ntr %||% as.character(n))
g_type  <- if (mode %in% c("run", "both")) args$g_type %||% stop("Need --g_type") else NA_character_

if (mode %in% c("run", "both")) {
  if (!(setting %in% c(1, 2, 3))) stop("--setting must be one of 1, 2, 3")
  if (!(g_type %in% c("null", "linear", "interact"))) {
    stop("--g_type must be one of: null, linear, interact")
  }
}

eval.pts <- seq(0.15, 0.85, length.out = 50)
assign("eval.pts", eval.pts, envir = .GlobalEnv)

cores <- as.integer(Sys.getenv("SIM_CORES", unset = "1"))
if (is.na(cores) || cores < 1) cores <- 1
message("using mc.cores = ", cores)

cmdArgs <- commandArgs(trailingOnly = FALSE)
fileArg <- cmdArgs[grep("^--file=", cmdArgs)]
script_dir <- normalizePath(dirname(sub("^--file=", "", fileArg)))
setwd(script_dir)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

source("./estimators.R")
source("./plot_functions.R")

# ---- simulation engine ------------------------------------------------------
# Embedded here so the share folder has only one runnable script.
library(nprobust)
`%||%` <- function(x, y) if (!is.null(x)) x else y
library(splines)
library(splines2)
library(gbm)
library(caret)
library(glmnet)
library(mvtnorm)
library(mgcv)
library(future.apply)

expit <- function(u) 1/(exp(-u) + 1)
logit <- function(u) log(u/(1-u))

# ============================================================
# Semi-oracle helpers
# ------------------------------------------------------------
# For each real estimator, refit a "semi-oracle" version using the
# TRUE propensity r and the realized smoothing parameter (df for
# sieves, h for lpoly), with NO bias / rho correction. The semi-oracle
# answers: "what would a clean OR fit have produced at this estimator's
# smoothing level?" The downstream coverage check uses the
# across-replicate mean of these semi-oracle estimates as the target.
# ============================================================
semi_or_sieve <- function(y, r, a, eval.pts, df, basis) {
  fit <- cpi.spline(y = y, rhat = r, a = a, eval = eval.pts,
                    der.est = 0, df.seq = df, basis = basis)
  list(est = as.numeric(fit$est),
       se  = as.numeric(fit$se.adj %||% fit$se),
       df.or.h = df)
}

semi_or_lpoly <- function(y, r, eval.pts, h, poly_deg = 1) {
  fit <- lprobust(y = y, x = r, p = poly_deg, eval = eval.pts,
                  h = h, kernel = "gau")
  list(est = as.numeric(fit$Estimate[, "tau.us"]),
       se  = as.numeric(fit$Estimate[, "se.us"]),
       df.or.h = h)
}

# Pull the realized df / h from an estimator's fitted object.
# Sieves: cpi.spline / bc.spline store df as $df.or.h.
# lpoly via lprobust (PI/CPI/etc.): mean(Estimate[,"h"]) is what gets stored.
# lpoly via bc.estimator (BC): $df.or.h is best.h (scalar).
get_dfh <- function(obj, basis) {
  if (basis == "lpoly") {
    if (!is.null(obj$df.or.h)) return(as.numeric(obj$df.or.h))
    if (!is.null(obj$Estimate)) return(mean(obj$Estimate[, "h"]))
  }
  as.numeric(obj$df.or.h %||% obj$best.df)
}

build_semi_or_block <- function(y, r, a, eval.pts, poly_deg,
                                fits) {
  # `fits` is a named list of list(fourier=, bs=, poly=, lpoly=) blocks,
  # one per estimator. Returns the same shape with semi-oracle results.
  bases_sieve <- c(fourier = "fourier", bs = "bspline", poly = "polynomial")
  out <- lapply(fits, function(blk) {
    res <- list()
    for (b in names(bases_sieve)) {
      df_b <- get_dfh(blk[[b]], b)
      res[[b]] <- semi_or_sieve(y, r, a, eval.pts, df_b, bases_sieve[[b]])
    }
    h <- get_dfh(blk$lpoly, "lpoly")
    res$lpoly <- semi_or_lpoly(y, r, eval.pts, h, poly_deg)
    res
  })
  out
}

# functions to compute the true E(g(x)|r(X) = t) and its derivative when \rho!= 0
truth_fun <- function(u, g_type, beta_r, beta_g, Sigma) {

  if (g_type == "interact"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    v   <- as.vector(Sigma %*% beta_r)                  # Σβ
    v12 <- v[1] * v[2]
    c12 <- Sigma[1, 2]
    s   <- logit(u)
   beta_g[1]*(c12 - v12/den + (v12/den^2) * s^2)
  }else if (g_type =="linear"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)  
    lambda <- as.numeric(t(beta_g) %*% Sigma %*% beta_r / den)
    lambda * logit(u)
  }else if (g_type =="null"){
    0
  }
}

deriv_fun <- function(t, g_type, beta_r, beta_g, Sigma){
  if(g_type =="interact"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    v   <- as.vector(Sigma %*% beta_r)
    v12 <- v[1] * v[2]
    s   <- logit(t)
    2 * beta_g[1]*(v12/den^2) * s * dlogit(t)
  }else if (g_type =="linear"){
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    lambda <- as.numeric(t(beta_g) %*% Sigma %*% beta_r / den)
    lambda * dlogit(t)
  }else if(g_type =="null"){
    0
  }
}


dlogit <- function(p) {
  if (any(p <= 0 | p >= 1)) {
    stop("p must be strictly between 0 and 1")
  }
  return(1 / (p * (1 - p)))
}

# d2logit = d/dt[1/(t(1-t))] = (2t-1) / (t(1-t))^2
d2logit <- function(t) (2*t - 1) / (t * (1 - t))^2

deriv2_fun <- function(t, g_type, beta_r, beta_g, Sigma) {
  if (g_type == "null") {
    0
  } else if (g_type == "linear") {
    den    <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    lambda <- as.numeric(t(beta_g) %*% Sigma %*% beta_r / den)
    lambda * d2logit(t)
  } else if (g_type == "interact") {
    den <- as.numeric(t(beta_r) %*% Sigma %*% beta_r)
    v   <- as.vector(Sigma %*% beta_r)
    v12 <- v[1] * v[2]
    # deriv_fun = 2*beta_g[1]*(v12/den^2) * logit(t) * dlogit(t)
    # deriv2   = 2*beta_g[1]*(v12/den^2) * (dlogit(t)^2 + logit(t)*d2logit(t))
    2 * beta_g[1] * (v12 / den^2) * (dlogit(t)^2 + logit(t) * d2logit(t))
  }
}

gen_data <- function(beta_r, beta_g, n, d, Sigma, r.x, mu.x) {
  x <- rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  r <- r.x(beta_r,x)
  a <- rbinom(n, p = r, size = 1)
  y <- rnorm(n, mean = mu.x(x), sd = 1)
  dat <- data.frame(y = y, a = a, r = r)
  return(list(dat = dat, x = x))
}

simu_rslt = function(seed, config){
  
  set.seed(seed)
  # extract everything from config
  n = config$n
  n.tr = config$n.tr
  df.seq.bs = config$df.seq.bs
  df.seq.poly = config$df.seq.poly
  df.seq.fourier = config$df.seq.fourier
  g_type = config$g_type
  der.m.t = config$der.m.t
  m.t = config$m.t
  r.x = config$r.x
  g.x = config$g.x
  mu.x = config$mu.x
  beta_g = config$beta_g
  beta_r = config$beta_r
  Sigma = config$Sigma
  d = config$d
  poly_deg = config$poly_deg
  eval.pts = config$eval.pts
  rate = config$rate
  est_r = config$est_r
  h.seq = config$h.seq
  fix.df.or.h = config$fix.df.or.h
  # generate auxiliary sample to estimate r.x of size 
  train_data <- gen_data(beta_r = beta_r, beta_g = beta_g, n = n.tr, d = d, Sigma = Sigma,
                          r.x = r.x, mu.x = mu.x)
  
  # generate sample for inference
  inference_data <- gen_data(beta_r = beta_r, beta_g = beta_g, n = n, d = d, Sigma = Sigma,
                             r.x = r.x, mu.x = mu.x)
  dat <- inference_data$dat
  y <- dat$y; r <- dat$r; a <- dat$a; x <- inference_data$x
  y.tr <- train_data$dat$y; a.tr <- train_data$dat$a; x.tr <- train_data$x
  # estimate rhat or generate rhat given the convergence rate
  if(est_r == TRUE){
    rhat_fit <- train_rhat(a=a.tr, x=x.tr, new.x=x)
    rhat <- rhat_fit$rhat
    rhat.tr <- rhat_fit$rhat_train
  }else{
    rhat <- expit(logit(r) + rnorm(n, mean =1 / (n)^rate, sd = 1 / (n)^rate)*sample(c(-1,1),size = n, replace = T))
    rhat.tr= expit(logit(train_data$dat$r) + rnorm(n.tr, mean = 1 / (n.tr)^rate, sd = 1 / (n.tr)^rate)*sample(c(-1,1),size = n.tr, replace = T))
  }
  
  #calibrated rhat
  rhat_star_fit = isoreg(y = a, x = rhat)
  rhat_star  = numeric(length(a))
  rhat_star[ order(rhat)] = rhat_star_fit$yf
  eps <- 1e-6
  rhat_cali <- pmin(pmax(rhat_star, eps), 1 - eps)

  # estimate derivative
  eval.pts2 <- seq(0, 1, length.out=150)
  tmp<- lprobust(y = train_data$dat$y, x = rhat.tr, eval = eval.pts2,
                      deriv = 1, kernel = "gau")$Estimate[, "tau.us"]
  #estimated derivative for rhat
  der.est <- approx(x = eval.pts2, y = tmp, xout = rhat,      rule = 2)$y
  #estimated derivative for calibrated rhat
  der.est_cali = approx(x = eval.pts2, y = tmp, xout = rhat_cali, rule = 2)$y
  # estimtate muhat
  muhat_fit <- train_muhat(y=y.tr, x=x.tr, new.x=x)
  muhat <- muhat_fit$muhat
  mu <- as.vector(mu.x(x))

  ## oracle estimator knowing r.x(x)
  or.fourier <-  cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                            der.est = 0, df.seq = df.seq.fourier, basis = "fourier")
  or.poly <-  cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                         der.est = 0, df.seq = df.seq.poly, basis = "polynomial")
  or.bs <-   cpi.spline(y = y, a = a, rhat = r, eval = eval.pts, 
                        der.est = 0, df.seq = df.seq.bs, basis = "bspline")
  or.lpoly <- lprobust(y = y, x = r, p = poly_deg, eval = eval.pts, 
                     kernel = "gau")

  ##PI
  pi.fourier <- cpi.spline(y = y, rhat = rhat, a = a, eval=eval.pts, der.est = 0,
                           df.seq = df.seq.fourier , basis = "fourier")
  pi.poly <- cpi.spline(y = y, rhat = rhat, a=a ,eval=eval.pts, der.est= 0,
                        df.seq = df.seq.poly , basis ="polynomial")
  pi.bs <-cpi.spline(y = y, rhat = rhat, a=a, eval=eval.pts, der.est = 0,
                     df.seq = df.seq.bs , basis = "bspline")
  pi.lpoly = lprobust(y = y, x = rhat, p = poly_deg, eval = eval.pts, 
                      kernel = "gau")

  if(fix.df.or.h == TRUE){
    df.seq.fourier = pi.fourier$best.df
    df.seq.bs = pi.bs$best.df
    df.seq.poly = pi.poly$best.df
    h.seq = mean(pi.lpoly$Estimate[,"h"])
  }

  ### Corrected plug-ins oracle (knowledge of derivative) ###
  f.der_p = der.m.t(rhat)
  cpi.or.fourier <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                               df.seq = df.seq.fourier , basis = "fourier")
  cpi.or.poly <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                            df.seq = df.seq.poly , basis = "polynomial")
  cpi.or.bs <- cpi.spline(y = y-f.der_p*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = f.der_p,
                          df.seq = df.seq.bs , basis = "bspline")
  cpi.or.lpoly <- lprobust(y = y-f.der_p*(a-rhat), x = rhat, p = poly_deg, eval = eval.pts, 
           kernel = "gau")
  
  ### Corrected plug-ins ###

  cpi.fourier <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                            df.seq = df.seq.fourier , basis = "fourier")
  cpi.poly <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                         df.seq = df.seq.poly , basis = "polynomial")
  cpi.bs <- cpi.spline(y = y-der.est*(a-rhat), rhat = rhat, a = a, eval=eval.pts, der.est = der.est,
                       df.seq = df.seq.bs , basis = "bspline")
  cpi.lpoly <- lprobust(y = y-der.est*(a-rhat), x = rhat, p = poly_deg, eval = eval.pts, 
                        kernel = "gau")
  
  ##PI with calibration (corrected outcome: same as CPI but using rhat_cali)
  y_cali <- y - der.est_cali * (a - rhat_cali)
  pi.cali.fourier <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                            df.seq = df.seq.fourier, basis = "fourier")
  pi.cali.poly <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                         df.seq = df.seq.poly, basis = "polynomial")
  pi.cali.bs <- cpi.spline(y = y_cali, rhat = rhat_cali, a = a, eval=eval.pts, der.est = der.est_cali,
                       df.seq = df.seq.bs, basis = "bspline")
  # first try automatic bandwidth
  pi.cali.lpoly <- try(
    lprobust(y = y_cali, x = rhat_cali, p = poly_deg, eval = eval.pts, kernel = "gau"),
    silent = TRUE
  )
  # if failed, retry with chosen bandwidth
  if (inherits(pi.cali.lpoly, "try-error")) {
    message("Automatic lprobust failed for pi.cali.lpoly; retrying with cpi bandwidth.")
    pi.cali.lpoly <- lprobust(y = y_cali, x = rhat_cali, p = poly_deg,
                              eval = eval.pts, kernel = "gau", h = mean(cpi.lpoly$Estimate[,"h"]))
  }
  
  ### IF-based bias correction oracle (knowledge of mu) ### Calibration+debias
  bc.fourier.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df.seq.fourier , basis = "fourier", g_type = g_type)
  bc.poly.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df.seq.poly, basis ="polynomial", g_type = g_type)
  bc.bs.or.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = mu, eval = eval.pts,
                  df.seq = df.seq.bs, basis = "bspline", g_type = g_type)
  bc.lpoly.or.cali <- bc.estimator(y = y, rhat = rhat_cali, a = a, muhat = mu,
                         eval = eval.pts, hseq = h.seq, deg = poly_deg, g_type = g_type)

  ### IF-based bias correction ### Calibration+debias
  bc.fourier.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.fourier , basis = "fourier", g_type = g_type)
  bc.poly.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.poly, basis ="polynomial", g_type = g_type)
  bc.bs.cali = bc.spline(y = y, rhat = rhat_cali, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.bs, basis = "bspline", g_type = g_type)

  bc.lpoly.cali <- bc.estimator(y = y, rhat = rhat_cali, a = a, muhat = muhat,
                                   eval = eval.pts, hseq = h.seq, deg = poly_deg, g_type = g_type)


  ## BC without cailbration
  bc.fourier = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.fourier , basis = "fourier", g_type = g_type)
  bc.poly = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.poly, basis ="polynomial", g_type = g_type)
  bc.bs = bc.spline(y = y, rhat = rhat, a = a, muhat = muhat, eval = eval.pts,
                  df.seq = df.seq.bs, basis = "bspline", g_type = g_type)
  bc.lpoly <- bc.estimator(y = y, rhat = rhat, a = a, muhat = muhat,
                                eval = eval.pts, hseq = h.seq, deg = poly_deg, g_type = g_type)

  ## BC.or without cailbration (with knowledga of mu)
  bc.fourier.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df.seq.fourier , basis = "fourier", g_type = g_type)
  bc.poly.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df.seq.poly, basis ="polynomial", g_type = g_type)
  bc.bs.or = bc.spline(y = y, rhat = rhat, a = a, muhat = mu, eval = eval.pts,
                      df.seq = df.seq.bs, basis = "bspline", g_type = g_type)
  bc.lpoly.or <- bc.estimator(y = y, rhat = rhat, a = a, muhat = mu,
                              eval = eval.pts, hseq = h.seq, deg = poly_deg, g_type = g_type)


  # ============================================================
  # Semi-oracle: refit with TRUE r, no correction, at each estimator's
  # realized df/h. We compute one block per estimator family.
  # ============================================================
  build_block <- function(fourier_obj, bs_obj, poly_obj, lpoly_obj) {
    list(
      fourier = semi_or_sieve(y, r, a, eval.pts,
                              get_dfh(fourier_obj, "fourier"), "fourier"),
      bs      = semi_or_sieve(y, r, a, eval.pts,
                              get_dfh(bs_obj, "bs"),      "bspline"),
      poly    = semi_or_sieve(y, r, a, eval.pts,
                              get_dfh(poly_obj, "poly"),  "polynomial"),
      lpoly   = semi_or_lpoly(y, r, eval.pts,
                              get_dfh(lpoly_obj, "lpoly"), poly_deg)
    )
  }

  semi_or_PI      <- build_block(pi.fourier,        pi.bs,        pi.poly,        pi.lpoly)
  semi_or_CPI     <- build_block(cpi.fourier,       cpi.bs,       cpi.poly,       cpi.lpoly)
  semi_or_CPI_or  <- build_block(cpi.or.fourier,    cpi.or.bs,    cpi.or.poly,    cpi.or.lpoly)
  semi_or_BC      <- build_block(bc.fourier,        bc.bs,        bc.poly,        bc.lpoly)
  semi_or_BC_or   <- build_block(bc.fourier.or,     bc.bs.or,     bc.poly.or,     bc.lpoly.or)
  semi_or_PI_cali <- build_block(pi.cali.fourier,   pi.cali.bs,   pi.cali.poly,   pi.cali.lpoly)
  semi_or_BC_cali <- build_block(bc.fourier.cali,   bc.bs.cali,   bc.poly.cali,   bc.lpoly.cali)
  semi_or_BC_or_cali <- build_block(bc.fourier.or.cali, bc.bs.or.cali,
                                    bc.poly.or.cali, bc.lpoly.or.cali)

  rslt <- list(
    OR = list(
      poly = or.poly,
      bs = or.bs,
      fourier = or.fourier,
      lpoly = list(est = or.lpoly$Estimate[,"tau.us"], se = or.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(or.lpoly$Estimate[,"h"]))

    ),
    
    PI = list(
      fourier = pi.fourier,
      poly = pi.poly,
      bs = pi.bs,
      lpoly  = list(est = pi.lpoly$Estimate[,"tau.us"], se  = pi.lpoly$Estimate[,"se.us"],
                    df.or.h = mean(pi.lpoly$Estimate[,"h"]))
    ),
    
    CPI_or = list(
      fourier =cpi.or.fourier,
      poly = cpi.or.poly,
      bs = cpi.or.bs,
      lpoly = list(est = cpi.or.lpoly$Estimate[,"tau.us"], se = cpi.or.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(cpi.or.lpoly$Estimate[,"h"]))

    ),
    
    CPI = list(
      fourier = cpi.fourier,
      poly = cpi.poly,
      bs = cpi.bs,
      lpoly = list(est = cpi.lpoly$Estimate[,"tau.us"], se  = cpi.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(cpi.lpoly$Estimate[,"h"]))

    ),
    
    PI_cali = list(
      fourier = pi.cali.fourier,
      bs = pi.cali.bs,
      poly = pi.cali.poly,
      lpoly = list(est = pi.cali.lpoly$Estimate[,"tau.us"], se  = pi.cali.lpoly$Estimate[,"se.us"],
                   df.or.h = mean(pi.cali.lpoly$Estimate[,"h"]))

    ),
    
    BC_or_cali = list(
      poly = bc.poly.or.cali,
      fourier = bc.fourier.or.cali,
      bs = bc.bs.or.cali,
      lpoly = bc.lpoly.or.cali
    ),
    
    BC_cali = list(
      poly = bc.poly.cali,
      fourier = bc.fourier.cali,
      bs = bc.bs.cali,
      lpoly = bc.lpoly.cali   
    ),
    
    BC = list(
      poly = bc.poly,
      fourier = bc.fourier,
      bs = bc.bs,
      lpoly = bc.lpoly
    ),
    
    BC_or = list(
      poly = bc.poly.or,
      fourier = bc.fourier.or,
      bs = bc.bs.or,
      lpoly = bc.lpoly.or
    ),

    semi_or = list(
      PI         = semi_or_PI,
      CPI        = semi_or_CPI,
      CPI_or     = semi_or_CPI_or,
      BC         = semi_or_BC,
      BC_or      = semi_or_BC_or,
      PI_cali    = semi_or_PI_cali,
      BC_cali    = semi_or_BC_cali,
      BC_or_cali = semi_or_BC_or_cali
    )
  )


  return(rslt)
}



suppressPackageStartupMessages({
  library(parallel); library(MASS); library(mvtnorm)
  library(splines); library(splines2)
})

set.seed(15)
d <- 200
Sigma <- 0.7^abs(row(diag(d)) - col(diag(d)))
poly_deg <- 1
beta_g <- rnorm(d, 0.5/(1:d), 0.5/(1:d))
get_g <- function(g_type, beta_g = NULL){
  if (g_type == "null")     return(function(beta_g, x) 0)
  if (g_type == "linear")   return(function(beta_g, x) as.vector(x %*% beta_g))
  if (g_type == "interact") return(function(beta_g, x) beta_g[1] * x[, 1] * x[, 2])
}
beta_r <- 0.5/(1:d)
r.x    <- function(beta_r, x) as.vector(expit(x %*% beta_r))
k1 <- 10; k2 <- 2

make_f_setting <- function(setting, ...) {
  dots <- list(...)
  k1 <- dots$k1; k2 <- dots$k2
  beta_r <- dots$beta_r; beta_g <- dots$beta_g; Sigma <- dots$Sigma
  g_type <- dots$g_type; r.x <- dots$r.x
  g.x <- get_g(dots$g_type, beta_g = dots$beta_g)
  switch(as.character(setting),
    "1" = list(
      m.t = function(t) 0.15*sin(k1*pi*t) + cos(k2*pi*t) +
        truth_fun(t, g_type, beta_r, beta_g, Sigma),
      der.m.t = function(t) 0.15*k1*pi*cos(k1*pi*t) - k2*pi*sin(k2*pi*t) +
        deriv_fun(t, g_type, beta_r, beta_g, Sigma),
      der2.m.t = function(t) -0.15*(k1*pi)^2*sin(k1*pi*t) - (k2*pi)^2*cos(k2*pi*t) +
        deriv2_fun(t, g_type, beta_r, beta_g, Sigma),
      mu.x = function(x) 0.15*sin(k1*pi*r.x(beta_r, x)) +
        cos(k2*pi*r.x(beta_r, x)) + g.x(beta_g, x)
    ),
    "2" = list(
      m.t = function(t) 4 + 3*t + 0.15*sin(4*pi*t) +
        truth_fun(t, g_type, beta_r, beta_g, Sigma),
      der.m.t = function(t) 3 + 0.6*pi*cos(4*pi*t) +
        deriv_fun(t, g_type, beta_r, beta_g, Sigma),
      der2.m.t = function(t) -0.15*(4*pi)^2*sin(4*pi*t) +
        deriv2_fun(t, g_type, beta_r, beta_g, Sigma),
      mu.x = function(x) 4 + 3*r.x(beta_r, x) +
        0.15*sin(4*pi*r.x(beta_r, x)) + g.x(beta_g, x)
    ),
    "3" = list(
      m.t = function(t) 1 + 0.1*sin(2*pi*t) +
        sin(2*pi*t)*exp(-40*(t-0.3)^2) +
        truth_fun(t, g_type, beta_r, beta_g, Sigma),
      mu.x = function(x) 1 + 0.1*sin(2*pi*r.x(beta_r, x)) +
        sin(2*pi*r.x(beta_r, x))*exp(-40*(r.x(beta_r, x)-0.3)^2) +
        g.x(beta_g, x),
      der.m.t = function(t) 0.2*pi*cos(2*pi*t) +
        2*pi*cos(2*pi*t)*exp(-40*(t-0.3)^2) +
        sin(2*pi*t)*(-80)*(t-0.3)*exp(-40*(t-0.3)^2) +
        deriv_fun(t, g_type, beta_r, beta_g, Sigma),
      der2.m.t = function(t) {
        g <- exp(-40*(t-0.3)^2)
        -0.4*(pi)^2*sin(2*pi*t) +
          (-4*pi^2*sin(2*pi*t) - 320*pi*(t-0.3)*cos(2*pi*t) +
             (-80 + 6400*(t-0.3)^2)*sin(2*pi*t)) * g +
          deriv2_fun(t, g_type, beta_r, beta_g, Sigma)
      }
    )
  )
}

run_one_config <- function(setting, g_type, n, n.tr, nsim, d, rate,
                           beta_r, beta_g, Sigma, r.x,
                           eval.pts, k1, k2, cores, outdir) {
  config <- list(
    n = n, n.tr = n.tr, rate = rate, est_r = TRUE,
    eval.pts = eval.pts, poly_deg = 1,
    d = d, beta_r = beta_r, beta_g = beta_g, Sigma = Sigma,
    r.x = r.x, g.x = NULL, mu.x = NULL,
    m.t = NULL, der.m.t = NULL, g_type = NULL,
    truth_fun = truth_fun,
    df.seq.fourier = 1:10, df.seq.poly = 1:10, df.seq.bs = 1:10,
    h.seq = seq(0.05, 0.35, length.out = 15),
    fix.df.or.h = FALSE
  )
  config$g_type <- g_type
  config$g.x <- get_g(beta_g = beta_g, g_type)

  fs <- make_f_setting(setting = setting,
                       k1 = k1, k2 = k2, g_type = g_type,
                       r.x = r.x, beta_r = beta_r, beta_g = beta_g,
                       Sigma = Sigma)
  config$m.t <- fs$m.t; config$mu.x <- fs$mu.x
  config$der.m.t <- fs$der.m.t; config$der2.m.t <- fs$der2.m.t

  message("Running [semi-oracle n.tr grid]: setting=", setting,
          " g_type=", g_type, " n=", n, " n.tr=", n.tr,
          " nsim=", nsim, " cores=", cores)

  res_list <- mclapply(
    X = seq_len(nsim),
    FUN = function(s) simu_rslt(seed = s, config = config),
    mc.cores = cores, mc.preschedule = FALSE
  )

  out <- list(result = res_list, config = config)
  stamp <- format(Sys.time(), "%m%d")
  outfile <- file.path(outdir,
              sprintf("sim_semior_setting%d_%s_n%d_ntr%d_nsim%d_%s.rds",
                      setting, g_type, n, n.tr, nsim, stamp))
  saveRDS(out, outfile)
  message("Saved: ", outfile)
  invisible(outfile)
}

if (mode %in% c("run", "both")) {
  run_one_config(
    setting = setting, g_type = g_type, n = n, n.tr = n.tr,
    nsim = nsim, d = d, beta_r = beta_r, beta_g = beta_g,
    Sigma = Sigma, r.x = r.x, eval.pts = eval.pts,
    k1 = k1, k2 = k2, cores = cores, outdir = outdir, rate = rate
  )
}

if (mode %in% c("plot", "both")) {
  generate_plots(result_dir = outdir, out_dir = file.path(outdir, "plots"))
}
