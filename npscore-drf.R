#!/usr/local/bin/Rscript

require(mvtnorm)
require(CVXR)
require(hal9001)
require(MASS)

###########################

logit <- function(x) {
  out <- log(x/(1-x))
  return(out)
}

expit <- function(x) {
  out <- exp(x)/(1 + exp(x))
  return(out)
}

qcqp <- function(tmle, H, pen, gamma.hat, tol = 1e-04) {
  
  n <- nrow(H)
  V <- t(H) %*% H/n # assuming here that H is mean-centered
  g <- gamma.hat
  
  lam.old <- 1/g
  lam.converge <- FALSE
  step <- 1e-06
  max.iter <- 1000
  iter <- 0

  while(!lam.converge) {
    iter <- iter + 1
    a <- solve(V + lam.old * pen) %*% tmle
    a.1 <- solve(V + (lam.old + step) * pen) %*% tmle
    a.2 <- solve(V + (lam.old - step) * pen) %*% tmle

    deriv.approx <-  (((t(a.1) %*% pen %*% a.1)/(t(a.1) %*% V %*% a.1)) -
                      ((t(a.2) %*% pen %*% a.2)/(t(a.2) %*% V %*% a.2)))/
                     (2 * step)
    lam <- lam.old - c((t(a) %*% pen %*% a)/(t(a) %*% V %*% a) - g)/c(deriv.approx)
    if(lam <= 0) {lam <- lam.old * .5}

    lam.converge <- abs(log(lam) - log(lam.old)) < tol
    lam.old <- lam
    if(iter > 1000) break
  }
  
  a <- (solve(V + lam.old * pen) %*% tmle)
  test.stat <-  c(abs(sum(tmle * a))/sqrt(t(a) %*% V %*% a))
  out <- list(test.stat = test.stat,
              lam = lam,
              coef = a/c(sqrt(t(a) %*% V %*% a)))
  return(out)
}

#' Sobolev basis construction
#'
#' Evaluation of Sobolev basis functions at a fixed value.
#' 
#' @param a0 Evaluation point, a scalar value.
#' @param a Predictor of interest, an n-dimensional vector.
#' @param d Dimension of basis expension for mean regression function, an integer.
#'                  
#' @return A vector containing the evaluation of each basis function at a0. 
#'         
#' @export
SobBasis <- function(a0, a, d = 20) {
  
  n <- length(a)
  a.scale <- a
  a0.scale <- a0
  
  a0.scale <- (a0.scale - ((n+1)/n) * min(a.scale))/
              (((n+1)/n) * (max(a.scale) - min(a.scale)))
  a.scale <- (a.scale - ((n+1)/n) * min(a.scale))/
             (((n+1)/n) * (max(a.scale) - min(a.scale)))
  
  phi <- numeric(d+2)
  phi[1] <- 1
  phi[2] <- a0.scale
  phi[2 + seq(1, d, 2)] <- 2^(1/2) * cos(2 * (1:(d/2)) * pi * (a0.scale))
  phi[3 + seq(1, d, 2)] <- 2^(1/2) * sin(2 * (1:(d/2)) * pi * (a0.scale))
  
  return(phi)
}

#' Outcome regression
#'
#' Estimate the conditional mean of the outcome given the exposure and covariates
#' using the highly adaptive lasso (HAL)
#' 
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#'
#'                  
#' @return A vector containing the evaluation of each basis function at x0. 
#'         
#' @export
OutcomeRegression <- function(y, a, w, continuous = TRUE) {
  n <- length(y)
  
  if(continuous) {
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  # Estimate conditional mean of y given w and a
  hal.fit <- fit_hal(X = cbind(a, w), Y = y, yolo = FALSE,
                     family = fam)
  outcome.regression <- predict(hal.fit, new_data= cbind(a,w))
  
  # Initial estimate of conditional mean at each observed a and w
  # Q as an n x n matrix with Q[i,j] as cond'l mean at (a[i], w[j])
  Q <- matrix(NA, n, n)
  for(i in 1:n) {
    ai <- a[i]
    Q[i,] <- predict(hal.fit,
                     new_data = cbind(rep(ai, n), w))
  }
  
  out <- Q
  return(out)
}

#' Propensity score estimation
#'
#' Estimate the generalized propensity score using kernel smoothing
#' 
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param no.folds Number of folds to be used for selection of bandwidth via cross-validation
#' @param no.eval Number of evaluation points (levels of predictor of interest) at which the conditional density will be estimated.
#'                Estimates at intermediate points are estimated using linear interpolation.
#' @param bw Value for bandwidth, if selected a priori. Set to NULL by default.
#'
#'  @return A list containing the following: 
#'         \describe{
#'            \item{cond.dens} Estimate of conditional density of predictor of interest, given covariates.
#'            \item{marg.dens} Estimate of marginal density of predictor of interest.
#'            \item{prop.score} Ratio of marginal density to conditional density.
#'          }             
#' @return An n-dimensional vector containing the conditional density, marginal density, and their ratior.
#'         
#' @export
EstPropScore <- function(a, w, no.folds = 10, no.eval = 20, bw = NULL) {
  n <- length(a)
  a.eval <- seq(min(a), max(a), length.out = no.eval)
  
  # use bandwidth that is optimal for kernel estimation of margianl density
  # is sensible choice if con'l density is about as smooth in a as marginal density
  # could try other choices via cross-validation, but this is comp'l cheaper
  # con'l density estimation is generally hard, and needs further study
  
  if(is.null(bw)) {
    # select bandwidth via cross-validation
    folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
    bw.seq <- exp(seq(log(.01 * sd(a)), log(25 * sd(a)), length.out = 100))
  
    pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
    for(i in 1:no.folds) {
      for(j in 1:100) {
        a.train <- a[folds != i]
        a.test <- a[folds == i]
        
        n.test <- length(a.test)
        fit.test <- numeric(n.test)
        
        for(k in 1:n.test) {
          kern <- dnorm(a[folds != i], mean = a.test[k], sd = bw.seq[j])
          fit.test[k] <- mean(kern)
        }
        
        pred.error[j,i] <- mean(-log(fit.test))
      }
    }
    avg.pred.error <- apply(pred.error, 1, mean)
    se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
    bw.opt <- max(bw.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
    # bw.opt <- bw.seq[which.min(apply(pred.error, 1, mean))]
   } else {
   # set bandwidth to a pre-specified value, if one is provided
    bw.opt <- bw
   }
  
  marg.dens <- numeric(n)
  for(i in 1:n) {
    marg.dens[i] <- mean(dnorm(a, mean = a[i], sd = 2 * bw.opt))
  }
  
  f.eval <- matrix(NA, nrow = n, ncol = no.eval)
  for(j in 1:no.eval) {
    kern.aj <- dnorm(a, mean = a.eval[j], sd = bw.opt)
    # hal.fit <- fit_hal(X = w, Y = kern.aj)
    hal.fit <- fit_hal(X = w, Y = kern.aj,
                       family = "poisson", return_x_basis = TRUE)
    f.eval[,j] <- predict(hal.fit, new_data = w)
    # f.eval[,j] <- exp(hal.fit$X.basis %*% hal.fit$coefs)
  }
  
  # set cond.dens[i,j] = cond.dens(a[i], w[j])
  cond.dens <- matrix(NA, n, n)
  for(i in 1:n) {
    k <- min(which(a.eval - a[i] >= 0))
    for(j in 1:n) {
      if(k == 1) {
        cond.dens[i,j] <- f.eval[j,1]
      } else if(k == no.eval) {
        cond.dens[i,j] <- f.eval[j,no.eval]
      } else {
        t <- (a.eval[k] - a[i])/(a.eval[k] - a.eval[k-1])
        cond.dens[i,j] <- f.eval[j,k] * (1-t) + f.eval[j,k-1] * t
      }
    }
  }
  
  marg.dens <- apply(cond.dens, 1, mean)
  prop.score <- diag(marg.dens) %*% (cond.dens^(-1))
  
  out <- list(cond.dens = cond.dens, marg.dens = marg.dens,
              prop.score = prop.score)
  return(out)
}


#' Test flat null using targeted minimum loss-based estimation.
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#' @param Q.1 Estimate of conditional mean of outcome given covariates and predictor of interest
#' @param prop.score Ratio of marginal density of exposure to conditional density, given covariates
#' @param d Dimension of basis expansion for dose response curve.
#' @param no.bs Number of multiplier bootstrap samples to be drawn.
#' @param no.folds Number of folds in cross-validation to determine smoothness of dose response curve.
#' @param gamma.hat Smoothness parameter for dose response curve, if determined a priori.
#'    
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{tmle} Estimate of expectation of drf times any given basis function.
#'            \item{eif} Estimate of efficient influence function for expectation of drf times any given basis function.
#'            \item{test.stat} Test statistic for flat null.
#'            \item{bs.samples} Bootstrap samples from null limiting distribution of test statistic.
#'            \item{p.val} p-value for test of flat null hypothesis.
#'          }
#'
#' @export
FlatNullTMLE <- function(y, a, w, continuous = TRUE,
                         Q.1, prop.score,
                         d = 20, no.bs = 1000,
                         no.folds = 10L, gamma.hat = NULL) {
  
  if(continuous) {
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  n <- length(a)
  if(d >= n) d <- round(sqrt(n))
  
  
  # Estimate of DRF at observed exposure levels
  drf.est <- apply(Q.1, 1, mean)
  
  # Construct set of basis functions
  a.scale <- a
  a.scale <- (a.scale - ((n+1)/n) * min(a.scale))/
             (((n+1)/n) * (max(a.scale) - min(a.scale)))
  
  H <- matrix(NA, nrow = n, ncol = d + 1)
  H[,1] <- a.scale
  eig <- numeric(d) # eigenvalues
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (a.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (a.scale))
    eig[2*i-1] <- (pi * i)^(-4)
    eig[2*i] <- (pi * i)^(-4)
  }
  H <- (diag(n) - matrix(1/n, n, n)) %*% H # center columns of H
  pen <- diag(c(0, 1/eig)) # penalty matrix
  
  ### Select smoothness tuning parameter
  ### Select lambda that minimizes CV error
   
  pseudo.outcome <- drf.est + diag(prop.score) * (y - diag(Q.1))
  pseudo.outcome <- pseudo.outcome - mean(pseudo.outcome)
   
  lam.max <- 2 * max(t(H[,-1] %*% diag(1/sqrt(eig))) %*% pseudo.outcome/n)
  lam.min <- lam.max * 1e-08 #lam.max * .00001
  lam.seq <- exp(seq(log(lam.min), log(lam.max), length.out = 100))
   
  no.folds <- 10
  folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
  pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
  for(i in 1:no.folds) {
    for(j in 1:100) {
      H.train <- H[folds != i,]
      H.test <- H[folds == i,]
      pseudo.train <- pseudo.outcome[folds != i]
      pseudo.test <- pseudo.outcome[folds == i]
        
      coef.train <- solve(t(H.train) %*% H.train + lam.seq[j] * pen) %*% t(H.train) %*% pseudo.train
      pred.error[j,i] <- mean((pseudo.test - H.test %*% coef.train)^2)
    }
  }
  avg.pred.error <- apply(pred.error, 1, mean)
  se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
  lam <- max(lam.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
  # lam <- lam.seq[which.min(avg.pred.error)]
   
  coef <- solve(t(H) %*% H + lam * pen) %*% t(H) %*% pseudo.outcome
  fitted <- H %*% coef
  hn <- c((fitted - mean(fitted))/sd(fitted))
  if(is.null(gamma.hat)) gamma.hat <- c(t(coef) %*% pen %*% coef)/c(var(fitted))
  # gamma.hat <- gamma.hat/c(var(fitted))
   
  # Construct clever covariate matrix  
  X <- diag(diag(prop.score)) %*% H
  
  eps <- sd((y - diag(Q.1)) * diag(prop.score))/500
  max.iter <- 5000
  iter <- 1
  
  Q.2 <- Q.1
  score.qcqp <- qcqp(tmle = c(t(y - diag(Q.2)) %*% X/n), H = H, pen = pen,
                     gamma.hat = gamma.hat)
  score.threshold <- sd((y - diag(Q.1)) * diag(prop.score) * hn)/sqrt(n * log(n))
  while(abs(score.qcqp$test.stat) >= score.threshold) {
     Q.2 <- Q.2 + eps * diag(c(H %*% score.qcqp$coef)) %*% prop.score
     score.qcqp <- qcqp(tmle = c(t(y - diag(Q.2)) %*% X/n), H = H, pen = pen,
                        gamma.hat = gamma.hat)
     iter <- iter + 1
     if(iter > max.iter) break
  }
  
  # Q.2s <- array(NA, dim = c(501, n, n))
  # Q.2s[1,,] <- Q.1
  # 
  # max.score <- numeric(500)
  # 
  # for(b in 1:500) {
  #   score.qcqp <- qcqp(tmle = c(t(y - diag(Q.2s[b,,])) %*% X/n), H = H, pen = pen,
  #                      gamma.hat = gamma.hat)
  #   max.score[b] <- score.qcqp$test.stat
  #   coef.score <- score.qcqp$coef
  #   Q.2s[b+1,,] <- Q.2s[b,,] + eps * diag(c(H %*% coef.score)) %*% prop.score
  # }
  # 
  # Q.2 <- Q.2s[min(which((max.score < sd((y - diag(Q.1)) * diag(prop.score))/sqrt(n * log(n))))),,]
  
   # Compute the TMLE
   drf.2 <- apply(Q.2, 1, mean)
   tmle <- t(H) %*% (drf.2 - mean(drf.2))/n
  
   # Estimate the influence function
   eif <- diag(drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))) %*% H +
          t(Q.1) %*% H/n - 2 * matrix(tmle, nrow = n, ncol = d + 1, byrow = TRUE)

   ### Compute the test statistic
   test.stat <- qcqp(tmle, H, pen, gamma.hat)$test.stat
   
   ### Approximate the limiting dist'n using bootstrap
   bs.samples <- rep(NA, no.bs)
   for(b in 1:no.bs) {
     mult <- rnorm(n)
     # mult <- mult - mean(mult)
     bs.tmle <- t(eif) %*% mult/n
     try(bs.samples[b] <- qcqp(bs.tmle, H, pen, gamma.hat)$test.stat)
   }
   
   ### Calculate p-value
   p.val <- mean(bs.samples > test.stat, na.rm = TRUE)
  

  out <- list(drf.est = drf.est, tmle = tmle, eif = eif,
              test.stat = test.stat, bs.samples = bs.samples,
              p.val = p.val)
  return(out)
}

#' Test flat null using one-step estimator
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#' @param Q.1 Estimate of conditional mean of outcome given covariates and predictor of interest
#' @param prop.score Ratio of marginal density of exposure to conditional density, given covariates
#' @param d Dimension of basis expansion for dose response curve.
#' @param no.bs Number of multiplier bootstrap samples to be drawn.
#' @param no.folds Number of folds in cross-validation to determine smoothness of dose response curve.
#' @param gamma.hat Smoothness parameter for dose response curve, if determined a priori.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{tmle} Estimate of expectation of drf times any given basis function.
#'            \item{eif} Estimate of efficient influence function for expectation of drf times any given basis function.
#'            \item{test.stat} Test statistic for flat null.
#'            \item{bs.samples} Bootstrap samples from null limiting distribution of test statistic.
#'            \item{p.val} p-value for test of flat null hypothesis.
#'          }
#'
#' @export
FlatNullOneStep <- function(y, a, w, continuous = TRUE,
                            Q.1, prop.score,
                            d = 20, no.bs = 1000,
                            no.folds = 10L, gamma.hat = NULL) {
  
  if(continuous) {
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  n <- length(a)
  if(d >= n) d <- round(sqrt(n))
  
  # Estimate of DRF at observed exposure levels
  drf.est <- apply(Q.1, 1, mean)
  
  # Construct set of basis functions
  a.scale <- a
  a.scale <- (a.scale - ((n+1)/n) * min(a.scale))/
             (((n+1)/n) * (max(a.scale) - min(a.scale)))
  
  H <- matrix(NA, nrow = n, ncol = d + 1)
  H[,1] <- a.scale
  eig <- numeric(d) # eigenvalues
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (a.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (a.scale))
    eig[2*i-1] <- (pi * i)^(-4)
    eig[2*i] <- (pi * i)^(-4)
  }
  H <- (diag(n) - matrix(1/n, n, n)) %*% H # center columns of H
  pen <- diag(c(0, 1/eig)) # penalty matrix
  
  ### Select smoothness tuning parameter
   pseudo.outcome <- drf.est + diag(prop.score) * (y - diag(Q.1))
   pseudo.outcome <- pseudo.outcome - mean(pseudo.outcome)
   
   # select lambda that minimizes CV error
   lam.max <- 2 * max(t(H[,-1] %*% diag(1/sqrt(eig))) %*% pseudo.outcome/n)
   lam.min <- lam.max * 1e-08 #lam.max * .00001
   lam.seq <- exp(seq(log(lam.min), log(lam.max), length.out = 100))
   
   no.folds <- 10
   folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
   pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
   for(i in 1:no.folds) {
     for(j in 1:100) {
       H.train <- H[folds != i,]
       H.test <- H[folds == i,]
       pseudo.train <- pseudo.outcome[folds != i]
       pseudo.test <- pseudo.outcome[folds == i]
        
       coef.train <- solve(t(H.train) %*% H.train + lam.seq[j] * pen) %*% t(H.train) %*% pseudo.train
       pred.error[j,i] <- mean((pseudo.test - H.test %*% coef.train)^2)
     }
   }
   avg.pred.error <- apply(pred.error, 1, mean)
   se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
   # lam <- lam.seq[which.min(avg.pred.error)]
   lam <- max(lam.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
   
   coef <- solve(t(H) %*% H + lam * pen) %*% t(H) %*% pseudo.outcome
   fitted <- H %*% coef
   if(is.null(gamma.hat)) gamma.hat <- c(t(coef) %*% pen %*% coef)/c(var(fitted))
   # gamma.hat <- gamma.hat/c(var(fitted))
  
   
   ### Compute one-step one-step
   plugin <- t(H) %*% (drf.est - mean(drf.est))/n
   one.step <- plugin + t(H) %*% (diag(prop.score) * (y - diag(Q.1)))/n
  
   # Estimate the influence function
   eif <- diag(drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))) %*% H +
          t(Q.1) %*% H/n - 2 * matrix(plugin, nrow = n, ncol = d + 1, byrow = TRUE)
   eif <- (diag(n) - matrix(1/n, n, n)) %*% eif

   ### Compute the test statistic
   test.stat <- qcqp(one.step, H, pen, gamma.hat)$test.stat
   
   ### Approximate the limiting dist'n using bootstrap
   bs.samples <- rep(NA, no.bs)
   for(b in 1:no.bs) {
     mult <- rnorm(n)
     mult <- mult - mean(mult)
     bs.one.step <- t(eif) %*% mult/n
     try(bs.samples[b] <- qcqp(bs.one.step, H, pen, gamma.hat)$test.stat)
   }
   
   ### Calculate p-value
   p.val <- mean(bs.samples > test.stat, na.rm = TRUE)
  

   out <- list(drf.est = drf.est, one.step = one.step, eif = eif,
               test.stat = test.stat, bs.samples = bs.samples,
               p.val = p.val)
  return(out)
}

#' Test flat null using one-step estimator
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#' @param Q.1 Estimate of conditional mean of outcome given covariates and predictor of interest
#' @param prop.score Ratio of marginal density of exposure to conditional density, given covariates
#' @param d Number of primitive functions to be estimated.
#' @param no.bs Number of multiplier bootstrap samples to be drawn.
#' @param no.folds Number of folds in cross-validation to determine smoothness of dose response curve.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{one.step} Estimate of expectation of any given primitive function.
#'            \item{eif} Estimate of efficient influence function for any given primitive function.
#'            \item{test.stat} Test statistic for flat null.
#'            \item{bs.samples} Bootstrap samples from null limiting distribution of test statistic.
#'            \item{p.val} p-value for test of flat null hypothesis.
#'          }
#'
#' @export
Westling <- function(y, a, w, continuous = TRUE,
                     Q.1, prop.score,
                     d = 20, no.bs = 1000, no_bins = 10L,
                     no.folds = 10L, no.eval = 20) {
  
  n <- length(a)
  if(d < n) d <- round(sqrt(n))
  
  
  # Estimate of DRF at observed exposure levels
  drf.est <- apply(Q.1, 1, mean)
  
  # Construct set of basis functions
  H <- matrix(NA, nrow = n, ncol = d)
  # a.cutoff <- quantile(a, (1:d)/(d+1))
  a.cutoff <- seq(min(a), max(a), length.out = d + 2)[-c(1, d+2)]
  for(i in 1:d) {
    H[,i] <- ifelse(a >= a.cutoff[i], 1, 0)
  }
  H <- (diag(n) - matrix(1/n, n, n)) %*% H # center columns of H
   
  ### Compute one-step
  plugin <- t(H) %*% (drf.est - mean(drf.est))/n
  one.step <- plugin + t(H) %*% (diag(prop.score) * (y - diag(Q.1)))/n
  
  # Estimate the influence function
  eif <- diag(drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))) %*% H +
          t(Q.1) %*% H/n - 2 * matrix(plugin, nrow = n, ncol = d, byrow = TRUE)
  eif <- (diag(n) - matrix(1/n, n, n)) %*% eif
   
  ####
  sd.eif <- apply(eif, 2, sd)
  test.stat <- max(abs(one.step))
   
  bs.samples <- rep(NA, no.bs)
  for(b in 1:no.bs) {
    mult <- rnorm(n)
    # mult <- mult - mean(mult)
    bs.eif <- diag(mult) %*% eif
    bs.one.step <- apply(bs.eif, 2, mean)
    bs.samples[b] <- max(abs(bs.one.step))
  }
   
   ### Calculate p-value
   p.val <- mean(bs.samples > test.stat, na.rm = TRUE)
  

  out <- list(drf.est = drf.est, one.step = one.step, eif = eif,
              test.stat = test.stat, bs.samples = bs.samples,
              p.val = p.val)
  return(out)
}

#' Construction of confidence band for dose response function using one-step estimation.
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param gamma.hat Smoothness parameter for dose response curve, if determined a priori.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#' @param Q.1 Estimate of conditional mean of outcome given covariates and predictor of interest
#' @param prop.score Ratio of marginal density of exposure to conditional density, given covariates
#' @param d Number of primitive functions to be estimated.
#' @param no.bs Number of multiplier bootstrap samples to be drawn.
#' @param no.folds Number of folds in cross-validation to determine smoothness of dose response curve.
#' @param alpha One minus the nominal coverage rate.
#' @param a0 Evaluation points for confidence band.
#' @param no.eval Number of evaluation points for confidence band, if a0 is not provided.
#'  
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{eif} Estimate of efficient influence function for any given primitive function.
#'            \item{fitted} Evaluation of drf estimate at evaluation points a0.
#'            \item{cb.lower} Lower bound of confidence band.
#'            \item{cb.upper} Upper bound of confidence band.
#'            \item{a0} Evaluation points for confidence band.
#'          }
#'
#' @export
ConfBandOneStep <- function(y, a, w, gamma.hat = NULL, continuous = TRUE,
                            Q.1, prop.score, no.folds = 10L,
                            d = 20, no.bs = 1000,
                            alpha = .05, a0 = NULL, no.eval = 50) {
  
  if(continuous) {
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  n <- length(a)
  if(d >= n) d <- round(sqrt(n))
  
  
  # Estimate of DRF at observed exposure levels
  drf.est <- apply(Q.1, 1, mean)
  
  # Construct set of basis functions
  a.scale <- a
  a.scale <- (a.scale - ((n+1)/n) * min(a.scale))/
             (((n+1)/n) * (max(a.scale) - min(a.scale)))
  
  H <- matrix(NA, nrow = n, ncol = d + 1)
  H[,1] <- a.scale
  eig <- numeric(d) # eigenvalues
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (a.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (a.scale))
    eig[2*i-1] <- (pi * i)^(-4)
    eig[2*i] <- (pi * i)^(-4)
  }
  H.means <- apply(H, 2, mean)
  H <- (diag(n) - matrix(1/n, n, n)) %*% H # center columns of H
  pen <- diag(c(0, 1/eig)) # penalty matrix
    
   
  ### Select smoothness tuning parameter
  pseudo.outcome <- drf.est + diag(prop.score) * (y - diag(Q.1))
  pseudo.outcome <- pseudo.outcome - mean(pseudo.outcome)
   
   # select lambda that minimizes CV error
   lam.max <- 2 * max(t(H[,-1] %*% diag(1/sqrt(eig))) %*% pseudo.outcome/n)
   lam.min <- lam.max * 1e-08 #lam.max * .00001
   lam.seq <- exp(seq(log(lam.min), log(lam.max), length.out = 100))
   
   no.folds <- 10
   folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
   pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
   for(i in 1:no.folds) {
     for(j in 1:100) {
       H.train <- H[folds != i,]
       H.test <- H[folds == i,]
       pseudo.train <- pseudo.outcome[folds != i]
       pseudo.test <- pseudo.outcome[folds == i]
        
       coef.train <- solve(t(H.train) %*% H.train + lam.seq[j] * pen) %*% t(H.train) %*% pseudo.train
       pred.error[j,i] <- mean((pseudo.test - H.test %*% coef.train)^2)
     }
   }
   avg.pred.error <- apply(pred.error, 1, mean)
   se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
   # lam <- lam.seq[which.min(avg.pred.error)]
   lam <- max(lam.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
   
   coef <- solve(t(H) %*% H + lam * pen) %*% t(H) %*% pseudo.outcome
   fitted <- H %*% coef
   if(is.null(gamma.hat)) gamma.hat <- c(t(coef) %*% pen %*% coef)/c(var(fitted))
   # gamma.hat <- gamma.hat/c(var(fitted))
  
   
   ### Compute one-step one-step
   plugin <- t(H) %*% (drf.est - mean(drf.est))/n
   one.step <- plugin + t(H) %*% (diag(prop.score) * (y - diag(Q.1)))/n
   
   # Estimate the influence function
   eif <- (diag(diag(prop.score) * (y - diag(Q.1))) + t(Q.1)/n) %*% H
   eif <- (diag(n) - matrix(1/n, n, n)) %*% eif

   ### Compute the test statistic
   lam <- c(qcqp(one.step, H, pen, 2 * gamma.hat)$lam)
   gram <- t(H) %*% H/n
   Pi <- H %*% solve(gram + lam * pen) %*% t(H)/n
   
   ### Approximate the limiting dist'n using bootstrap
   A <- solve(gram + lam * pen)
   bs.samples <- rep(NA, no.bs)
   for(b in 1:no.bs) {
     mult <- rnorm(n)
     mult <- mult - mean(mult)
     bs.one.step <- t(eif) %*% mult/n
     bs.samples[b] <- c(t(bs.one.step) %*% A %*% bs.one.step)
   }
   
   t.star <- quantile(bs.samples, 1-alpha)
   roughness <- c(gamma.hat * var(fitted)) + 10
   
   if(is.null(a0)) a0 <- seq(min(a), max(a), length.out = no.eval)
   cb.upper <- rep(NA, no.eval)
   cb.lower <- rep(NA, no.eval)
   
   t.star <- t.star * n
   pen <- pen/roughness * t.star
   roughness <- t.star
   
   pseudo.y <- drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))
   A <- t(H) %*% Pi %*% H # constraints for opt problems
   B <- 2 * t(pseudo.y) %*% Pi %*% H
   C <- t(pseudo.y) %*% Pi %*% pseudo.y
   
   for(j in 1:no.eval) {
     a0j <- a0[j]
     phi <- SobBasis(a0 = a0j, a = a, d = d)[-1] - H.means
     
     # upper band
     cb.coef <- Variable(d + 1)
     objective <- Minimize(-t(phi) %*% cb.coef)
     prob <- Problem(objective,
                     list(quad_form(cb.coef, pen) - roughness <= 0,
                          quad_form(cb.coef, A) - B %*% cb.coef + C - t.star <= 0))
     CVXR.result <- solve(prob)
     cb.upper[j] <- ifelse(CVXR.result$status == "solver_error",
                           NA, -CVXR.result$value)
    
     # lower band
     cb.coef <- Variable(d + 1)
     objective <- Minimize(t(phi) %*% cb.coef)
     prob <- Problem(objective,
                     list(quad_form(cb.coef, pen) - roughness <= 0,
                          quad_form(cb.coef, A) - B %*% cb.coef + C - t.star <= 0))
     CVXR.result <- solve(prob)
     cb.lower[j] <- ifelse(CVXR.result$status == "solver_error",
                           NA, CVXR.result$value)
   }

  out <- list(drf.est = drf.est, eif = eif, fitted = fitted,
              cb.lower = cb.lower, cb.upper = cb.upper, a0 = a0)
  return(out)
}

#' Construction of confidence band for dose response function using  targeted minimum loss-based estimation.
#'
#' @param y Outcome variable, an n-dimesional vector
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param gamma.hat Smoothness parameter for dose response curve, if determined a priori.
#' @param continuous Indicator of whether the outcome is continuous or binary, valued TRUE for continuous outcomes.
#' @param Q.1 Estimate of conditional mean of outcome given covariates and predictor of interest
#' @param prop.score Ratio of marginal density of exposure to conditional density, given covariates
#' @param d Number of primitive functions to be estimated.
#' @param no.bs Number of multiplier bootstrap samples to be drawn.
#' @param no.folds Number of folds in cross-validation to determine smoothness of dose response curve.
#' @param alpha One minus the nominal coverage rate.
#' @param a0 Evaluation points for confidence band.
#' @param no.eval Number of evaluation points for confidence band, if a0 is not provided.
#'  
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{eif} Estimate of efficient influence function for any given primitive function.
#'            \item{fitted} Evaluation of drf estimate at evaluation points a0.
#'            \item{cb.lower} Lower bound of confidence band.
#'            \item{cb.upper} Upper bound of confidence band.
#'            \item{a0} Evaluation points for confidence band.
#'          }
#'
#' @export
ConfBandTMLE <- function(y, a, a0 = NULL, w, gamma.hat = NULL, continuous = TRUE,
                         Q.1, prop.score, alpha = .05,
                         d = 20, no.bs = 1000,
                         no.folds = 10L, no.eval = 50) {
  
  if(continuous) {
    fam <- "gaussian"
  } else {
    fam <- "binomial"
  }
  
  n <- length(a)
  if(d >= n) d <- round(sqrt(n))
  
  
  # Estimate of DRF at observed exposure levels
  drf.est <- apply(Q.1, 1, mean)
  
  # Construct set of basis functions
  a.scale <- a
  a.scale <- (a.scale - ((n+1)/n) * min(a.scale))/
             (((n+1)/n) * (max(a.scale) - min(a.scale)))
  
  H <- matrix(NA, nrow = n, ncol = d + 1)
  H[,1] <- a.scale
  eig <- numeric(d) # eigenvalues
  for(i in 1:(d/2)) {
    H[,2*i] <- 2^(1/2) * cos(2 * i * pi * (a.scale))
    H[,2*i+1] <- 2^(1/2) * sin(2 * i * pi * (a.scale))
    eig[2*i-1] <- (pi * i)^(-4)
    eig[2*i] <- (pi * i)^(-4)
  }
  H.means <- apply(H, 2, mean)
  H <- (diag(n) - matrix(1/n, n, n)) %*% H # center columns of H
  pen <- diag(c(0, 1/eig)) # penalty matrix
    
   
  ### Select smoothness tuning parameter
  pseudo.outcome <- drf.est + diag(prop.score) * (y - diag(Q.1))
  pseudo.outcome <- pseudo.outcome - mean(pseudo.outcome)
   
   # select lambda that minimizes CV error
   lam.max <- 2 * max(t(H[,-1] %*% diag(1/sqrt(eig))) %*% pseudo.outcome/n)
   lam.min <- lam.max * 1e-08 #lam.max * .00001
   lam.seq <- exp(seq(log(lam.min), log(lam.max), length.out = 100))
   
   no.folds <- 10
   folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
   pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
   for(i in 1:no.folds) {
     for(j in 1:100) {
       H.train <- H[folds != i,]
       H.test <- H[folds == i,]
       pseudo.train <- pseudo.outcome[folds != i]
       pseudo.test <- pseudo.outcome[folds == i]
        
       coef.train <- solve(t(H.train) %*% H.train + lam.seq[j] * pen) %*% t(H.train) %*% pseudo.train
       pred.error[j,i] <- mean((pseudo.test - H.test %*% coef.train)^2)
     }
   }
   avg.pred.error <- apply(pred.error, 1, mean)
   se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
   # lam <- lam.seq[which.min(avg.pred.error)]
   lam <- max(lam.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
   
   coef <- solve(t(H) %*% H + lam * pen) %*% t(H) %*% pseudo.outcome
   fitted <- H %*% coef
   if(is.null(gamma.hat)) {
     gamma.hat <- c(t(coef) %*% pen %*% coef)/c(var(fitted))
     hn <- c((fitted - mean(fitted))/sd(fitted))
   } else {
     gammas <- numeric(length(lam.seq))
     for(j in 1:length(gammas)) {
       coef.j <- solve(t(H) %*% H + lam.seq[j] * pen) %*% t(H) %*% pseudo.outcome
       gammas[j] <- c(t(coef.j) %*% pen %*% coef.j)/c(var(H %*% coef.j))
     }
     oracle.lam <- lam.seq[which.min(abs(gammas - gamma.hat))]
     oracle.coef <- solve(t(H) %*% H + oracle.lam * pen) %*% t(H) %*% pseudo.outcome
     fitted.oracle <- H %*% oracle.coef
     hn <- c((fitted.oracle - mean(fitted.oracle))/sd(fitted.oracle))
   }
   # gamma.hat <- gamma.hat/c(var(fitted))
  
   
   ### Compute TMLE
   # Construct clever covariate matrix  
   X <- diag(diag(prop.score)) %*% H
   max.iter <- 5000
   iter <- 1
   
   eps <- sd((y - diag(Q.1)) * diag(prop.score))/500
   
   Q.2 <- Q.1
   score.qcqp <- qcqp(tmle = c(t(y - diag(Q.2)) %*% X/n), H = H, pen = pen,
                        gamma.hat = gamma.hat)
   score.threshold <- sd((y - diag(Q.1)) * diag(prop.score) * hn)/sqrt(n * log(n))
   while(abs(score.qcqp$test.stat) >= score.threshold) {
     Q.2 <- Q.2 + eps * diag(c(H %*% score.qcqp$coef)) %*% prop.score
     score.qcqp <- qcqp(tmle = c(t(y - diag(Q.2)) %*% X/n), H = H, pen = pen,
                        gamma.hat = gamma.hat)
     iter <- iter + 1
     if(iter > max.iter) break
   }
  
   # Compute the TMLE
   drf.2 <- apply(Q.2, 1, mean)
   tmle <- t(H) %*% (drf.2 - mean(drf.2))/n
  
   # Estimate the influence function
   eif <- diag(drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))) %*% H +
          t(Q.1) %*% H/n - 2 * matrix(tmle, nrow = n, ncol = d + 1, byrow = TRUE)

   ### Compute the test statistic
   lam <- c(qcqp(tmle, H, pen, 2 * gamma.hat)$lam)
   gram <- t(H) %*% H/n
   Pi <- solve(gram + lam * pen)
   
   ### Approximate the limiting dist'n using bootstrap
   bs.samples <- rep(NA, no.bs)
   for(b in 1:no.bs) {
     mult <- rnorm(n)
     mult <- mult - mean(mult)
     bs.tmle <- t(eif) %*% mult/n
     bs.samples[b] <- c(t(bs.tmle) %*% Pi %*% bs.tmle)
   }
   
   t.star <- quantile(bs.samples, 1-alpha)
   roughness <- c(gamma.hat * var(fitted)) + 10
   
   if(is.null(a0)) a0 <- seq(min(a), max(a), length.out = no.eval)
   cb.upper <- rep(NA, no.eval)
   cb.lower <- rep(NA, no.eval)
   
   t.star <- t.star * n
   tmle <- tmle * sqrt(n)
   pen <- pen/roughness * t.star
   roughness <- t.star
   
   pseudo.y <- drf.est - mean(drf.est) + diag(prop.score) * (y - diag(Q.1))
   A <- (1/n) * t(H) %*% H %*% Pi %*% t(H) %*% H # constraints for opt problems
   B <- (2/sqrt(n)) * t(tmle) %*% Pi %*% t(H) %*% H
   C <- t(tmle) %*% Pi %*% tmle
   
   for(j in 1:no.eval) {
     a0j <- a0[j]
     phi <- SobBasis(a0 = a0j, a = a, d = d)[-1] - H.means
     
     # upper band
     cb.coef <- Variable(d + 1)
     objective <- Minimize(-t(phi) %*% cb.coef)
     prob <- Problem(objective,
                     list(quad_form(cb.coef, pen) - roughness <= 0,
                          quad_form(cb.coef, A) - B %*% cb.coef + C - t.star <= 0))
     CVXR.result <- solve(prob)
     cb.upper[j] <- ifelse(CVXR.result$status == "solver_error",
                           NA, -CVXR.result$value)
    
     # lower band
     cb.coef <- Variable(d + 1)
     objective <- Minimize(t(phi) %*% cb.coef)
     prob <- Problem(objective,
                     list(quad_form(cb.coef, pen) - roughness <= 0,
                          quad_form(cb.coef, A) - B %*% cb.coef + C - t.star <= 0))
     CVXR.result <- solve(prob)
     cb.lower[j] <- ifelse(CVXR.result$status == "solver_error",
                           NA, CVXR.result$value)
   }

  out <- list(drf.est = drf.est, eif = eif, fitted = fitted,
              cb.lower = cb.lower, cb.upper = cb.upper, a0 = a0)
  return(out)
}

################
### Example
################



#######

# dose response curve
# theta0 <- function(a) {
#   a <- (a + .2)/(1.25)
#   out <- 0/3 * sign(a) * a^2 + 2 * sin(pi * (a) * exp(a))
#   return(out)
# }
# 
# # theta0 <- function(a) {
# #   out <- sin(pi * sign(a) * abs(a)^.5) + sign(a) * (a/2)^2 + exp(a*1.25)
# #   return(out)
# # }
# # 
# # theta0 <- function(x){-sin(3 * pi * x/2) * (1 + 18*x^2 * (sign(x) + 1))^(-1)}
# 
# # cond'l density of a given w
# g0 <- function(a, beta) {
#   out <- 2 * beta * expit(2 * beta * a)/(log(exp(2 * beta) + 1) - log(exp(-(2 * beta)) + 1))
#   out[beta == 0] <- .5
#   return(out)
# }
# 
# # generate from cond'l density
# rg0 <- function(n, beta) {
#   u <- runif(n)
#   out <- (1/(2 * beta)) * log(exp((log(exp(2 * beta) + 1) - log(exp(-(2 * beta)) + 1)) * u +
#                                   log(exp(- 2 * beta) + 1)) - 1)
#   out[beta == 0] <- runif(sum(beta == 0), -1, 1) #if beta == 0, generate from a uniform dist'n
#   return(out)
# }
# 
# ### generate data
# 
# n <- 1000 #sample size
# 
# # generate confounder
# sigma <- diag(2) + .5 * matrix(1, 2, 2) * (1 - diag(2))
# # sigma <- sigma * 1/4
# w <- rmvnorm(n, mean = c(0,0), sigma = sigma)
# 
# # generate treatment
# beta <- 2 * (expit(2 * w[,1] + 2 * w[,2]) - .5) #-  atan(2 * (w[,2] + w[,3]))/pi
# a <- rg0(n, beta)
# 
# # generate outcome
# y <- theta0(a) * (1/2 - 2 * sin(w[,1])) - 4 * (w[,1] + w[,2])  +  rnorm(n)
# 
# 
# # Some plots to have a look at a synthetic dataset
# plot(a, y)
# lines(sort(a), theta0(sort(a))/2, col = "red", lwd = 2)
# smoother <- smooth.spline(x = a, y = y)
# lines(smoother$x, smoother$y, col = "blue", lwd = 2)
