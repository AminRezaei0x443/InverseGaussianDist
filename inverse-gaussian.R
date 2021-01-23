# Title     : Inverse Gaussian Distribution
# Objective : Utility Functions over Inverse Gaussian Dist. (Sampling + Density)
# Created by: Amin Rezaei
# Created on: 1/17/2021

library(zoo)

dInvGauss <- function(x, lambda, mu) {
  sqrt(lambda / (2 * pi * (x * x * x))) * exp(-1 *
                                                lambda *
                                                (x - mu) *
                                                (x - mu) *
                                                (1 / (2 * mu * mu * x)))
}

dInvGaussGen <- function(lambda, mu) {
  return(function(x) {
    dInvGauss(x, lambda, mu)
  })
}

pInvGauss <- function(q, lambda, mu) {
  t1 <- sqrt(lambda / q)
  t2 <- lambda / mu
  t3 <- q / mu
  return(pnorm(t1 * (t3 - 1)) + exp(2 * t2) * pnorm(-1 * t1 * (t3 + 1)))
}

rInvGaussMH <- function(n, lambda, mu) {
  alpha <- lambda / mu
  beta <- lambda / (mu * mu)
  scale <- 1/beta
  X <- numeric(n)
  X[1] <- rgamma(1, shape=alpha, scale=scale) # initialize the chain
  dGauss <- function (x) dInvGauss(x, lambda, mu)
  dGamma <- function (x) dgamma(x, shape=alpha, scale=scale)
  for (i in 2:n) {
    Y <- rgamma(1, shape=alpha, scale=scale)
    rho <- (dGauss(Y) * dGamma(X[i - 1])) / (dGauss(X[i - 1]) * dGamma(Y))
    if (rho > 1) {
      rho <- 1
    }
    X[i] <- X[i - 1] + (Y - X[i - 1]) * (runif(1) < rho)
  }
  return(X)
}


findXNearZero <- function(targetFunc, threshold = 0.05) {
  seenHigher <- FALSE
  i <- 0.001
  while (TRUE) {
    lastV <- targetFunc(i)
    if (lastV > threshold && !seenHigher) {
      seenHigher <- TRUE
    }
    if (lastV < threshold && seenHigher) {
      break
    }
    i <- i + 0.04
  }
  return(i)
}

rInvGaussAR <- function(n, lambda, mu) {
  d <- findXNearZero(dInvGaussGen(lambda, mu), threshold = 0.01)
  k <- optimize(f=dInvGaussGen(lambda, mu), interval = c(0,d), maximum=T)$objective
  k <- k * 1.4
  needed <- n
  X <- c()
  while (needed > 0){
    m <- as.integer(k * needed) + 1
    y <- runif(m, 0, d)
    u <- runif(m, 0, 1)
    accepted <- y[(dInvGauss(y, lambda, mu) / k) >= u]
    X <- c(X, accepted)
    needed <- needed - length(accepted)
  }
  return(X[1:n])
}

rInvGaussT <- function(n, lambda, mu) {
  v <- rnorm(n)
  v <- v * v
  w <- v * mu
  c <- mu / (2. * lambda)
  x <- mu + c * (w - sqrt(w * (4 * lambda + w)))
  z <- runif(n)
  truth <- as.integer(z <= (mu / (x + mu)))
  R <- truth * x + (1 - truth) * (mu * mu / x)
  return(R)
}

rInvGauss <- function(n, lambda, mu, method = "transform") {
  if (method == "metropolis-hastings") {
    return(rInvGaussMH(n, lambda, mu))
  }else if (method == "accept-reject") {
    return(rInvGaussAR(n, lambda, mu))
  }else if (method == "transform") {
    return(rInvGaussT(n, lambda, mu))
  }
}

invGaussResults <- function(n) {
  par(mfcol = c(2, 3))
  muC <- c(1, 1, 1, 3, 3, 2)
  lambdaC <- c(0.2, 1, 3, 0.2, 1, 0.2)
  methods <- c("accept-reject", "metropolis-hastings", "transform")
  mColor <- c("#3DB1F5", "#9EE493", "#FFF689")
  lColor <- c("#191F24", "#191F24", "#191F24")
  mIdx <- 1
  for (m in methods) {
    for (i in seq_along(muC)) {
      mu <- muC[i]
      lambda <- lambdaC[i]
      X <- rInvGauss(n, lambda, mu, m)
      cH <- hist(X, plot = FALSE)
      dMax <- max(cH$density)
      i <- seq(0.001, max(X), .01)
      dY <- dInvGauss(i, lambda, mu)
      hist(X, prob = TRUE, col = mColor[mIdx], ylim = c(0, max(dMax, dY)), main=sprintf("X ~ IG(%.2f, %.2f)", mu, lambda))
      lines(i, dY, type = "l", lw = 2, col =lColor[mIdx])
      t <- -1
      if(length(unique(X)) == length(X)){
        t <- ks.test(X, pInvGauss, lambda=lambda, mu=mu)$p.value
      }

      breaks_cdf <- pInvGauss(cH$breaks, lambda=lambda, mu=mu)
      null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])
      a <- chisq.test(cH$counts, p=null.probs, rescale.p=TRUE, simulate.p.value=TRUE)
      if(t != -1){
        print(sprintf("K-S p-value - IG(%.2f,%.2f) %s method: %.4f", mu, lambda, m, t))
      }else{
        print(sprintf("KS p-value - IG(%.2f,%.2f) %s method: NaN (has ties)", mu, lambda, m))
      }
      print(sprintf("ChiSq p-value - IG(%.2f,%.2f) %s method: %.4f", mu, lambda, m, a$p.value))
    }
    mIdx <- mIdx + 1
    mtext(sprintf("Method: %s", m), outer = TRUE, cex = 1, line = -1)
  }

}

set.seed(9712017)
invGaussResults(1000)