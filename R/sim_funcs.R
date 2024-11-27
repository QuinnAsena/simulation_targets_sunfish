# library(mvtnorm) # for rmvnorm(

sim_sin_env_func <- function(seed = NULL, bX = 4, bsd = 0.5, Tlen = 100, Tmin = 100) {
  X <- bX * sin(2*pi*(1:(Tmin + Tlen)/Tlen)) + rnorm(n = Tmin + Tlen, sd = bsd)
  X <- t(as.matrix((X - mean(X[-(1:Tmin)]))))
}

sim_walk_env_func <- function(seed = NULL, stdev = 1, mu = 0, Tlen = 200) {
  X <- cumsum(rnorm(Tlen, mean = mu, sd = stdev))
  X <- t(X)
}


sim_prox_func <- function(seed = NULL, n = 3, size = 100, Tmax = 100,
                          Tstart = 100, X = NULL,
                          sigma = 0.2,
                          # r = 0.9,
                          C = matrix(0, n, n),
                          B = matrix(0, n, 1),
                          B0 = matrix(rep(0, n), nrow = n, ncol = 1),
                          V = diag(n),
                          print_plot = TRUE, ...) {

  require(mvtnorm)
  parameters <- c(as.list(environment()), list(...))
  set.seed(seed = seed)

  ### P = nrow X or nrow X + 1?!
  if (is.null(X)) {p = 1} else {p = nrow(X)}
  if (is.null(B)) {B = matrix(0, n, p)}

  C.fixed <- C
  C.fixed[C.fixed != 0] <- NA
  # C.fixed[1, 1] <- 0

  if (sum(is.na(C.fixed[upper.tri(C.fixed)])) > sum(is.na(C.fixed[lower.tri(C.fixed)])) ) {
    C.fixed[lower.tri(C.fixed)] <- C.fixed[upper.tri(C.fixed)]
  } else {
    if (sum(is.na(C.fixed[upper.tri(C.fixed)])) < sum(is.na(C.fixed[lower.tri(C.fixed)])) ) {
      C.fixed[upper.tri(C.fixed)] <- C.fixed[lower.tri(C.fixed)]
    }
  }

  C.start <- C.fixed
  C.start[is.na(C.start)] = 0.1
  # C.start[upper.tri(C.start)][which(C.start[upper.tri(C.start)] != 0)] <- 5
  C.start[upper.tri(C.start) & C.start != 0] <- 0.01
  C.start[lower.tri(C.start) & C.start != 0] <- 0.01

  V.fixed = V
  V.fixed[V.fixed != 0] <- NA
  V.fixed[1] <- 1
  # V.fixed[V.fixed != 0][V.fixed[V.fixed != 0] != 1] <- NA
  # diag(V.fixed) <- 1

  V.tmp <- matrix(c(1,0,0,0,1,-.01, 0, -.01,.2),3,3)

  V.start = V.fixed
  # V.start[is.na(V.start)] <- 0.1
  V.start[is.na(V.start)] <- V.tmp[is.na(V.start)]


  # Create data in the "normal" space of the Process Equation
  y <- matrix(0, nrow = n, ncol = Tmax + Tstart)

  if (!is.null(X)) {
    for (t in 2:(Tmax + Tstart)) {
      y[, t] <-  B0 + C %*% (y[, t - 1, drop = F] - B0 - B %*% X[, t-1, drop = F]) + B %*% X[, t, drop = F] + t(rmvnorm(n = 1, sigma = sigma ^2 * V))
    }
  } else {
    for (t in 2:(Tmax + Tstart)) {
      y[, t] <- B0 + C %*% (y[, t - 1, drop = F] - B0) + t(rmvnorm(n = 1, sigma = sigma ^2 * V))
    }
  }

  y <- t(y)
  y <- y[(Tstart + 1):(Tmax + Tstart), ]
  y <- exp(y)
  mu <- matrix(NA, Tmax, n)
  for (i in 1:Tmax) {
    mu[i,] <- y[i,] / sum(y[i,])
  }
  mu[is.na(mu)] <- 0
  # sample from a multinomial distribution for the Updating Equation
  Y <- matrix(NA, Tmax, n)
  for (t in 1:Tmax) {
    if (!all(mu[t, ] == 0)) {
    Y[t,] <- rmultinom(n = 1, size = size, prob = mu[t,])
    } else {
    # print(mu[t, ])
      Y[t, ] <- mu[t, ]
    }
  }

  if (!is.null(X) & Tstart > 0) {
  X <- t(X[,-(1:Tstart), drop = F])
  }


  if (print_plot == TRUE & !is.null(X)) {
    # par(mfrow = c(2, 1), mai = c(.1, .1, .1, .1))
    # matplot(X, typ = "l")
    matplot(Y, typ = "l", col = c("grey", "red", "blue"))
  }

  if (print_plot == TRUE & is.null(X)) {
    par(mfrow = c(1, 1))
    matplot(Y, typ = "l", col = c("grey", "red", "blue"))
  }

  return(list(Y = Y, X = X,
              C.fixed = C.fixed,
              C.start = C.start,
              V.fixed = V.fixed,
              V.start = V.start,
              params = parameters
  )
  )
}
