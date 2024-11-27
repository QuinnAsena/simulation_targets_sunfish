library(multinomialTS)
# This object is an initial fit of the model containing the necessary data
ss_mod_list_ll <- readRDS("./data/sunfish_loop_ll.rds")

### Pinus-Quercus-Fagus

Y <- ss_mod_list_ll[[1]]$Y
X <- ss_mod_list_ll[[1]]$X
Tsample <- ss_mod_list_ll[[1]]$Tsample

p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

B0.start <- ss_mod_list_ll[[1]]$B0
B.start <- ss_mod_list_ll[[1]]$B

sigma.start <- ss_mod_list_ll[[1]]$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start <- ss_mod_list_ll[[1]]$V

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

C.start <- ss_mod_list_ll[[1]]$C
C.start[3,2]=C.start[4,2]=C.start[2,3]=C.start[2,4] = -0.05

C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA
print(C.fixed)


ssmPQF <- multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                           C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                           V.fixed = V.fixed, V.start = V.start,
                           B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+07)

### Pinus-Fagus

C.start <- ss_mod_list_ll[[1]]$C
C.start[3,2]=C.start[2,3] = -0.05

C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA
print(C.fixed)

ssmPF <- multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                           C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                           V.fixed = V.fixed, V.start = V.start,
                           B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+07)

### Quercus-Fagus

C.start <- ss_mod_list_ll[[1]]$C
C.start[4,2]=C.start[2,4] = -0.05

C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA
print(C.fixed)

ssmQF <- multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                             C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                             V.fixed = V.fixed, V.start = V.start,
                             B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+07)


# Function for refitting the model to check change in AIC
refit_func <- function(mod, n_refit = 10) {

  Tsample <- mod$Tsample
  Y <- mod$Y
  X <- mod$X
  p <- ncol(X) + 1 # Number of independent variables plus intercept
  n <- ncol(Y)

  ss_seq_list <- vector(mode = "list", length = n_refit)
  ss_seq_list[[1]] <- mod
  idx <- 1

  while (idx < length(ss_seq_list)) {

    B0.start <- ss_seq_list[[idx]]$B0
    B.start <- ss_seq_list[[idx]]$B

    sigma.start <- ss_seq_list[[idx]]$sigma

    V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
    V.fixed[1] = 1

    V.start <- ss_seq_list[[idx]]$V

    B.fixed <- matrix(NA, ncol(X), n)
    B.fixed[,1] <- 0
    B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

    C.start <- ss_seq_list[[idx]]$C
    C.fixed <- C.start
    C.fixed[C.fixed != 0] <- NA
    print(C.start)
    print(C.fixed)

    ssm <- multinomialTS::mnTS(Y = Y, X = X, Tsample = Tsample, B0.start = B0.start, B.start = B.start,
                               C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                               V.fixed = V.fixed, V.start = V.start,
                               B.fixed = B.fixed, dispersion.fixed = 1, maxit.optim = 1e+07)
    print(ssm$AIC)
    idx <- idx + 1
    print(idx)
    ss_seq_list[[idx]] <- ssm

  }
  return(ss_seq_list)
}

ssmPQF_refit <- refit_func(ssmPQF)
ssmPF_refit <- refit_func(ssmPF)
ssmQF_refit <- refit_func(ssmQF)

# saveRDS(ssmPQF_refit, "./data/ssmPQF_refit.rds")
# saveRDS(ssmPF_refit, "./data/ssmPF_refit.rds")
# saveRDS(ssmQF_refit, "./data/ssmQF_refit.rds")

sapply(ssmPQF_refit, \(x) {x$AIC})


ssmPQF_refit_boot <- ssmPQF_refit[[9]]
summary(ssmPQF_refit_boot)

ssmPQF_refit <- readRDS("./data/ssmPQF_refit.rds")
ssmPF_refit <- readRDS("./data/ssmPF_refit.rds")
ssmQF_refit <- readRDS("./data/ssmQF_refit.rds")


library(future)
library(furrr)
# Bootstrap models in parallel
plan(multisession, workers = 3)
start_time <- Sys.time()

bootstrap_results <- future_map(list(ssmPQF_refit[[1]], ssmPF_refit[[1]], ssmQF_refit[[2]]),
                     multinomialTS::boot.mnTS, reps = 1500, .options = furrr_options(seed = TRUE))

end_time <- Sys.time()
end_time - start_time
saveRDS(bootstrap_results, "./data/bootstrap_results.rds")
bootstrap_results <- readRDS("./data/bootstrap_results.rds")


lapply(bootstrap_results, \(pars) {
  most_pars <- pars$all_mods_pars
  most_pars <- most_pars[most_pars[, colnames(most_pars)%in% c("opt.convergence")] == 0, ]
  print(dim(most_pars)[1])
  most_pars_mean <- matrixStats::colMeans2(most_pars, na.rm = T)
  most_pars_sd <- matrixStats::colSds(most_pars, na.rm = T)
  most_pars_upper_68 <- most_pars_mean + most_pars_sd
  most_pars_lower_68 <- most_pars_mean - most_pars_sd
  t_scores <- most_pars_mean/most_pars_sd
  p_vals <- 2 * pnorm(q = abs(t_scores), lower.tail = F)

  res <- cbind(boot_mean = most_pars_mean, se = most_pars_sd,
               t = t_scores, P = p_vals, upper68 = most_pars_upper_68,
               lower68 = most_pars_lower_68)
  res[rownames(res)[!grepl("^v.*", rownames(res))], ]
})



