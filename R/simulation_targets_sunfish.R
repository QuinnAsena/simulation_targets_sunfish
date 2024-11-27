# file <- c("../data/sunfish_counts_with_dates.rds")
# library(tidyverse)
# library(patchwork)
# library(MultinomialStateSpace)
# source("./R/sim_funcs.R")

get_data <- function(file) {
  sunfish_counts_with_dates <- readRDS(file)
  sunfish_counts_with_dates <- sunfish_counts_with_dates[nrow(sunfish_counts_with_dates):1, , drop = F]

  sunfish_years <- sunfish_counts_with_dates[ ,1, drop = F]
  X <- as.matrix(scale(sunfish_years))

  sunfish_time_forward <- abs(abs(sunfish_counts_with_dates[ ,1, drop = F]) - max(sunfish_counts_with_dates[ ,1, drop = F]))+1

  sunfish_tsample <- floor(sunfish_time_forward / 100)
  sunfish_tsample <- sunfish_tsample[, , drop = T] + 1 # +1 for index not to start at 0

  sunfish_spp_long <- pivot_longer(sunfish_counts_with_dates, cols = -sunfish_years, names_to = "variablename") %>%
    mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus.*"), "Pinus"),
           variablename = replace(variablename, stringr::str_detect(variablename, "Acer*"), "Acer")) %>%
    group_by(variablename, sunfish_years) %>%
    summarise(value = sum(value), .groups = 'keep')

  sunfish_PT <- sunfish_spp_long %>%
    filter(variablename %in% c("Pinus", "Tsuga"))

  sunfish_other <- sunfish_spp_long %>%
    filter(!variablename %in% c("Pinus", "Tsuga")) %>%
    mutate(variablename = "other") %>%
    group_by(variablename, sunfish_years) %>%
    summarise(value = sum(value), .groups='keep')

  sunfish_spp_wide <- bind_rows(sunfish_other, sunfish_PT) %>%
    pivot_wider(id_cols = sunfish_years, names_from = variablename, values_from = value) %>%
    arrange(desc(sunfish_years)) %>%
    as.data.frame()

  sunfish_spp <- sunfish_spp_wide[,-1]

  age_range_sim <- seq(from = 1, to = sunfish_tsample[length(sunfish_tsample)], length.out = 150)
  X_long_sim <- scale(age_range_sim)

  years_mat <- matrix(NA, nrow = sunfish_tsample[length(sunfish_tsample)], ncol = 1)
  sunfish_years <- as.matrix(sunfish_years)
  for (i in seq_along(sunfish_tsample)) {
    years_mat[sunfish_tsample[i], 1] <- sunfish_years[i]
  }
  X_long <- forecast::na.interp(years_mat)
  X_long <- abs(abs(X_long) - max(X_long))+1
  X_long <- as.numeric(scale(X_long))


  data <- list(
    Y = sunfish_spp,
    X = X,
    X_long = X_long,
    X_long_sim = X_long_sim,
    sunfish_years = sunfish_years,
    sunfish_tsample = sunfish_tsample
  )
  return(data)

}


# Fit with Y_agg and X_long cv --------------------------------------------

fit_model <- function(data) {

  source("./../cpp_testing_groud/simulate_func_TEMP.R")
  library(Rcpp)
  library(RcppArmadillo)
  library(minqa)
  sourceCpp("./../R/source_multinomialSS.cpp")

  list2env(data, globalenv())

  n <- ncol(Y)
  p <- ncol(X) + 1 # Number of independent variables plus intercept

  V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
  V.fixed[1] <- 1

  # V.fixed = diag(NA, n) # Covariance matrix of environmental variation in process eq
  # V.fixed[2, 3] <- NA
  # V.fixed[3, 2] <- NA
  # V.fixed[1] <- 1

  dispersion.fixed <- 1
  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(Y = Y, X = X, B.start = B.start, B.fixed = B.fixed,
                              V.fixed = V.fixed, dispersion.fixed = dispersion.fixed)
  summary(glmm_mod)

  B0.start <- glmm_mod$B[1, , drop = F]
  B.start <- glmm_mod$B[2, , drop = F]

  C.fixed <- diag(NA, n)
  C.fixed[2, 3] <- NA
  C.fixed[3, 2] <- NA
  C.fixed[C.fixed == F] <- 0

  C.start = C.fixed
  C.start[is.na(C.start)] = .01
  diag(C.start) <- 0.5

  sigma.start <- glmm_mod$sigma

  V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
  V.fixed[1] = 1

  V.start = V.fixed
  V.start <- glmm_mod$V
  V.start <- diag(diag(V.start))

  B.fixed <- matrix(NA, ncol(X), n)
  B.fixed[,1] <- 0
  B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  start_time <- Sys.time()

  ss_mod <- tryCatch(
    {
      multinomialTS::mnTS(Y = Y, X = X_long, B0.start = B0.start, B.start = B.start,
                    C.start = C.start, C.fixed = C.fixed, B0.fixed = B0.fixed,
                    V.fixed = V.fixed, V.start = V.start,
                    B.fixed = B.fixed, dispersion.fixed = 1,
                    Tsample = sunfish_tsample, sigma.start = sigma.start)
    },
    error=function(cond) {
      message("Here's the original error message:")
      message(conditionMessage(cond))
      return(cond)
    }
  )
  end_time <- Sys.time()

  summary(ss_mod)

  model <- list(glmm_mod =  glmm_mod,
                ss_mod = ss_mod,
                ss_time = end_time - start_time,
                Y = Y,
                X = X,
                X_long = X_long,
                X_long_sim = X_long_sim,
                sunfish_years = sunfish_years)
  return(model)
}


# simulate Y with custom inputs -------------------------------------------

Y_sim <-
  function(model,
           reps = 3,
           C = matrix(c(.5, 0, 0, 0, .9, 0, 0, -.3, .6), 3, 3),
           V = matrix(c(1, 0.7, 0.7, 0.7, .6, .5, 0.7, .5, .6), 3, 3),
           B = c(0, 0.25, -0.5),
           X_type = c("time", "pulse")) {

    list2env(model, globalenv())

    n <- ncol(Y)
    seed <- 1984
    set.seed(seed)

    if (X_type == "pulse") {
      X_pulse <- c(rep(0, 30),
                   rep(1, 30),
                   rep(0, 30),
                   rep(1, 30),
                   rep(0, 30))
      X <- t(X_pulse)
    } else {

      if (X_type == "time") {
        X <- t(X_long_sim)
      }
    }

    Tmax <- ncol(X)

    par_list <- list(
      sim_1 = list(C = diag(diag(C)),
                   B = matrix(B),
                   V = diag(diag(V))),

      sim_2 = list(C = diag(diag(C)),
                   B = matrix(B),
                   V = V),

      sim_3 = list(C = C,
                   B = matrix(B),
                   V = diag(diag(V))),

      sim_4 = list(C = C,
                   B = matrix(B),
                   V = V),

      sim_5 = list(C = C,
                   B = matrix(0, ncol = 1, nrow = 3),
                   V = V),

      sim_6 = list(C = diag(diag(C)),
                   B = matrix(0, ncol = 1, nrow = 3),
                   V = V)
      )

    sim_list <- vector(mode = "list", length = length(par_list))
    for (i in seq_along(par_list)) {

      sim_list[[i]] <- replicate(sim_prox_func(n = n, size = 100, Tmax = Tmax,
                                               Tstart = 0, X = X,
                                               sigma = 0.2,
                                               C = par_list[[i]]$C,
                                               B = par_list[[i]]$B,
                                               B0 = matrix(rep(0, n), nrow = n, ncol = 1),
                                               V = par_list[[i]]$V,
                                               print_plot = T), n = reps)
    }

    return(sim_list)
  }



###
fit_sim <- function(sim, data, sub_sample = TRUE) {

  ss_mod_fit_list <- vector(mode = "list", length = length(sim))

  for (i in seq_along(sim)) {
    ss_summary <- vector(mode = "list", length = ncol(sim[[i]]) )

      for (j in seq_along(ss_summary) ) {
      Y <- sim[[i]][1, j]$Y[11:nrow(sim[[i]][1, j]$Y), ]
      X <- t(sim[[i]][2, j]$X[ ,11:ncol(sim[[i]][2, j]$X), drop = F])
      n <- ncol(Y)
      p <- ncol(X) + 1 # Number of independent variables plus intercept

      if (sub_sample) {
        Y <- Y[data$sunfish_tsample, ]
      }

      B.fixed <- matrix(NA, ncol(X), n)
      B.fixed[,1] <- 0

      print(i)
      print(j)

      B0.start <- c(t(sim[[i]][7, j]$params$B0))
      B.start <- c(t(sim[[i]][7, j]$params$B))

      V.fixed <- matrix(NA, n, n)
      V.fixed[1] <- 1
      print(V.fixed)

      # V.start <- ssm$V
      V.start <- sim[[i]][7, j]$params$V
      print(V.start)

      ## to mirror sim
      # C.fixed  <- sim[[i]][3, j]$C.fixed

      C.fixed <- diag(NA, n)
      C.fixed[2, 3] <- NA
      C.fixed[3, 2] <- NA
      C.fixed[C.fixed == F] <- 0
      print(C.fixed)

      C.start  <- sim[[i]][7, j]$params$C
      print(C.start)
      # C.start = ssm$C

      t0 <- Sys.time()
        ss_summary[[j]] <- tryCatch(
          {
            ss_mod <- multinomialTS::mnTS(Y = Y, X = X, B0.start = B0.start, B.start = B.start,
                                    C.fixed = C.fixed, C.start = C.start,
                                    V.fixed = V.fixed, V.start = V.start,
                                    B.fixed = B.fixed, dispersion.fixed = 1,
                                    Tsample = data$sunfish_tsample)
            show(summary(ss_mod))
            t1 <- Sys.time()
            ss_mod <- list(ss_mod = ss_mod,
                           time_taken = t1 - t0)
          },
          error=function(cond) {
            message("Here's the original error message:")
            message(conditionMessage(cond))
            return(list(err = cond, Y = Y))
          }
        )
      # t1 <- Sys.time()
    }
    ss_mod_fit_list[[i]] <- ss_summary
  }
  # end_time <- Sys.time()

  return(ss_mod_fit_list)
}



fit_sim_ll <- function(sim, ssms, data, sub_sample = TRUE, B_pos = 2) {

  ss_mod_fit_list <- vector(mode = "list", length = length(sim))

  for (i in seq_along(sim)) {
    ss_summary <- vector(mode = "list", length = ncol(sim[[i]]) )

    for (j in seq_along(ss_summary) ) {
      Y <- sim[[i]][1, j]$Y[11:nrow(sim[[i]][1, j]$Y), ]
      X <- t(sim[[i]][2, j]$X[ ,11:ncol(sim[[i]][2, j]$X), drop = F])
      n <- ncol(Y)
      p <- ncol(X) + 1

      if (sub_sample) {
        Y <- Y[data$sunfish_tsample, ]
      }

      B.fixed <- matrix(NA, ncol(X), n)
      B.fixed[ ,1] <- 0
      B.fixed[ ,B_pos] <- 0
      print(B.fixed)

      print(i)
      print(j)

      B0.start <- ssms[[i]][[j]][[1]]$B0
      B.start <- ssms[[i]][[j]][[1]]$B
      print(B.start)
      B.start[ ,B_pos] <- 0
      print(B.start)

      V.fixed <- matrix(NA, n, n)
      V.fixed[1] <- 1
      print(V.fixed)

      # V.start <- ssm$V
      V.start <- ssms[[i]][[j]][[1]]$V
      print(V.start)

      C.fixed <- diag(NA, n)
      C.fixed[2, 3] <- NA
      C.fixed[3, 2] <- NA
      C.fixed[C.fixed == F] <- 0
      print(C.fixed)

      C.start  <- ssms[[i]][[j]][[1]]$C
      print(C.start)
      # C.start = ssm$C

      t0 <- Sys.time()
      ss_summary[[j]] <- tryCatch(
        {
          ss_mod <- multinomialTS::mnTS(Y = Y, X = X, B0.start = B0.start, B.start = B.start,
                                  C.fixed = C.fixed, C.start = C.start,
                                  V.fixed = V.fixed, V.start = V.start,
                                  B.fixed = B.fixed, dispersion.fixed = 1)
          show(summary(ss_mod))
          t1 <- Sys.time()
          p_ll <- pchisq(2 * (ssms[[i]][[j]][[1]]$logLik - ss_mod$logLik), df = 1, lower.tail = F)
          ss_mod <- list(ss_mod = ss_mod,
                         time_taken = t1 - t0,
                         p_ll = p_ll)
        },
        error=function(cond) {
          message("Here's the original error message:")
          message(conditionMessage(cond))
          return(list(err = cond, Y = Y))
        }
      )
    }
    ss_mod_fit_list[[i]] <- ss_summary
  }
  return(ss_mod_fit_list)
}




###
boot_sim <- function(ss_mod_fit_list, n_boot = 5, n_sim = 10, case = NULL) {

  if(!is.null(case)) {
    ss_mod_fit_list <- ss_mod_fit_list[case]
  }

  ssms <- lapply(ss_mod_fit_list, \(sce) {
    ss_without_error <- lapply(sce, \(rep) {
      ssm <- rep[[1]]
    })
    ss_mods <- ss_without_error[lapply(ss_without_error, class) == "multinomialSS"]
    ss_mods
  })

  ssms <- lapply(ssms, \(filt) {
    Filter(function(x) x$opt.convergence == 0, filt)
  })

  ss_mod_boot_list <- vector(mode = "list", length = length(ssms))

  for (i in seq_along(ssms)) {
    ss_rep <- vector(mode = "list", length = n_sim)

      counter <- 1
      success_count <- 0
      while (success_count < n_sim) {
        print(paste0("ss_case_i = ", i, "  counter = ", counter, "  success counter = ", success_count))
        ss_mod <- tryCatch(
          {
            multinomialTS::boot.mnTS(ssms[[i]][[counter]], reps = n_boot, maxit.optim = 1e+05)
          },
          error=function(cond) {
            message("Here's the original error message:")
            message(conditionMessage(cond))
            return(NULL)
          }
        )
        if (!is.null(ss_mod)) {
          success_count <- success_count + 1
          mod_coef <- coef(ssms[[i]][[counter]])
          full_mod <- ssms[[i]][[counter]]
          ss_rep[[success_count]] <- c(list(boot_mod = ss_mod), mod_coef, list(full_mod = full_mod))
          counter <- counter + 1
        } else {
          counter <- counter + 1
        }
      }
    ss_mod_boot_list[[i]] <- ss_rep
  }
  return(ss_mod_boot_list)

}
