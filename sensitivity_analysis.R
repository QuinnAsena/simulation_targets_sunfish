if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(targets, tidyverse, ggtext, patchwork, future, furrr, future.apply,
               multinomialTS, viridis, RColorBrewer)

tar_load(data_sim_time_c5)

sim_dat <- data_sim_time_c5[[4]][,1]

X <- t(sim_dat$X)
X <- X[11:nrow(X), , drop = FALSE]
Y <- sim_dat$Y[11:nrow(sim_dat$Y), ]

rm(data_sim_time_c5, sim_dat)

prop_dat <- seq(1, 0.2, -0.1)

y_names <- paste0("Y_", prop_dat)

Y_list <- lapply(prop_dat, \(i) {
  dat <- Y
  dat[-unique(round(seq(1, nrow(dat), length.out = (nrow(dat)*i)))), ] <- NA
  return(dat)
})
names(Y_list) <- y_names


p <- ncol(X) + 1
n <- ncol(Y)


plan(multisession, workers = length(Y_list))
res_par <- furrr::future_map(Y_list, \(x) {
  Y <- x
  Tsample <- which(rowSums(Y) != 0)
  V.fixed = diag(n)
  # V.fixed = matrix(NA, n, n)
  # V.fixed[1] = 1

  B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
  B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

  glmm_mod <- multinomialTS::mnGLMM(Y = Y[Tsample, ],
                                    X = X[Tsample, ,drop = F],
                                    B.start = B.start, B.fixed = B.fixed,
                                    V.fixed = V.fixed)

  B0.start <- glmm_mod$B[1, , drop = F]
  B.start <- glmm_mod$B[2:p, , drop = F]

  sigma.start <- glmm_mod$sigma

  V.fixed = matrix(NA, n, n)
  V.fixed[1] = 1

  V.start <- glmm_mod$V
  # V.start <- diag(diag(V.start))

  B.fixed <- matrix(NA, ncol(X), n)
  B.fixed[,1] <- 0
  B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

  # Set-up C with interactions
  C.start = .5 * diag(n)
  C.start[2,3] = C.start[3,2] = -0.5
  C.fixed <- C.start
  C.fixed[C.fixed != 0] <- NA

  # Model with interactions
  start_time <- Sys.time()
  mnTS_mod <- mnTS(Y = Y[Tsample, ],
                   X = X, Tsample = Tsample,
                   B0.start = B0.start, B0.fixed = B0.fixed,
                   B.start = B.start, B.fixed = B.fixed,
                   C.start = C.start, C.fixed = C.fixed,
                   V.start = V.start, V.fixed = V.fixed,
                   dispersion.fixed = 1, maxit.optim = 1e+6)

  mnTS_mod_boot <- multinomialTS::boot(mnTS_mod, reps = 1000)

  end_time <- Sys.time()
  end_time - start_time
  return(mnTS_mod_boot)
}, .options = furrr_options(seed = 1984))

# saveRDS(res_par, "./sensitivity_results/samp_reduction.rds")

# Plot temporal resolution sensitivity ------------------------------------

samp_sensitivity <- readRDS("./sensitivity_results/samp_reduction.rds")

mods_boot <- map(samp_sensitivity, ~ {
  as_tibble(.x[[2]]) |>
    pivot_longer(-c(logLik, opt.convergence))
}) |>
  bind_rows(.id = "resolution")

mods_boot_68 <- mods_boot |>
  #  filter(opt.convergence == 0) |> # option to filter based on convergence
  group_by(resolution, name) |>
  summarise(boot_mean = mean(value),
            boot_med = median(value),
            boot_sd = sd(value),
            upper_68 = quantile(value, probs = 0.84),
            lower_68 = quantile(value, probs = 0.16)) |>
  mutate(t_scores = boot_mean / boot_sd,
         p_vals = 2 * pnorm(q = abs(t_scores), lower.tail = F),
         sig = p_vals < 0.05)
#

mods_boot_table <- mods_boot_68 |>
  mutate(targets = case_when(name == "x1.y2" ~ 0.25,
                        name == "x1.y3" ~ -0.5,
                        name == "sp.y1.y1" ~ 0.5,
                        name == "sp.y2.y2" ~ 0.9,
                        name == "sp.y3.y2" ~ 0,
                        name == "sp.y2.y3" ~ -0.5,
                        name == "sp.y3.y3" ~ 0.6,
                        .default = 0))

mods_boot_68_B <- mods_boot_table |>
  filter(name %in% c("x1.y2", "x1.y3")) |>
  mutate(resolution = fct(resolution, levels = (paste0("Y_", rev(prop_dat)))))

dat_labs <- paste0("Y = ", rev(prop_dat))
names(dat_labs) <- levels(mods_boot_68_B$resolution)

X_labs <- c("x1.y2" = "B<sub>1, 2</sub>",
            "x1.y3" = "B<sub>1, 3</sub>")

mod_plots_B <- mods_boot_68_B |>
  ggplot(aes(x = name, y = boot_med)) +
  geom_point() +
  geom_point(data = mods_boot_68_B, aes(x = name, y = targets), colour = "red", shape = 5) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                width = .2, alpha = 0.5) +
  scale_x_discrete(labels = X_labs) +
  facet_wrap(~ resolution, nrow = 1,
             labeller = labeller(resolution = dat_labs)) +
  labs(x = NULL, y = "B Coefficient") +
  theme_bw() +
  theme(
    axis.text = element_markdown(size = 10, angle = 45, hjust = 1),
    axis.title = element_text(size = 10),
    strip.background = element_rect(fill = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
  )


## Bootstrap plotting C ---------------------------------------------------
C_labs2 <- c("sp.y1.y1" = "C<sub>1, 1</sub>",
             "sp.y2.y2" = "C<sub>2, 2</sub>",
             "sp.y3.y3" = "C<sub>3, 3</sub>",
             "sp.y2.y3" = "C<sub>2, 3</sub>",
             "sp.y3.y2" = "C<sub>3, 2</sub>")

mods_boot_68_C <- mods_boot_table |>
  filter(name %in% c("sp.y1.y1", "sp.y2.y2", "sp.y3.y3", "sp.y2.y3", "sp.y3.y2")) |>
  mutate(name = fct_relevel(name, "sp.y1.y1", "sp.y2.y2", "sp.y3.y3",
                                "sp.y2.y3", "sp.y3.y2"))

mod_plots_C <- mods_boot_68_C |>
  ggplot(aes(x = name, y = boot_med)) +
  geom_point() +
  geom_point(data = mods_boot_68_C, aes(x = name, y = targets), colour = "red", shape = 5) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                width = .2, alpha = 0.5) +
  labs(x = NULL, y = "C coefficient") +
  scale_x_discrete(labels = C_labs2) +
  facet_wrap(~resolution, nrow = 1, labeller = labeller(resolution = dat_labs)) +
  theme_bw() +
  theme(
    axis.text = element_markdown(size = 10, angle = 45, hjust = 1),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = NA)
  )

bootstrap_data_reduction <- mod_plots_B / mod_plots_C

ggsave(
  filename = "./figures/bootstrap_data_reduction.png",
  plot = bootstrap_data_reduction,
  device = "png",
  dpi = 300,
  width = 10,
  height = 6,
  units = "in"
)

ggsave(
  filename = "./figures/bootstrap_data_reduction.svg",
  plot = bootstrap_data_reduction,
  device = "svg",
  width = 10,
  height = 6,
  units = "in"
)

# 500 reps ----------------------------------------------------------------

prop_dat <- seq(1, 0.2, -0.1)
names(prop_dat) <- paste("prop_", prop_dat)

sim_dat_reps <- data_sim_time_c5[[4]]
# sim_dat_reps <- sim_dat_reps[,1:3]

plan(multisession, workers = length(prop_dat))
prop_res <- future_lapply(prop_dat, \(prop) {
  prop <- prop
  ss_summary <- vector(mode = "list", length = ncol(sim_dat_reps))
  for (j in 1:ncol(sim_dat_reps)) {
    Y <- sim_dat_reps[1, j]$Y[11:nrow(sim_dat_reps[1, j]$Y), ]
    X <- t(sim_dat_reps[2, j]$X[ ,11:ncol(sim_dat_reps[2, j]$X), drop = F])
    n <- ncol(Y)
    p <- ncol(X) + 1 # Number of independent variables plus intercept
    Y[-unique(round(seq(1, nrow(Y), length.out = (nrow(Y)*prop)))), ] <- NA
    print(Y[1:10, ])
    Tsample <- which(rowSums(Y) != 0)

    B.fixed <- matrix(NA, ncol(X), n)
    B.fixed[,1] <- 0

    B0.start <- c(t(sim_dat_reps[7, j]$params$B0))
    B.start <- c(t(sim_dat_reps[7, j]$params$B))

    V.fixed <- matrix(NA, n, n)
    V.fixed[1] <- 1

    V.start <- sim_dat_reps[7, j]$params$V

    ## to mirror sim_dat_reps
    # C.fixed  <- sim_dat_reps[3, j]$C.fixed

    C.fixed <- diag(NA, n)
    C.fixed[2, 3] <- NA
    C.fixed[3, 2] <- NA
    C.fixed[C.fixed == F] <- 0
    print(C.fixed)

    C.start  <- sim_dat_reps[7, j]$params$C
    print(C.start)
    # C.start = ssm$C

    ss_summary[[j]] <- tryCatch(
      {
        multinomialTS::mnTS(Y = Y, X = X, B0.start = B0.start, B.start = B.start,
                                      C.fixed = C.fixed, C.start = C.start,
                                      V.fixed = V.fixed, V.start = V.start,
                                      B.fixed = B.fixed, dispersion.fixed = 1,
                                      Tsample = Tsample, maxit.optim = 1e+6)
      },
      error=function(cond) {
        message("Here's the original error message:")
        message(conditionMessage(cond))
        return(list(err = cond, Y = Y))
      }
    )
  }
  return(ss_summary)
})

saveRDS(prop_res, "./sensitivity_results/samp_reduction500reps.rds")

prop_res <- readRDS("./sensitivity_results/samp_reduction500reps.rds")

pars <- lapply(prop_res, \(case) {
  all_pars <- lapply(case, \(rep) {
    v.pars <- names(rep$par)[grep(names(rep$par), pattern = "v.")]
    par.names.noshow.se <- c("dispersion","sigma", v.pars)
    par.show.se <- rep$par[!is.element(names(rep$par), par.names.noshow.se)]
    par.all <- c(par.show.se, logLik = rep$logLik, opt.convergence = rep$opt.convergence)
  })
  bind_rows(all_pars, .id = "rep")
})

pars_wide_bs <- bind_rows(pars, .id = "prop") %>%
  select(c(prop, rep, opt.convergence), starts_with("x"))


pars_long_bs <- pars_wide_bs %>%
  pivot_longer(cols = -c(prop, rep, opt.convergence)) %>%
  mutate(
    bs = case_when(name == "x1.y2" ~ 0.25,
                   name == "x1.y3" ~ -0.5,
                   .default = NA)
  )

pars_long_bs %>%
  group_by(prop, name, bs) %>%
  summarise(med = median(value),
            quant25 = quantile(value, probs = 0.25),
            quant75 = quantile(value, probs = 0.75)) %>%
  print(n = 24)


X_labs <- c("x1.y2" = "B<sub>1, 2</sub>",
            "x1.y3" = "B<sub>1, 3</sub>")


prop_labs <- paste0("Proportion of data = ", prop_dat)
names(prop_labs) <-  names(prop_dat)

B_procession_plot <- ggplot(data = pars_long_bs) +
  geom_violin(aes(x = name, y = value, alpha = 0.5),
              draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3,
              fill = brewer.pal(11, "BrBG")[3]) +
  geom_point(aes(x = name, y = bs), colour = viridis(10)[10], size = 1) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = X_labs) +
  facet_wrap(~prop, labeller = labeller(prop = prop_labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_markdown(size = 12),
        strip.text = element_markdown(size = 11),
        strip.text.y = element_markdown(size = 12),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1)
  )

ggsave(
  filename = "./figures/replicate_data_reductionB.png",
  plot = B_procession_plot,
  device = "png",
  dpi = 300,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  filename = "./figures/replicate_data_reductionB.svg",
  plot = B_procession_plot,
  device = "svg",
  width = 8,
  height = 8,
  units = "in"
)

# C sensitivity -----------------------------------------------------------

pars_wide_cs <- bind_rows(pars, .id = "prop") %>%
  select(c(prop, rep, opt.convergence), starts_with("sp."))


pars_long_cs <- pars_wide_cs %>%
  pivot_longer(cols = -c(prop, rep, opt.convergence)) %>%
  mutate(
    cs = case_when(name == "sp.y1.y1" ~ 0.5,
                   name == "sp.y2.y2" ~ 0.9,
                   name == "sp.y3.y2" ~ 0,
                   name == "sp.y2.y3" ~ -0.5,
                   name == "sp.y3.y3" ~ 0.6,
                   .default = NA),
    name = fct_relevel(name, "sp.y1.y1", "sp.y2.y2", "sp.y3.y3",
                       "sp.y2.y3", "sp.y3.y2")
  )

pars_long_group_c <- pars_long_cs %>%
  group_by(prop, name, cs) %>%
  summarise(med = median(value),
            quant25 = quantile(value, probs = 0.25),
            quant75 = quantile(value, probs = 0.75))



C_labs <- c("sp.y1.y1" = "sp.1, sp.1",
            "sp.y2.y2" = "sp.2, sp.2",
            "sp.y3.y3" = "sp.3, sp.3",
            "sp.y2.y3" = "sp.2, sp.3",
            "sp.y3.y2" = "sp.3, sp.2"
)

C_labs2 <- c("sp.y1.y1" = "C<sub>1, 1</sub>",
             "sp.y2.y2" = "C<sub>2, 2</sub>",
             "sp.y3.y3" = "C<sub>3, 3</sub>",
             "sp.y2.y3" = "C<sub>2, 3</sub>",
             "sp.y3.y2" = "C<sub>3, 2</sub>"
)

C_procession_plot_lin <- ggplot(data = pars_long_cs) +
  geom_violin(aes(x = name, y = value, alpha = 0.5),
              draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3,
              fill = brewer.pal(11, "BrBG")[3]) +
  geom_point(aes(x = name, y = cs), colour = viridis(10)[10], size = 1) +
  lims(y = c(-2, 2)) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(label = C_labs2) +
  facet_wrap(~ prop,
             labeller = labeller(prop = prop_labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_markdown(size = 12),
        strip.text = element_markdown(size = 11),
        strip.text.y = element_markdown(size = 12),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1)
  )

ggsave(
  filename = "./figures/replicate_data_reductionC.png",
  plot = C_procession_plot_lin,
  device = "png",
  dpi = 300,
  width = 8,
  height = 8,
  units = "in"
)

ggsave(
  filename = "./figures/replicate_data_reductionC.svg",
  plot = C_procession_plot_lin,
  device = "svg",
  width = 8,
  height = 8,
  units = "in"
)
