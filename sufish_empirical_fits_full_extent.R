if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(tidyverse, multinomialTS, forecast, mgcv, future, future.apply, furrr, ggtext, patchwork)

# Sometimes re-fitting the model with different starting conditions gives
# a better fit. This function refits the mnTS model using the output of
# the previous model fit. This is ok to do because it is only adjusting
# starting conditions.
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


# The data for Sunfish Pond are from Johnson et al., in prep
# a bit of wrandling is required
sunfish_pollen <- read_csv("./data/Sunfish_585_pollen_counts_bchron062722B_cleaned.csv")
sunfish_ll <- read_csv("./data/Sunfish_LL_Recon_Feb22.csv")

# Pull pollen columns only
sunfish_pollen_only <- sunfish_pollen[-c(1:8)]
# Remove spp with less than 20 count in whole core (recommended by lead author)
filter_5_names <- names(which(colSums(sunfish_pollen_only) <= 20))

# Remove some aquatic/semi-aquatic spp (recommended by lead author)
aquatic_spp <- c(
  "Brasenia", "Nuphar",
  "Nym", "Nym.cell",
  "Potamoge","Myrica",
  "C.stolon", "Botrich",
  "Botrioc")

filter_5_names[which(aquatic_spp %in% filter_5_names)]
spp_remove <- c(aquatic_spp, setdiff(filter_5_names, aquatic_spp))
colSums(sunfish_pollen_only[ ,spp_remove])

# clean things up to get depth, age, and spp counts
sunfish_pollen_clean <- sunfish_pollen |>
  select(-c(all_of(spp_remove), DepthCM.Vania, Depth.predict,
             `_2.5`, `_10`, `_90`, `_975`)) |>
  rename("depth" = "Depth.Correct",
          "age" = "Age.bchron062722B")

dim(sunfish_pollen_clean)


# Distribute undifined pine counts between P.strobus and P.diplo based on their
# relative proportions (recommended by lead author)
sunfish_pollen_assigned <- sunfish_pollen_clean |>
  mutate(
    total_pollen = P.strobu + P.diplo,
    P.strobu = if_else(total_pollen == 0, P.strobu, round(P.strobu + (Pinus * (P.strobu / total_pollen)))),
    P.diplo = if_else(total_pollen == 0, P.diplo, round(P.diplo + (Pinus * (P.diplo / total_pollen))))
  ) |>
  select(-c(total_pollen, Pinus))

sum(colSums(sunfish_pollen_clean[ ,c("Pinus", "P.strobu", "P.diplo")]))
sum(colSums(sunfish_pollen_assigned[ ,c("P.strobu", "P.diplo")]))
colnames(sunfish_pollen_assigned)

sunfish_spp_long <- pivot_longer(sunfish_pollen_assigned, cols = -c(age, depth),
                                 names_to = "variablename")

# P.diplo now being aggregated with other
target_spp <- c("P.strobu", "Tsuga", "Fagus", "Quercus", "Betula")
colSums(sunfish_pollen_assigned[ ,target_spp])
colSums(sunfish_pollen_assigned[, !colnames(sunfish_pollen_assigned) %in% c(target_spp, "age", "depth"), drop = FALSE])


# Pull out target species
sunfish_target_spp <- sunfish_spp_long |>
  filter(variablename %in% target_spp)

# Group all other spp
sunfish_other <- sunfish_spp_long |>
  filter(!variablename %in% target_spp) |>
  mutate(variablename = "other") |>
  group_by(variablename, age, depth) |>
  summarise(value = sum(value), .groups='keep')

sunfish_grouped_long <- bind_rows(sunfish_other, sunfish_target_spp) |>
  mutate(variablename = fct(variablename, levels = c("other", target_spp)))

# Plot spp
sunfish_grouped_long |>
  arrange(age) |>
  ggplot(aes(x = age, y = value)) +
    geom_point() +
    geom_line() +
    scale_x_reverse() +
    facet_wrap(~ variablename, scales = "free")

# Stack and pivot to wide
# arrange from oldest to youngest
sunfish_spp_wide <- bind_rows(sunfish_other, sunfish_target_spp) |>
  pivot_wider(id_cols = age, names_from = variablename, values_from = value) |>
  arrange(desc(age))

# Create bins
bin_width <- 100
sunfish_bins <-
  cut(
    sunfish_spp_wide$age,
    breaks = seq(
      from = min(sunfish_spp_wide$age),
      to = max(sunfish_spp_wide$age + bin_width),
      by = bin_width
    ), include.lowest = T, labels = F)

# Check if data fall within a bin and need to be summed
diff(sunfish_bins)

nrow(sunfish_spp_wide)

# Clipping to Holocene
sunfish_spp_binned <- bind_cols(bins = sunfish_bins, sunfish_spp_wide) |>
  group_by(bins) |> # Group the data by the bins so that we calculate per time bin
  summarise(across(c(other, all_of(target_spp)), \(x) sum(x ,na.rm = TRUE)),
            age = mean(age, na.rm = TRUE)) |>
              arrange(desc(age))
  # filter(bins <= 122) # for clipping to holocene only


ll_bins <-
  cut(
    sunfish_ll$Age,
    breaks = seq(
      from = min(sunfish_spp_wide$age),
      to = max(sunfish_spp_wide$age + bin_width),
      by = bin_width
    ), include.lowest = T, labels = F)

ll_mean <- cbind(bins = ll_bins, sunfish_ll) |>
  drop_na(bins) |>
  select(bins, LakeElev.cm.below.mod) |>
  group_by(bins) |>
  summarise(mean_ll = mean(LakeElev.cm.below.mod)) |>
  arrange(desc(bins))

# join up the data by bin
sunfish_all <- sunfish_spp_binned |>
  full_join(ll_mean) |>
  arrange(desc(bins))

# Final time window does not have a lake-level observation
which(is.na(sunfish_all$mean_ll))

# using simple linear interpolation to fill final bin
sunfish_all_interp <- sunfish_all |>
  mutate(across(c(mean_ll), forecast::na.interp))

Y <- sunfish_all_interp |>
  select(other, all_of(target_spp)) |>
  as.matrix()


Tsample <- which(rowSums(Y) != 0)

# set up X
X <- sunfish_all_interp |>
  select(mean_ll) |>
  as.matrix() |>
  scale() # necessary when using multiple X variabes

matplot(X, type = 'l')

# set up parameters
p <- ncol(X) + 1 # Number of independent variables plus intercept
n <- ncol(Y)

V.fixed = diag(n) # Covariance matrix of environmental variation in process eq
# V.fixed = matrix(NA, n, n) # other fomrms of V possible
# V.fixed[1] = 1

B.fixed <- matrix(c(rep(0,p),rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0,p),rep(.01, (n - 1) * p)), p, n)

# Run glmm
glmm_mod <- multinomialTS::mnGLMM(Y = Y[which(rowSums(Y) != 0),],
                                  X = X[which(rowSums(Y) != 0), ,drop = F],
                                  B.start = B.start, B.fixed = B.fixed,
                                  V.fixed = V.fixed)
summary(glmm_mod)

# set up parameters for mnTS
B0.start <- glmm_mod$B[1, , drop = F]
B.start <- glmm_mod$B[2:p, , drop = F]

sigma.start <- glmm_mod$sigma

V.fixed = matrix(NA, n, n) # Covariance matrix of environmental variation in process eq
V.fixed[1] = 1

V.start <- glmm_mod$V
# V.start <- diag(diag(V.start))

B.fixed <- matrix(NA, ncol(X), n)
B.fixed[,1] <- 0
B0.fixed = matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)

# Set-up C without interactions
C.start.diag = .5 * diag(n)
C.fixed.diag <- C.start.diag
C.fixed.diag[C.fixed.diag != 0] <- NA

# Model with no interactions
start_time <- Sys.time()
mnTS_mod <- mnTS(Y = Y[which(rowSums(Y) != 0),],
                 X = X, Tsample = Tsample,
                 B0.start = B0.start, B0.fixed = B0.fixed,
                 B.start = B.start, B.fixed = B.fixed,
                 C.start = C.start.diag, C.fixed = C.fixed.diag,
                 V.start = V.start, V.fixed = V.fixed,
                 dispersion.fixed = 1, maxit.optim = 1e+6)
# maxit.optim is the max number of iterations the optimiser will complete before stopping.
# increase maxit.optim if the model needs a lot of time to fit.
end_time <- Sys.time()
end_time - start_time
summary(mnTS_mod)
saveRDS(mnTS_mod, "./empirical_results/mnts_mod.rds")

# res <- multinomialTS::boot(mnTS_mod, 1500)
# saveRDS(res, "./empirical_results/sunfish_mnts_ll_bootstraps_1000.rds")

# With interactions PQF ---------------------------------------------------

C.start.pqf = .5 * diag(n)
C.start.pqf[5,2]=C.start.pqf[4,2]=C.start.pqf[2,5]=C.start.pqf[2,4] = -0.05

C.fixed.pqf <- C.start.pqf
C.fixed.pqf[C.fixed.pqf != 0] <- NA


start_time <- Sys.time()
mnTS_mod_pqf <- mnTS(Y = Y[Tsample, ],
                     X = X, Tsample = Tsample,
                     B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                     B.start = mnTS_mod$B, B.fixed = B.fixed,
                     C.start = C.start.pqf, C.fixed = C.fixed.pqf,
                     V.start = mnTS_mod$V, V.fixed = V.fixed,
                     dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_pqf_refit <- refit_func(mnTS_mod_pqf, 5)
lapply(mnTS_mod_pqf_refit, coef)

# res_pqf <- multinomialTS::boot(mnTS_mod_pqf_refit[[5]], 1500)
# saveRDS(res, "./empirical_results/sunfish_mnts_ll_pqf_bootstraps_1000.rds")

# With interactions PF ----------------------------------------------------

C.start.pf = .5 * diag(n)
C.start.pf[4,2]=C.start.pf[2,4] = -0.05

C.fixed.pf <- C.start.pf
C.fixed.pf[C.fixed.pf != 0] <- NA


start_time <- Sys.time()
mnTS_mod_pf <- mnTS(Y = Y[Tsample, ],
                     X = X, Tsample = Tsample,
                     B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                     B.start = mnTS_mod$B, B.fixed = B.fixed,
                     C.start = C.start.pf, C.fixed = C.fixed.pf,
                     V.start = mnTS_mod$V, V.fixed = V.fixed,
                     dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_pf_refit <- refit_func(mnTS_mod_pf, 5)
lapply(mnTS_mod_pf_refit, coef)

# res_pf <- multinomialTS::boot(mnTS_mod_pf, 1500)
# saveRDS(res, "./empirical_results/sunfish_mnts_ll_pf_bootstraps_1000.rds")

# With interactions qf ----------------------------------------------------

C.start.qf = .5 * diag(n)
C.start.qf[5,4]=C.start.qf[4,5] = -0.05

C.fixed.qf <- C.start.qf
C.fixed.qf[C.fixed.qf != 0] <- NA


start_time <- Sys.time()
mnTS_mod_qf <- mnTS(Y = Y[Tsample, ],
                    X = X, Tsample = Tsample,
                    B0.start = mnTS_mod$B0, B0.fixed = B0.fixed,
                    B.start = mnTS_mod$B, B.fixed = B.fixed,
                    C.start = C.start.qf, C.fixed = C.fixed.qf,
                    V.start = mnTS_mod$V, V.fixed = V.fixed,
                    dispersion.fixed = 1, maxit.optim = 1e+6)

end_time <- Sys.time()
end_time - start_time

mnTS_mod_qf_refit <- refit_func(mnTS_mod_qf, 5)
lapply(mnTS_mod_qf_refit, coef)

# res_qf <- multinomialTS::boot(mnTS_mod_qf, 1500)
# saveRDS(res, "./empirical_results/sunfish_mnts_ll_qf_bootstraps_1000.rds")


# Bootstrap parallel ------------------------------------------------------

mods_list <- list(
  mnTS_mod = mnTS_mod,
  mnTS_mod_pqf_refit = mnTS_mod_pqf_refit[[5]],
  mnTS_mod_pf = mnTS_mod_pf,
  mnTS_mod_qf = mnTS_mod_qf
)


plan(multisession, workers = length(mods_list))
res_par <- furrr::future_map(
  mods_list, multinomialTS::boot, rep = 1500,
  .options = furrr_options(seed = 1984))

saveRDS(res_par, "./empirical_results/sunfish_res_par_ll_bootstraps_1500.rds")


## Bootstrap plotting B ---------------------------------------------------

res_par <- readRDS("./empirical_results/sunfish_res_par_ll_bootstraps_1500.rds")

lapply(res_par, \(x) {
  length(x[[3]])
})

X_names_list <- c(
  mean_ll ="Lake level"
)

mods_boot <- map(res_par, ~ {
  as_tibble(.x[[2]]) |>
    pivot_longer(-c(logLik, opt.convergence))
}) |>
  bind_rows(.id = "mod")

mods_boot_68 <- mods_boot |>
  #  filter(opt.convergence == 0) |>
  group_by(mod, name) |>
  summarise(boot_mean = mean(value),
            boot_sd = sd(value),
            upper_68 = quantile(value, probs = 0.84),
            lower_68 = quantile(value, probs = 0.16)) |>
  mutate(t_scores = boot_mean / boot_sd,
         p_vals = 2 * pnorm(q = abs(t_scores), lower.tail = F),
         sig = p_vals < 0.05)
#

mods_boot_table <- mods_boot_68 |>
  mutate(name = str_replace_all(name,
                                pattern = "y1|y2|y3|y4|y5|y6",
                                replacement = function(x) case_when(
                                  x == "y1" ~ "Other",
                                  x == "y2" ~ "P.strobus",
                                  x == "y3" ~ "Tsuga",
                                  x == "y4" ~ "Fagus",
                                  x == "y5" ~ "Quercus",
                                  x == "y6" ~ "Betula",
                                  TRUE ~ x))) |>
  filter(!str_detect(name, "v."))

mods_boot_68_B <- mods_boot_68 |>
  filter(grepl(paste(names(X_names_list), collapse = "|"), name)) |>
  separate_wider_delim(cols = name, delim = ".", names = c("cov", "name")) |>
  mutate(name = str_replace_all(name,
                                pattern = "y2|y3|y4|y5|y6",
                                replacement = function(x) case_when(
                                  x == "y2" ~ "_P.strobus_",
                                  x == "y3" ~ "_Tsuga_",
                                  x == "y4" ~ "_Fagus_",
                                  x == "y5" ~ "_Quercus_",
                                  x == "y6" ~ "_Betula_",
                                  TRUE ~ x)))

mod_plots_B <- mods_boot_68_B |>
  mutate(panel = "Lake level effect") |>
  ggplot(aes(x = name, y = boot_mean, colour = as_factor(sig))) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
               width = .2, alpha = 0.5) +
  scale_color_manual(name = "Significance", labels = c("> 0.05", "< 0.05"),
                    values = c("#202020", "#d80000")) +
  facet_wrap(~mod) +
  labs(x = "Taxa", y = "Coefficient") +
  theme_bw() +
  theme(
   axis.text = element_markdown(size = 10, angle = 45, hjust = 1),
   axis.title = element_text(size = 10),
   strip.background = element_rect(fill = NA),
   legend.position = "inside",
   legend.position.inside = c(.11, .91),
   legend.text = element_text(size = 8),
   legend.title = element_text(size = 8),
   legend.background = element_rect(fill = NA)
  )


## Bootstrap plotting C ---------------------------------------------------

names_factor <- c(map2_vec(c("Other", target_spp), c("Other", target_spp), \(x, y) paste(x, y, sep = ".")),
                  c("Fagus.Quercus", "Quercus.Fagus", "Betula.Tsuga",
                    "Tsuga.Betula", "Tsuga.Fagus", "Fagus.Tsuga",
                    "P.strobu.Tsuga", "Tsuga.P.strobu", "Tsuga.Quercus", "Quercus.Tsuga",
                    "P.strobu.Fagus", "Fagus.P.strobu", "P.strobu.Quercus", "Quercus.P.strobu"))


mods_boot_68_C <- mods_boot_68 |>
  filter(!grepl(paste(c(names(X_names_list), "v."), collapse = "|"), name),
         !name %in% c("y2", "y3", "y4", "y5", "y6")) |>
  mutate(name = str_remove(name, "sp.")) |>
  # separate_wider_delim(cols = name, delim = ".", names = c("cov", "name")) |>
  mutate(name = str_replace_all(name,
                                pattern = "y1|y2|y3|y4|y5|y6",
                                replacement = function(x) case_when(
                                  x == "y1" ~ "Other",
                                  x == "y2" ~ "P.strobu",
                                  x == "y3" ~ "Tsuga",
                                  x == "y4" ~ "Fagus",
                                  x == "y5" ~ "Quercus",
                                  x == "y6" ~ "Betula",
                                  TRUE ~ x)),
         name = fct(name, levels = names_factor))
         # name = str_replace_all(name, pattern = "P\\.strobu", replacement = "_P.strobus_"),
         # name = str_replace_all(name, pattern = "Tsuga", replacement = "_Tsuga_"))


mod_plots_C <- mods_boot_68_C |>
  # filter(name %in% c("_P.strobus_._Tsuga_", "_Tsuga_._P.strobus_")) |>
  mutate(panel = "Species interactions") |>
  ggplot(aes(x = name, y = boot_mean)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
  geom_errorbar(aes(ymin = lower_68, ymax = upper_68),
                width = .2, alpha = 0.5) +
  labs(x = "Taxa", y = NULL) +
  facet_wrap(~mod) +  # adds the strip
  theme_bw() +
  theme(
    axis.text = element_markdown(size = 10, angle = 45, hjust = 1),
    axis.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    strip.background = element_rect(fill = NA)
  )

sunfish_B_C_ESA <- mod_plots_B + mod_plots_C + plot_layout(width = c(0.7, 0.3))

ggsave(
  filename = "../esa-2025/figures/sunfish_B_C_ESA.svg",
  plot = sunfish_B_C_ESA,
  device = "svg",
  width = 8,
  height = 6,
  units = "in"
)

sunfish_grouped_long_holocene <- sunfish_spp_wide |>
  filter(age <= 12070) |>
  pivot_longer(-c(age))

spp_count_plot <- ggplot(sunfish_grouped_long_holocene, aes(x = age, y = value)) +
  geom_area(colour = "grey90") +
  geom_segment(data = sunfish_grouped_long_holocene,
               aes(x = age, xend = age,
                   y = 0, yend = value), colour = "grey30", linewidth = 0.6) +
  scale_x_reverse() +
  coord_flip() +
  labs(y = "Pollen counts", x = "Time (ybp)") +
  facet_wrap(~name,
             nrow = 1) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
  )

ll_plot <- ggplot(sunfish_ll |> filter(Age <= 12000),
                  aes(x = Age, y = LakeElev.cm.below.mod,
                      ymin = LowerBound, ymax = UpperBound)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(alpha = 0.3) +
  coord_flip() +
  scale_x_reverse()+
  labs(x = NULL, y = "Lake elevation \n (cm relative to modern)") +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_blank()
  )


ll_plot <- sunfish_ll |>
  filter(Age <= 12000) |>
  mutate(panel = "Lake level") |>  # fake facet variable
  ggplot(aes(x = Age, y = LakeElev.cm.below.mod,
             ymin = LowerBound, ymax = UpperBound)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(alpha = 0.3) +
  coord_flip() +
  scale_x_reverse() +
  labs(x = NULL, y = "Lake elevation \n (cm relative to modern)") +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_blank(),
  ) +
  facet_wrap(~panel)

sunfish_esa <- spp_count_plot + ll_plot + plot_layout(widths = c(0.85, 0.15))

ggsave(
  filename = "../images/sunfish_esa.svg",
  plot = sunfish_esa,
  device = "svg",
  width = 10.67,
  height = 6,
  units = "in"
)

