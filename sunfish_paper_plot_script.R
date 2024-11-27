# this script is a bit long due to

if (!require("pacman")) install.packages("pacman", repos="http://cran.r-project.org")
pacman::p_load(targets, tidyverse, ggtext, patchwork, forecast, MultinomialStateSpace, viridis, viridisLite, RColorBrewer)

tar_load(c(data_sim_time_c5, data_sim_pulse_c5, data_fit))
tar_load(c(sim_fit_time_sub_c5, sim_fit_pulse_sub_c5))
ss_mod_list <- readRDS("./data/ss_seq_list.rds")

# Data and model plot
sunfish_years <- data_fit$sunfish_years
sunfish_tsample <- data_fit$ss_mod$Tsample
ssm <- ss_mod_list[[5]]
Y_data <- ssm$Y

Y <- matrix(NA, nrow = sunfish_tsample[length(sunfish_tsample)], ncol = ncol(Y_data))
for (i in seq_along(sunfish_tsample)) {
  Y[sunfish_tsample[i], ] <- as.matrix(Y_data[i, ])
}
colnames(Y) <- colnames(Y_data)

sunfish_years_interp <- matrix(NA, nrow = sunfish_tsample[length(sunfish_tsample)], ncol = 1)
for (i in seq_along(sunfish_tsample)) {
  sunfish_years_interp[sunfish_tsample[i], ] <- as.matrix(sunfish_years[i])
}
sunfish_years_interp <- as.matrix(na.interp(sunfish_years_interp))
colnames(sunfish_years_interp) <- "sunfish_years_interp"

propY <- Y/rowSums(Y)
colnames(propY) <- colnames(ssm$Y)

se_y <- ssm$se.y.fitted
colnames(se_y) <- colnames(propY)
se_y <- se_y/rowSums(se_y)

ss_se_y_plot <- cbind(sunfish_years_interp, se_y) %>%
  as_tibble() %>%
  mutate(cat = "se_y") %>%
  pivot_longer(-c(sunfish_years_interp, cat))

propY_plot <- cbind(sunfish_years_interp, propY) %>%
  as_tibble %>%
  mutate(cat = "prop") %>%
  pivot_longer(-c(sunfish_years_interp, cat)) %>%
  mutate(se = ss_se_y_plot$value,
         name = factor(name, levels = c("other", "Pinus", "Tsuga", "Fagus", "Quercus")))

ggplot(propY_plot, aes(x = sunfish_years_interp, y = value, colour = name)) +
  geom_point() +
  geom_path() +
  scale_x_reverse() +
  facet_wrap(~name)


se_upper <- ssm$se.mu.upper.fitted
colnames(se_upper) <- colnames(propY)

ss_se_upper_plot <- cbind(sunfish_years_interp, se_upper) %>%
  as_tibble %>%
  mutate(cat = "se_upper") %>%
  pivot_longer(-c(sunfish_years_interp, cat))

se_lower <- ssm$se.mu.lower.fitted
colnames(se_lower) <- colnames(propY)

ss_se_lower_plot <- cbind(sunfish_years_interp, se_lower) %>%
  as_tibble %>%
  mutate(cat = "se_lower") %>%
  pivot_longer(-c(sunfish_years_interp, cat))

mu <- ssm$mu
colnames(mu) <- colnames(propY)

ss_mu_plot <- cbind(sunfish_years_interp, mu) %>%
  as_tibble %>%
  mutate(cat = "mu") %>%
  pivot_longer(-c(sunfish_years_interp, cat)) %>%
  mutate(se_upper = ss_se_upper_plot$value,
         se_lower = ss_se_lower_plot$value,
         name = factor(name, levels = c("other", "Pinus", "Tsuga", "Fagus", "Quercus")))


names_list <- c(
  other="other",
  Pinus="_Pinus spp._",
  Tsuga="_Tsuga spp._",
  Betula="_Betula spp._",
  Quercus="_Quercus spp._",
  Fagus="_Fagus spp._"
)

all_spp <- bind_rows(ss_mu_plot, propY_plot)

all_spp_plot <- ggplot(all_spp, aes(x = sunfish_years_interp, y = value,
                                        ymin = se_lower, ymax = se_upper,
                                        colour = cat, fill = cat)) +
  geom_point(size = 1) +
  geom_ribbon(alpha= 0.5, colour = NA) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~name, labeller = as_labeller(names_list), nrow = 1) +
  scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 4)) +
  scale_colour_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired", labels = c("Fitted", "Observed")) +
  coord_flip(ylim = c(0, 1)) +
  guides(color = "none") +
  labs(x = "Time (ybp)", y = "Relative abundances", fill = NULL) +
  theme_minimal() +
  theme(
    strip.text = element_markdown(size = 9),
    text = element_text(size = 9),
    axis.text.x = element_text(),
    legend.position = "bottom"
  )


prop_plot <- ggplot(propY_plot, aes(x = sunfish_years_interp, y = value)) +
  geom_area(colour = "grey90") +
  geom_col() +
  scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
  coord_flip() +
  # ylim(0, 0.5) +
  labs(y = "Relative abundances", x = "Time (ybp)") +
  facet_wrap(~name, labeller = as_labeller(names_list)) +
  theme_minimal() +
  theme(
    strip.text = element_markdown(size = 9),
    text = element_text(size = 9),
  )

prop_plot + all_spp_plot

ggsave(plot = all_spp_plot, filename = "../../mss-paper-resources/images/palaeo_fitted_plot.svg", width = 6.3, height = 6.2, units = "in")


# Data and model plot
ssm_sim <- sim_fit_pulse_sub_c5[[1]][[7]][[1]]
sim_X <- data_sim_pulse_c5[[1]][2, 7]$X[11:150]
Y_sim <- data_sim_pulse_c5[[1]][1, 7]$Y[11:150, ]
Y_sim[-sunfish_tsample, ] <- NA
colnames(Y_sim) <- c("y1", "y2", "y3")

sim_time <- 1:138

propY_sim <- Y_sim/rowSums(Y_sim)
colnames(propY_sim) <- colnames(Y_sim)

se_y_sim <- ssm_sim$se.y.fitted
colnames(se_y_sim) <- colnames(propY_sim)
se_y_sim <- se_y_sim/rowSums(se_y_sim)

ss_se_y_plot_sim <- cbind(sim_time, se_y_sim) %>%
  as_tibble() %>%
  mutate(cat = "se_y") %>%
  pivot_longer(-c(sim_time, cat))

propY_plot_sim <- cbind(sim_time, propY_sim[1:138, ]) %>%
  as_tibble() %>%
  mutate(cat = "prop") %>%
  pivot_longer(-c(sim_time, cat)) %>%
  mutate(se = ss_se_y_plot_sim$value,
         name = factor(name, levels = c("y1", "y2", "y3")))

ggplot(propY_plot_sim, aes(x = sim_time, y = value, colour = name)) +
  geom_point() +
  geom_path() +
  scale_x_reverse() +
  facet_wrap(~name)

se_upper_sim <- ssm_sim$se.mu.upper.fitted
colnames(se_upper_sim) <- colnames(propY_sim)

ss_se_upper_plot_sim <- cbind(sim_time, se_upper_sim) %>%
  as_tibble %>%
  mutate(cat = "se_upper") %>%
  pivot_longer(-c(sim_time, cat))

se_lower_sim <- ssm_sim$se.mu.lower.fitted
colnames(se_lower_sim) <- colnames(propY_sim)

ss_se_lower_plot_sim <- cbind(sim_time, se_lower_sim) %>%
  as_tibble %>%
  mutate(cat = "se_lower") %>%
  pivot_longer(-c(sim_time, cat))

mu_sim <- ssm_sim$mu
colnames(mu_sim) <- colnames(propY_sim)

ss_mu_plot_sim <- cbind(sim_time, mu_sim) %>%
  as_tibble() %>%
  mutate(cat = "mu") %>%
  pivot_longer(-c(sim_time, cat)) %>%
  mutate(se_upper = ss_se_upper_plot_sim$value,
         se_lower = ss_se_lower_plot_sim$value,
         name = factor(name, levels = c("y1", "y2", "y3")))


names_list_sim <- c(
  y1="sp. 1",
  y2="sp. 2",
  y3="sp. 3"
)


all_spp_sim <- bind_rows(ss_mu_plot_sim, propY_plot_sim)

all_spp_plot_sim <- ggplot(all_spp_sim, aes(x = sim_time, y = value,
                                    ymin = se_lower, ymax = se_upper,
                                    colour = cat, fill = cat)) +
  geom_point(size = 1) +
  geom_ribbon(alpha= 0.5, colour = NA) +
  geom_line(linewidth = 0.3) +
  facet_wrap(~name, labeller = as_labeller(names_list_sim)) +
  scale_colour_brewer(palette="Paired") +
  scale_fill_brewer(palette="Paired", labels = c("Fitted", "Observed")) +
  coord_flip() +
  guides(color = "none") +
  labs(x = NULL, y = "Relative abundances", fill = NULL) +
  theme_minimal() +
  theme(
    text = element_text(size = 10),
    axis.text.x = element_text(),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  )


x_data <- tibble(sim_time = 1:140,
                 X = sim_X)

x_plot <- ggplot(x_data, aes(x = sim_time, y = X)) +
  geom_line(linewidth = 0.3) +
  coord_flip() +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
  scale_y_continuous(breaks = c(0,1)) +
  labs(y = "X", x = "Simulation time") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    text = element_text(size = 10),
  )

x_plot + all_spp_plot_sim + plot_layout(widths = c(.1, .9))
sim_plot <- x_plot + all_spp_plot_sim + plot_layout(widths = c(.1, .9))
ggsave(plot = sim_plot, filename = "../../mss-paper-resources/images/sim_plot.svg", width = 6.3, height = 6.5, units = "in")



sim_list <- list(
  sim_fit_time_sub_c5,
  sim_fit_pulse_sub_c5
)
names(sim_list) <- c("time_subsam", "pulse_subsam")

ssms <- lapply(sim_list, \(sce) {
  lapply(sce, \(rep) {
    ss_without_error <- lapply(rep, \(ssm){
      ss <- ssm[[1]]
    })
    ss_mods <- ss_without_error[lapply(ss_without_error, class) == "multinomialSS"]
    print(length(ss_mods))
    ss_mods
  })
})


ssms <- unlist(ssms, recursive = F, use.names = T)

pars <- lapply(ssms, \(case) {
  all_pars <- lapply(case, \(rep) {
    v.pars <- names(rep$par)[grep(names(rep$par), pattern = "v.")]
    par.names.noshow.se <- c("dispersion","sigma", v.pars)
    par.show.se <- rep$par[!is.element(names(rep$par), par.names.noshow.se)]
    par.all <- c(par.show.se, logLik = rep$logLik, opt.convergence = rep$opt.convergence)
  })
  bind_rows(all_pars, .id = "rep")
})

pars_wide_bs <- bind_rows(pars, .id = "scenario") %>%
  separate(scenario, into = c("scenario", "case"), sep = "(?<=\\D)(?=\\d)", convert = F, remove = FALSE) %>%
  mutate(across(where(is.character), forcats::fct)) %>%
  select(c(scenario, case, rep, opt.convergence), starts_with("x"))


pars_long_bs <- pars_wide_bs %>%
  pivot_longer(cols = -c(scenario, case, rep, opt.convergence)) %>%
  mutate(
   bs = case_when(case %in% c(1:4) & name == "x1.y2" ~ 0.25,
                  case %in% c(1:4) & name == "x1.y3" ~ -0.5,
                  case %in% c(5:6) & name %in% c("x1.y2", "x1.y3") ~ 0,
                  .default = NA)
  )

pars_long_bs %>%
  group_by(scenario, case, name, bs) %>%
  summarise(med = median(value),
            quant25 = quantile(value, probs = 0.25),
            quant75 = quantile(value, probs = 0.75)) %>%
  print(n = 24)

h_line_dat <- tibble(case = factor(1:6),
                         intb1 = c(rep(0.25, 4), rep(0, 2)),
                         intb2 = c(rep(-0.5, 4), rep(0, 2))
                         )

case_labs <- c(
  `1`="C<sub>0,</sub> V<sub>0</sub>",
  `2`="C<sub>0,</sub> V<sub>cor</sub>",
  `3`="C<sub>comp,</sub> V<sub>0</sub>",
  `4`="C<sub>comp,</sub> V<sub>cor</sub>",
  `5`="C<sub>comp,</sub> V<sub>cor,</sub> B<sub>null</sub>",
  `6`="C<sub>0,</sub> V<sub>cor,</sub> B<sub>null</sub>"
)

scenario_labs <- c(
  time_subsam = "X<sub>lin</sub>",
  pulse_subsam="X<sub>dis</sub>"
)


X_labs <- c("x1.y2" = "B<sub>1, 2</sub>",
            "x1.y3" = "B<sub>1, 3</sub>")


B_procession_plot <- ggplot(data = pars_long_bs) +
  geom_violin(aes(x = name, y = value, fill = scenario, alpha = 0.5),
              draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3) +
  geom_point(aes(x = name, y = bs), colour = viridis(10)[10], size = 0.6) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = brewer.pal(11, "BrBG")[c(3, 10)]) +
  scale_x_discrete(labels = X_labs) +
  facet_grid(vars(scenario), vars(case), scales = "free",
             labeller = labeller(case = case_labs,
                                 scenario = scenario_labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_markdown(angle = 45),
        strip.text = element_markdown(size = 8),
        strip.text.y = element_markdown(size = 12),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1)
        )

ggsave(plot = B_procession_plot, filename = "../../mss-paper-resources/images/B_plot.svg", width = 6.3, height = 4, units = "in")


pars_wide_cs <- bind_rows(pars, .id = "scenario") %>%
  separate(scenario, into = c("scenario", "case"), sep = "(?<=\\D)(?=\\d)", convert = F, remove = FALSE) %>%
  mutate(across(where(is.character), forcats::fct)) %>%
  select(c(scenario, case, rep, opt.convergence), starts_with("sp."))


pars_long_cs <- pars_wide_cs %>%
  pivot_longer(cols = -c(scenario, case, rep, opt.convergence)) %>%
  mutate(
    cs = case_when(name == "sp.y1.y1" ~ 0.5,
                   name == "sp.y2.y2" ~ 0.9,
                   name == "sp.y3.y2" ~ 0,
                   name == "sp.y2.y3" ~ -0.5,
                   name == "sp.y3.y3" ~ 0.6,
                   name == "sp.y2.y3" & case == 1 ~ 0,
                   name == "sp.y2.y3" & case == 2 ~ 0,
                   name == "sp.y2.y3" & case == 6 ~ 0,
                   .default = NA),
    cs = case_when(name == "sp.y2.y3" & case == 1 ~ 0,
                   name == "sp.y2.y3" & case == 2 ~ 0,
                   name == "sp.y2.y3" & case == 6 ~ 0,
                   .default = as.numeric(cs)),
    name = fct_relevel(name, "sp.y1.y1", "sp.y2.y2", "sp.y3.y3",
                             "sp.y2.y3", "sp.y3.y2")
  )

pars_long_group_c <- pars_long_cs %>%
  group_by(scenario, case, name, cs) %>%
  summarise(med = median(value),
            quant25 = quantile(value, probs = 0.25),
            quant75 = quantile(value, probs = 0.75))


pars_long_cs %>%
  mutate(
    cs = case_when(
                   name == "sp.y2.y3" & case == 1 ~ 0,
                   name == "sp.y2.y3" & case == 2 ~ 0,
                   name == "sp.y2.y3" & case == 6 ~ 0,
                   .default = NA))

pars_long_cs_reordered <- pars_long_cs %>%
  mutate(
    case_reord = case_when(
      case == 3 ~ 6,
      case == 6 ~ 3,
      .default = as.numeric(case)),
    case_reord = as.factor(case_reord))

pars_long_cs_reordered <- pars_long_cs %>%
  mutate(case = fct_relevel(case, "1", "2", "6", "3", "4", "5"))


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

C_procession_plot_lin <- ggplot(data = pars_long_cs_reordered %>% filter(scenario == "time_subsam")) +
  geom_violin(aes(x = name, y = value, fill = scenario, alpha = 0.5),
              draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3) +
  geom_point(aes(x = name, y = cs), colour = viridis(10)[10], size = 0.6) +
  lims(y = c(-2, 2)) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = brewer.pal(11, "BrBG")[c(3, 10)]) +
  scale_x_discrete(label = C_labs2) +
  facet_wrap(~ case, scales = "free",
             labeller = labeller(case = case_labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_markdown(angle = 45),
        strip.text = element_markdown(),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1)
  )

C_procession_plot_lin
ggsave(plot = C_procession_plot_lin, filename = "../../mss-paper-resources/images/Clin_plot.svg", width = 6.3, height = 4.8, units = "in")



C_procession_plot_dis <- ggplot(data = pars_long_cs_reordered %>% filter(scenario == "pulse_subsam")) +
  geom_violin(aes(x = name, y = value, fill = scenario, alpha = 0.5),
              draw_quantiles = c(0.25, 0.5, 0.75), linewidth = 0.3) +
  geom_point(aes(x = name, y = cs), colour = viridis(10)[10], size = 0.6) +
  lims(y = c(-2, 2)) +
  labs(x = NULL, y = NULL) +
  scale_fill_manual(values = brewer.pal(11, "BrBG")[c(10)]) +
  scale_x_discrete(label = C_labs2) +
  facet_wrap(~ case, scales = "free",
             labeller = labeller(case = case_labs)) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_markdown(angle = 45),
        strip.text = element_markdown(),
        panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.1)
  )

C_procession_plot_dis
ggsave(plot = C_procession_plot_dis, filename = "../../mss-paper-resources/images/Cdis_plot.svg", width = 6.3, height = 4.8, units = "in")


###
blanding_ll <- read_csv("./bryan_data/Blanding_LLrecon_041416_ShumanBurrell17.csv")
blanding_ll[ ,2:4] <- -blanding_ll[ ,2:4]
X_ll <- data.frame(y = blanding_ll$LakeElev.cm.below.mod, x = blanding_ll$Age)
X_ll_gam <- mgcv::gam(y ~ s(x, bs = "bs", k = nrow(X_ll)), method = "REML", data =  X_ll[nrow(X_ll):1, , drop = F])
pred <- predict(X_ll_gam, newdata = data.frame(x = sunfish_counts_with_dates$sunfish_years[63:69]), type = "link", se.fit = TRUE)

upr <- pred$fit + (2 * pred$se.fit)
lwr <- pred$fit - (2 * pred$se.fit)
X_ll_gam$family$linkinv(upr)


sunfish_counts_with_dates <- readRDS("../data/sunfish_counts_with_dates.rds")

blanding_gam <- data.frame(
  Age = sunfish_counts_with_dates$sunfish_years[63:69],
  LakeElev.cm.below.mod = pred$fit,
  UpperBound = 0,
  LowerBound = 0
)

blanding_ll_gam <- bind_rows(blanding_ll, blanding_gam, .id = "mod")



sunfish_spp_long <- pivot_longer(sunfish_counts_with_dates, cols = -sunfish_years, names_to = "variablename") %>%
  mutate(variablename = replace(variablename, stringr::str_detect(variablename, "Pinus.*"), "Pinus"),
         variablename = replace(variablename, stringr::str_detect(variablename, "Acer*"), "Acer")) %>%
  group_by(variablename, sunfish_years) %>%
  summarise(value = sum(value), .groups = 'keep')

sunfish_PT <- sunfish_spp_long %>%
  filter(variablename %in% c("Pinus", "Tsuga", "Betula", "Fagus", "Quercus"))

sunfish_other <- sunfish_spp_long %>%
  filter(!variablename %in% c("Pinus", "Tsuga", "Betula", "Fagus", "Quercus")) %>%
  mutate(variablename = "other") %>%
  group_by(variablename, sunfish_years) %>%
  summarise(value = sum(value), .groups='keep')

sunfish_spp_all <- bind_rows(sunfish_other, sunfish_PT) %>%
  group_by(sunfish_years) %>%
  mutate(pollencount = sum(value, na.rm = TRUE)) %>%
  group_by(variablename) %>%
  mutate(prop = value / pollencount) %>%
  arrange(desc(sunfish_years)) %>%
  mutate(variablename = forcats::fct(variablename,
                                     levels = c("other", "Betula", "Fagus", "Pinus", "Quercus", "Tsuga")))


names_list5 <- c(
  other="other",
  Pinus="_Pinus_",
  Tsuga="_Tsuga_",
  Betula="_Betula_",
  Quercus="_Quercus_",
  Fagus="_Fagus_"
)

spp_plot <- ggplot(sunfish_spp_all %>% filter(!variablename %in% c("Betula")), aes(x = sunfish_years, y = value)) +
  geom_area(fill = "grey20") +
  geom_segment(data = sunfish_spp_all %>% filter(!variablename %in% c("Betula")),
               aes(x = sunfish_years, xend = sunfish_years,
                   y = 0, yend = value), colour = "grey30", linewidth = 0.6) +
  scale_x_reverse(breaks = scales::breaks_pretty(n = 6)) +
  coord_flip() +
  labs(y = "Pollen counts", x = "Time (ybp)") +
  facet_wrap(~variablename, labeller = as_labeller(names_list5),
             nrow = 1) +
  theme_minimal() +
  theme(
    strip.text = element_markdown(size = 9),
    text = element_text(size = 9),
  )
spp_plot


lake_level_plot <-
ggplot(blanding_ll, aes(x = Age, y = LakeElev.cm.below.mod,
                        ymin = LowerBound, ymax = UpperBound)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(alpha = 0.3) +
  scale_x_reverse(limits = c(14000, 0), breaks = seq(0, 14000, 2000)) +
  coord_flip() +
  labs(y = "Lake elevation \n (cm relative to modern)") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    text = element_text(size = 9)
  )

spp_plot + lake_level_plot + plot_layout(nrow = 1, widths = c(1, 0.3), heights = c(0.9, 0.1))
palaeo_plot <- spp_plot + lake_level_plot + plot_layout(nrow = 1, widths = c(1, 0.3), heights = c(0.9, 0.1))

ggsave(plot = palaeo_plot, filename = "../../mss-paper-resources/images/palaeo_plot.svg", width = 6.3, height = 6.2, units = "in")





