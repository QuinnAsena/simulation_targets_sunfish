# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(crew)
# library(tarchetypes) # Load other packages as needed. # nolint
tar_option_set(error = "null")

controller <- crew::crew_controller_local(
  name = "my_controller",
  workers = 4,
  seconds_idle = 10
)

# Set target options:
tar_option_set(
  packages = c("multinomialTS", "tidyverse", "patchwork", "forecast", "mgcv"),
  format = "rds", # default storage format
  storage = "worker", retrieval = "worker",
  controller = controller
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
future::plan(future.callr::callr)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
library(future)
library(future.callr)
plan(callr)
list(
  tar_target(file, "../data/sunfish_counts_with_dates.rds", deployment = "main"),
  tar_target(data, get_data(file), deployment = "main"),
  tar_target(data_fit, fit_model(data), deployment = "main"),

  tar_target(data_sim_time_c5,
             Y_sim(data_fit,
                   reps = 500,
                   C = matrix(c(.5, 0, 0, 0, .9, 0, 0, -.5, .6), 3, 3),
                   V = matrix(c(1, 0.7, 0.7, 0.7, .6, .5, 0.7, .5, .6), 3, 3),
                   B = c(0, 0.25, -0.5),
                   X = "time"), deployment = "main"),

  tar_target(data_sim_pulse_c5,
             Y_sim(data_fit,
                   reps = 500,
                   C = matrix(c(.5, 0, 0, 0, .9, 0, 0, -.5, .6), 3, 3),
                   V = matrix(c(1, 0.7, 0.7, 0.7, .6, .5, 0.7, .5, .6), 3, 3),
                   B = c(0, 0.25, -0.5),
                   X = "pulse"), deployment = "main"),

  tar_target(sim_fit_time_sub_c5, fit_sim(data_sim_time_c5, data = data, sub_sample = TRUE), deployment = "worker"),
  tar_target(sim_fit_pulse_sub_c5, fit_sim(data_sim_pulse_c5, data = data, sub_sample = TRUE), deployment = "worker"),

  tar_target(boot_time_sub5_c5, boot_sim(sim_fit_time_sub_c5, n_boot = 100, n_sim = 200, case = c(5)), deployment = "worker"),
  tar_target(boot_time_sub6_c5, boot_sim(sim_fit_time_sub_c5, n_boot = 100, n_sim = 200, case = c(6)), deployment = "worker"),

  tar_target(boot_pulse_sub5_c5, boot_sim(sim_fit_pulse_sub_c5, n_boot = 100, n_sim = 200, case = c(5)), deployment = "worker"),
  tar_target(boot_pulse_sub6_c5, boot_sim(sim_fit_pulse_sub_c5, n_boot = 100, n_sim = 200, case = c(6)), deployment = "worker")
)
