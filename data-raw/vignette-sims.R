#' Vignette simulations
#' Run offline (outside of CRAN) and save to ./data-raw

p_thresh <- seq(0.95, 0.995, 0.005)
tab <- NULL

for (i in 1:length(p_thresh)) {
  fit_p <- multi_trial(
    sens_true = 0.7,
    spec_true = 0.963,
    prev_true = 0.20,
    endpoint = "sens",
    sens_pg = 0.7,
    spec_pg = NULL,
    prior_sens = c(0.1, 0.1),
    prior_spec = c(0.1, 0.1),
    prior_prev = c(0.1, 0.1),
    succ_sens = p_thresh[i],
    n_at_looks = seq(100, 600, 50),
    n_mc = 10000,
    n_trials = 5000,
    ncores = 8L)

  out <- summarise_trials(fit_p, min_pos = 35, fut = 0.05)
  tab <- rbind(tab, out)

}

power <- multi_trial(
  sens_true = 0.824,
  spec_true = 0.963,
  prev_true = 0.20,
  endpoint = "sens",
  sens_pg = 0.7,
  spec_pg = NULL,
  prior_sens = c(0.1, 0.1),
  prior_spec = c(0.1, 0.1),
  prior_prev = c(0.1, 0.1),
  succ_sens = 0.985,
  n_at_looks = seq(100, 600, 50),
  n_mc = 10000,
  n_trials = 5000,
  ncores = 8L)

save(tab, power, file = "vignettes/vignette-sims.rda")
