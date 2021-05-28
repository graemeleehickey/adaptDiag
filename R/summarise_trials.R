summarise_trials <- function(data, min_pos, fut) {

  sims <- data$sims
  call <- data$call
  args <- as.list(call)
  args$fut <- fut

  n_looks <- length(args$n_at_looks)

  out <- by(sims, sims$trial, evaluate_trial, args = args)
  out <- do.call("rbind", out)

  summ <- data.frame(
    power         = mean(out$decision %in% c("early win", "late win")),
    stop_futility = mean(out$decision == "stop for futility"),
    n_avg         = mean(out$n),
    sens          = mean(out$sens_hat),
    spec          = mean(out$spec_hat),
    mean_pos      = mean(out$tp + out$fn))

  return(summ)

}

summarise_trials(data, fut = 0.05)

data <- multi_trial(
  sens_true = 0.9,
  spec_true = 0.95,
  prev_true = 0.1,
  endpoint = "both",
  sens_pg = 0.8,
  spec_pg = 0.8,
  prior_sens = c(1, 1),
  prior_spec = c(1, 1),
  prior_prev = c(1, 1),
  succ_sens = 0.95,
  succ_spec = 0.95,
  n_at_looks = c(200, 400, 600, 800, 1000),
  n_mc = 10000,
  n_trials = 200,
  ncores = 8
)

evaluate_trial <- function(x, args, fut) {

  n_looks <- nrow(x)
  pass <- 0
  futile <- 0

  for (i in 1:n_looks) {
    if (args$endpoint == "both") {
      if ((x$pp_sens[i] >= args$succ_sens) & (x$pp_spec[i] >= args$succ_spec)) {
        decision <- ifelse(i > 0 & i < n_looks, "early win", "late win")
        break
      } else if ((x$ppp_succ_both[i] < args$fut) & (i < n_looks)) {
        decision <- "stop for futility"
        break
      } else {
        decision <- "no stopping"
      }
    } else if (args$endpoint == "sens") {
      if (x$pp_sens[i] >= args$succ_sens) {
        decision <- ifelse(i > 0 & i < n_looks, "early win", "late win")
        break
      } else if ((x$ppp_succ_sens[i] < args$fut) & (i < n_looks)) {
        decision <- "stop for futility"
        break
      } else {
        decision <- "no stopping"
      }
    } else {
      if (x$pp_spec[i] >= args$succ_spec) {
        decision <- ifelse(i > 0 & i < n_looks, "early win", "late win")
        break
      } else if ((x$ppp_succ_spec[i] < args$fut) & (i < n_looks)) {
        decision <- "stop for futility"
        break
      } else {
        decision <- "no stopping"
      }
    }
  }

  x$decision <- decision
  x[i, ]

}
