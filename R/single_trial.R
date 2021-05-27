#' @title Simulate and analyse a single trial
#'
#' @description A single trial is simulated and analysed up to the final
#'   analysis stage, irrespective of whether it would have been stopped for
#'   early success or expected futility. The output of the trial is handled
#'   elsewhere.
#'
#' @inheritParams multi_trial
#'
#' @noRd
single_trial <- function(
  sens_true, spec_true, prev_true,
  sens_pg, spec_pg,
  prior_sens, prior_spec, prior_prev,
  succ_sens, succ_spec,
  n_at_looks,
  n_mc
) {

  # TODO: allow for NULL sens or spec (not both)
  # TODO: simplify inputs for sens and spec

  n_stages <- length(n_at_looks)

  sim_dat <- simulate_data(n          = max(n_at_looks),
                           n_at_looks = n_at_looks,
                           sens_true  = sens_true,
                           spec_true  = spec_true,
                           prev_true  = prev_true)

  out <- NULL
  for (k in 1:length(n_at_looks)) {

    trial <- analysis(data = sim_dat,
                      k = k,
                      sens_pg = sens_pg,
                      spec_pg = spec_pg,
                      prior_sens = prior_sens,
                      prior_spec = prior_spec,
                      prior_prev = prior_prev,
                      succ_sens = succ_sens,
                      succ_spec = succ_spec,
                      n_at_looks = n_at_looks,
                      n_mc = n_mc)

    out <- rbind(out, trial)

  }

  out$n <- n_at_looks
  rownames(out) <- NULL

  return(out)

}
