#' @title Simulate and analyse a single trial
#'
#' @description A single trial is simulated and analysed up to the final
#'   analysis stage, irrespective of whether it would have been stopped for
#'   early success or expected futility. The output of the trial is handled
#'   elsewhere.
#'
#' @param sens_true scalar. True assumed sensitivity (must be between 0 and 1).
#' @param spec_true scalar. True assumed specificity (scalar 1).
#' @param prev_true scalar. True assumed prevalence as measured by the
#'   gold-standard reference test (must be between 0 and 1).
#' @param sens_pg scalar. Performance goal (PG) for the sensitivity endpoint,
#'   such that the the posterior probability that the PG is exceeded is
#'   calculated. Must be between 0 and 1.
#' @param spec_pg scalar. Performance goal (PG) for the specificity endpoint,
#'   such that the the posterior probability that the PG is exceeded is
#'   calculated. Must be between 0 and 1.
#' @param prior_sens vector. A vector of length 2 with the prior shape
#'   parameters for the sensitivity Beta distribution.
#' @param prior_spec vector. A vector of length 2 with the prior shape
#'   parameters for the specificity Beta distribution.
#' @param prior_prev vector. A vector of length 2 with the prior shape
#'   parameters for the prevalence Beta distribution.
#' @param succ_sens scalar. Probability threshold for the sensitivity to exceed
#'   in order to declare a success. Must be between 0 and 1.
#' @param succ_spec scalar. Probability threshold for the specificity to exceed
#'   in order to declare a success. Must be between 0 and 1.
#' @param fut scalar. TBC.
#' @param n_at_looks vector. Sample sizes for each interim look. The final value
#'   (or only value if no interim looks are planned) is the maximum allowable
#'   sample size for the trial.
#' @param n_mc integer. Number of Monte Carlo draws to use for sampling from the
#'   Beta-Binomial distribution.
#'   
#' @noRd
single_trial <- function(
  sens_true, spec_true, prev_true,
  sens_pg, spec_pg,
  prior_sens, prior_spec, prior_prev, 
  succ_sens, succ_spec,
  fut,
  n_at_looks,
  min_pos,
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
                      n_at_looks = n_at_looks)
    
    out <- rbind(out, trial)
    
  }
  
  out$n <- n_at_looks
  rownames(out) <- NULL
  
  return(out)
  
}