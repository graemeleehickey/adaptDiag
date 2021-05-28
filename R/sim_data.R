#' @title Simulate a complete dataset
#'
#' @description Simulate a single dataset for a diagnostic test under the
#'   assumption of a gold standard reference and potentially adaptive
#'   group-sequential design.
#'
#' @param n integer. Maximum sample size for the trial.
#' @inheritParams multi_trial
#'
#' @importFrom stats rbinom
#'
#' @return Simulated data frame with columns for:
#'
#' - Subject ID
#' - Reference test outcome (\code{1 =} positive, \code{0 =} negative)
#' - Diagnostic test outcome (\code{1 =} positive, \code{0 =} negative)
#' - Stage of the adaptive analysis (integer)
#'
#' @noRd
simulate_data <- function(n, n_at_looks, sens_true, spec_true, prev_true) {

  id <- 1:n
  stage <- vector(length = n)
  if (!is.null(n_at_looks)) {
    for (i in length(n_at_looks):1) {
      stage[id <= n_at_looks[i]] <- i
    }
  } else {
    stage <- rep(1, n)
  }
  ref_pos  <- rbinom(n, 1, prev_true)
  test_pos <- rbinom(n, 1, ifelse(ref_pos == 1, sens_true, 1 - spec_true))
  dat <- data.frame(id, ref_pos, test_pos, stage)

  return(dat)

}
