#' @title Summarise results of multiple simulated trials to give the operating
#'   characteristics
#'
#' @param data list. Output from the \code{\link{multi_trial}} function.
#' @param min_pos integer. The minimum number of reference positive cases before
#'   stopping is allowed. Default is \code{min_pos = 1}.
#' @param fut scalar. A probability threshold at which the posterior predictive
#'   probability of eventual success is compared to. If the probability is less
#'   than \code{fut}, the trial stops for binding futility. Default is \code{fut
#'   = 0}, which corresponds to no stopping for futility.
#'
#' @return A data frame of row length 1, with the following columns:
#'
#' \itemize{
#'     \item{\code{power}:} Power is defined as the proportion of trials that
#'     result in success, irrespective of whether it is an early stop for
#'     success or not. Trials that stop for futility, but which subsequently go
#'     on to be successful, are not considered as a success. In other words, the
#'     futility decision is binding, and in practice, if a trial triggered a
#'     futility rule, the sponsor would not see the eventual outcome if the
#'     trial were to continue enrolling. When the performance goals are set
#'     equal to the respective true values, the power returned is the type I
#'     error.
#'     \item{\code{stop_futility}:} The proportion of trials that stopped early
#'     for expected futility.
#'     \item{\code{n_avg}:} The average sample size for trials at the stage they
#'     stopped.
#'     \item{\code{sens}:} The average sensitivity for trials at the stage they
#'     stopped.
#'     \item{\code{spec}:} The average specificity for trials at the stage they
#'     stopped.
#'     \item{\code{mean_pos}:} The average number of reference positive cases
#'     for trials at the stage they stopped.
#' }
#'
#' @export
#'
#' @examples
#' data <- multi_trial(
#'     sens_true = 0.9,
#'     spec_true = 0.95,
#'     prev_true = 0.1,
#'     endpoint = "both",
#'     sens_pg = 0.8,
#'     spec_pg = 0.8,
#'     prior_sens = c(1, 1),
#'     prior_spec = c(1, 1),
#'     prior_prev = c(1, 1),
#'     succ_sens = 0.95,
#'     succ_spec = 0.95,
#'     n_at_looks = c(200, 400, 600, 800, 1000),
#'     n_mc = 10000,
#'     n_trials = 20,
#'     ncores = 1
#'     )
#'
#' summarise_trials(data, fut = 0.05, min_pos = 10)
summarise_trials <- function(data, min_pos = 1, fut = 0) {

  sims <- data$sims
  args <- data$args
  args$fut <- fut
  args$min_pos <- min_pos

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

  print(with(out, table(decision, n, useNA = "i")))
  cat("\n")

  return(summ)

}

#' @title Evaluate a single trial
#'
#' @description For a single trial that has been fully enumerated irrespective
#'   of whether stopping rules were triggered, the \code{evaluate_trial}
#'   function will process the output alongside the probability success
#'   threshold(s), the futility futility threshold, and any other constraints to
#'   indicate the decision taken and the stopping time.
#'
#' @param x data frame. Simulated trial data from the \code{\link{multi_trial}}
#'   function.
#' @param args list. Arguments passed to the \code{\link{multi_trial}} function
#'   with additional constrained concatenated.
#'
#' @return Data frame of simulated trials, with 1 row per trial, which extracts
#'   the row that met the stopping criteria. In addition to the columns returned
#'   by \code{\link{multi_trial}}, an additional column is appended named
#'   \code{decision}. This can take the following values:
#'
#' \itemize{
#'   \item{\code{stop for futility}: Expected futility threshold was crossed.}
#'   \item{\code{early win}: Early success threshold was crossed.}
#'   \item{\code{late win}: Success threshold was crossed, but only at the final
#'   look.}
#'   \item{\code{no stopping}:} The trial progressed all the way to the final
#'   sample size look and did not trigger any stopping rules or other
#'   constraints.
#'   \item{\code{no stopping}:} The trial progressed all the way to the final
#'   sample size look and did not trigger any stopping rules, however the number
#'   of reference positive cases was less than the minimum constrain
#'   (\code{min_pos}).
#' }
#'
#' @noRd
evaluate_trial <- function(x, args) {

  n_looks <- nrow(x)
  pass <- 0
  futile <- 0

  for (i in 1:n_looks) {
    if (args$min_pos > (x$tp[i] + x$fn[i])) {
      if (i == n_looks) {
        decision <- "no stopping - insufficient positive cases"
        break
      } else {
        next
      }
    }
    if (args$endpoint == "both") {
      if ((x$pp_sens[i] >= args$succ_sens) & (x$pp_spec[i] >= args$succ_spec)) {
        decision <- ifelse(i < n_looks, "early win", "late win")
        break
      } else if ((x$ppp_succ_both[i] < args$fut) & (i < n_looks)) {
        decision <- "stop for futility"
        break
      } else {
        decision <- "no stopping"
      }
    } else if (args$endpoint == "sens") {
      if (x$pp_sens[i] >= args$succ_sens) {
        decision <- ifelse(i < n_looks, "early win", "late win")
        break
      } else if ((x$ppp_succ_sens[i] < args$fut) & (i < n_looks)) {
        decision <- "stop for futility"
        break
      } else {
        decision <- "no stopping"
      }
    } else {
      if (x$pp_spec[i] >= args$succ_spec) {
        decision <- ifelse(i < n_looks, "early win", "late win")
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
