#' @title Simulate and analyse multiple trials
#'
#' @description Multiple trials and simulated and analysed up to the final
#'   analysis stage, irrespective of whether it would have been stopped for
#'   early success or expected futility. The output of the trials is handled
#'   elsewhere.
#'
#' @param sens_true scalar. True assumed sensitivity (must be between 0 and 1).
#' @param spec_true scalar. True assumed specificity (must be between 0 and 1).
#' @param prev_true scalar. True assumed prevalence as measured by the
#'   gold-standard reference test (must be between 0 and 1).
#' @param endpoint character. The endpoint(s) that must meet a performance goal
#'   criterion. The default is `code = "both"`, which means that the endpoint is
#'   based simultaneously on sensitivity and specificity. Alternative options
#'   are to specify `code = "sens"` or `code = "spec"` for sensitivity and
#'   specificity, respectively. If only a single endpoint is selected (e.g.
#'   sensitivity), then the PG and success probability threshold of the other
#'   statistic are set to 1, and ignored for later analysis.
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
#' @param fut TBA
#' @param n_at_looks vector. Sample sizes for each interim look. The final value
#'   (or only value if no interim looks are planned) is the maximum allowable
#'   sample size for the trial.
#' @param n_mc integer. Number of Monte Carlo draws to use for sampling from the
#'   Beta-Binomial distribution.
#' @param n_trials integer. The number of clinical trials to simulate overall,
#'   which will be used to evaluate the operating characteristics.
#' @param ncores integer. The number of cores to use for parallel processing. If
#'   `ncores` is missing, it defaults to the maximum number of cores available
#'   (spare 1).
#'
#' @details
#'
#' @section Simulation routine:
#'
#' TBA
#'
#' @section Parallelization:
#' To use will multiple cores (where available), the argument \code{ncores}
#'   can be increased from the default of 1. Note: on Windows machines, it is
#'   not possible to use the \code{\link[parallel]{mclapply}} function with
#'   \code{ncores} \eqn{>1}.
#'
#' @return A vector with:
#' - The trial stage (`stage`)
#' - The posterior probability of exceeding the performance goal for
#'   sensitivity (`pp_sens`)
#' - The posterior probability of exceeding the performance goal for
#'   specificity (`pp_spec`)
#' - The true positive count (`tp`)
#' - The true negative count (`tn`)
#' - The false positive count (`fp`)
#' - The false negative count (`fn`)
#'
#' @examples
#'
#' multi_trial(
#'   sens_true = 0.9,
#'   spec_true = 0.95,
#'   prev_true = 0.1,
#'   endpoint = "both",
#'   sens_pg = 0.8,
#'   spec_pg = 0.8,
#'   prior_sens = c(1, 1),
#'   prior_spec = c(1, 1),
#'   prior_prev = c(1, 1),
#'   succ_sens = 0.95,
#'   succ_spec = 0.95,
#'   fut = 0.05,
#'   n_at_looks = c(200, 400, 600, 800, 1000),
#'   min_pos = 1,
#'   n_mc = 10000,
#'   n_trials = 2,
#'   ncores = 1
#' )
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#'
#' @export
multi_trial <- function(
  sens_true,
  spec_true,
  prev_true,
  endpoint = "both",
  sens_pg = 0.8,
  spec_pg = 0.8,
  prior_sens = c(0.1, 0.1),
  prior_spec = c(0.1, 0.1),
  prior_prev = c(0.1, 0.1),
  succ_sens = 0.95,
  succ_spec = 0.95,
  fut = 0.05,
  n_at_looks,
  min_pos = 1,
  n_mc = 10000,
  n_trials = 1000,
  ncores
) {

  Call <- match.call()

  # Check: missing 'ncores' defaults to maximum available (spare 1)
  if (missing(ncores)) {
    ncores <- max(1, parallel::detectCores() - 1)
  }

  # Check: cannot specify <1 core
  if (ncores < 1) {
    warning("Must use at least 1 core... setting ncores = 1")
  }

  # Check: if Windows and if ncores = 1
  if (.Platform$OS.type == "Windows" & ncores > 1L) {
    message("On Windows machines it is required that ncores = 1L")
    ncores <- 1
  }

  # Check: endpoint selection
  if (endpoint == "both") {
    # Both
    if (is.null(sens_pg) | is.null(spec_pg) |
        missing(sens_pg) | missing(spec_pg)) {
      stop("Missing performance goal argument")
    }
    if (is.null(succ_sens) | is.null(succ_spec) |
        missing(succ_sens) | missing(succ_spec)) {
      stop("Missing probability threshold argument")
    }
  } else if (endpoint == "sens") {
    # Sensitivity only
    if (is.null(sens_pg) | missing(sens_pg) | is.na(sens_pg)) {
      stop("Missing performance goal argument")
    }
    if (!is.null(spec_pg)) {
      warning("spec_pg is being ignored")
      spec_pg <- 1 # can never exceed this
      succ_spec <- 1 # can never exceed this
    }
  } else if (endpoint == "spec") {
    # Specificity only
    if (is.null(spec_pg) | missing(spec_pg) | is.na(spec_pg)) {
      stop("Missing performance goal argument")
    }
    if (!is.null(sens_pg)) {
      warning("sens_pg is being ignored")
      sens_pg <- 1 # can never exceed this
      succ_sens <- 1 # can never exceed this
    }
  } else {
    stop("endpoint should be either 'both', 'sens', or 'spec'")
  }

  # Check: futility bound
  if (is.null(fut) | missing(fut) | is.na(fut)) {
    fut <- 1 # never stop for futility
  }

  # Check: true values specified
  if (missing(sens_true) | missing(spec_true) | missing(prev_true)) {
    stop("True values must be provided for for sensitivity, specificity, and prevalence")
  }

  # Check: prior distributions specified
  if (missing(prior_sens) | missing(prior_spec) | missing(prior_prev) |
      is.null(prior_sens) | is.null(prior_spec) | is.null(prior_prev)) {
    stop("Prior distribution parameters must be provided for sensitivity, specificity, and prevalence")
  }

  single_trial_wrapper <- function(x) {
    single_trial(
      sens_true  = sens_true,
      spec_true  = spec_true,
      prev_true  = prev_true,
      sens_pg    = sens_pg,
      spec_pg    = spec_pg,
      prior_sens = prior_sens,
      prior_spec = prior_spec,
      prior_prev = prior_prev,
      succ_sens  = succ_sens,
      succ_spec  = succ_spec,
      n_at_looks = n_at_looks,
      n_mc       = n_mc)
  }

  sims <- pbmclapply(1:n_trials, single_trial_wrapper, mc.cores = ncores)

  sims <- do.call("rbind", sims)
  sims$trial <- rep(1:n_trials, each = length(n_at_looks))
  out <- list(sims = sims, call = Call)

  invisible(out)

}
