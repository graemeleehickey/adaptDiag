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
#'   criterion. The default is \code{code = "both"}, which means that the
#'   endpoint is based simultaneously on sensitivity and specificity.
#'   Alternative options are to specify \code{code = "sens"} or \code{code =
#'   "spec"} for sensitivity and specificity, respectively. If only a single
#'   endpoint is selected (e.g. sensitivity), then the PG and success
#'   probability threshold of the other statistic are set to 1, and ignored for
#'   later analysis.
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
#' This function simulates multiple trials and analyses each stage of the trial
#' (i.e. at each interim analysis sample size look) irrespective of whether a
#' stopping rule was triggered or not. The operating characteristics are handled
#' by a separate function, which accounts for the stopping rules and any other
#' trial constraints. By enumerating each stage of the trial, additional
#' insights can be gained such as: for a trial that stopped early for futility,
#' what is the probability that it would eventually go on to be successful if
#' the trial had not stopped. The details on how each trial are simulated here
#' are described below.
#'
#' \strong{Simulating a single trial}
#'
#' Given true values for the test sensitivity (\code{sens_true}), specificity
#' (\code{spec_true}), and the prevalence (\code{prev_true}) of disease, along
#' with a sample size look strategy (\code{n_at_looks}), it is straightforward
#' to simulate a complete dataset using the binomial distribution. That is, a
#' data frame with true disease status (reference test), and the new diagnostic
#' test result.
#'
#' \strong{Posterior probability of exceeding PG at current look}
#'
#' At a given sample size look, the posterior probability of an endpoint (e.g.
#' sensitivity) exceeding the pre-specified PG (\code{sens_pg}) can be
#' calculated as follows.
#'
#' If we let \eqn{\theta} be the test property of interest (e.g. sensitivity),
#' and if we assume a prior distribution of the form
#'
#' \deqn{\theta ~ Beta(\alpha, \beta),}
#'
#' then with \eqn{X | \theta \sim Bin(n, \theta)}, where \eqn{X} is the number
#' of new test positive cases from the reference positive cases, the posterior
#' distribution of \eqn{\theta} is
#'
#' \deqn{\theta | X=x ~ Beta(\alpha + x, \beta + n - x).}
#'
#' The posterior probability of exceeding the PG is then calculated as
#'
#' \eqn{P[\theta \ge sens_pg | X = x, n]}.
#'
#' A similar calculation can be performed for the specificity, with
#' corresponding PG, \code{spec_pg}.
#'
#' \strong{Posterior predictive probability of eventual success}
#'
#' When at an interim sample size that is less the maximum
#' (i.e. \code{max(n_at_looks)}), we can calculate the probability that the trial
#' will go on to eventually meet the success criteria.
#'
#' At the \eqn{j}-th look, we have observed \eqn{n_j} tests, with \eqn{n_j^* =
#' n_{max} - n_j} subjects yet to be enrolled for testing. For the \eqn{n_j^*}
#' subjects remaining, we can simulate the number of reference positive results,
#' \eqn{y_j^*}, using the posterior predictive distribution for the prevalence
#' (reference positive tests), which is off the form
#'
#' \deqn{y_j^* | y_j, n_j, n_j^* ~ Beta-Bin(n_j^*, \alpha_0 + y_j, \beta + n_j - y_j),}
#'
#' where \eqn{y_j} is the observed number of reference positive cases.
#' Conditional on the number of subjects with a positive reference test in the
#' remaining sample together with \eqn{n_j^*}, one can simulate the complete 2x2
#' contingency table by using the posterior predictive distributions for
#' sensitivity and specificity, each of which has a Beta-Binomial form.
#' Combining the observed \eqn{n_j} subjects' data with a sample of the
#' \eqn{n_j^*} subjects' data drawn from the predictive distribution, one can
#' then calculate the posterior probability of trial success (exceeding a PG)
#' for a specific endpoint. Repeating this many times and calculating the
#' proportion of probabilities that exceed the probability success threshold
#' yields the probability of eventual trial success at the maximum sample size.
#'
#' As well as calculating the predictive posterior probability of eventual
#' success for sensitivity and specificity, separately, we can also calculate
#' the probability for both endpoints simultaneously.
#'
#' @section Parallelization:
#'
#' To use multiple cores (where available), the argument \code{ncores} can be
#' increased from the default of 1. On UNIX machines (including macOS),
#' parallelization is performed using the \code{\link[parallel]{mclapply}}
#' function with \code{ncores} \eqn{>1}. On Windows machines, parallel
#' processing is implemented via the \code{\link[foreach]{foreach}} function.
#'
#' @return A list containing a data frame with rows for each stage of the trial
#'   (i.e. each sample size look), irrespective of whether the trial meets the
#'   stopping criteria. Multiple trial simulations are stacked longways and
#'   indicated by the `trial` column. The data frame has the following columns:
#'
#' \itemize{
#'     \item{\code{stage}:} Trial stage.
#'     \item{\code{pp_sens}:} Posterior probability of exceeding the performance
#'     goal for sensitivity.
#'     \item{\code{pp_spec}:} Posterior probability of exceeding the performance
#'     goal for specificity.
#'     \item{\code{ppp_succ_sens}:} Posterior predictive probability of eventual
#'     success for sensitivity at the maximum sample size.
#'     \item{\code{ppp_succ_spec}:} Posterior predictive probability of eventual
#'     success for specificity at the maximum sample size.
#'     \item{\code{ppp_succ_both}:} Posterior predictive probability of eventual
#'     success for *both* sensitivity and specificity at the maximum sample
#'     size.
#'     \item{\code{tp}:} True positive count.
#'     \item{\code{tn}:} True negative count.
#'     \item{\code{fp}:} False positive count.
#'     \item{\code{fn}:} False negative count.
#'     \item{\code{sens_hat}:} Posterior median estimate of the test
#'     sensitivity.
#'     \item{\code{sens_CrI2.5}:} Lower bound of the 95% credible interval of
#'     the test sensitivity.
#'     \item{\code{sens_CrI97.5}:} Upper bound of the 95% credible interval of
#'     the test sensitivity.
#'     \item{\code{spec_hat}:} Posterior median estimate of the test
#'     specificity.
#'     \item{\code{spec_CrI2.5}:} Lower bound of the 95% credible interval of
#'     the test specificity.
#'     \item{\code{spec_CrI97.5}:} Upper bound of the 95% credible interval of
#'     the test specificity.
#'     \item{\code{n}:} The sample size at the given look for the row.
#'     \item{\code{trial}:} The trial number, which will range from 1 to
#'     `n_trials`.
#' }
#'
#' The list also contains the arguments used and the call.
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
#'   prior_sens = c(0.1, 0.1),
#'   prior_spec = c(0.1, 0.1),
#'   prior_prev = c(0.1, 0.1),
#'   succ_sens = 0.95,
#'   succ_spec = 0.95,
#'   n_at_looks = c(200, 400, 600, 800, 1000),
#'   n_mc = 10000,
#'   n_trials = 2,
#'   ncores = 1
#' )
#'
#' @importFrom parallel detectCores
#' @importFrom pbmcapply pbmclapply
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach registerDoSEQ '%dopar%'
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
  n_at_looks,
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
    }
    spec_pg <- 1   # can never exceed this
    succ_spec <- 1 # can never exceed this
  } else if (endpoint == "spec") {
    # Specificity only
    if (is.null(spec_pg) | missing(spec_pg) | is.na(spec_pg)) {
      stop("Missing performance goal argument")
    }
    if (!is.null(sens_pg)) {
      warning("sens_pg is being ignored")
    }
    sens_pg <- 1   # can never exceed this
    succ_sens <- 1 # can never exceed this
  } else {
    stop("endpoint should be either 'both', 'sens', or 'spec'")
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

  if (.Platform$OS.type == "Windows") {
    # Windows systems
    doParallel::registerDoParallel(cores = ncores)
    sims <- foreach(x = 1:n_trials, .packages = 'adaptDiag',
                    .combine = rbind) %dopar% {
      single_trial_wrapper()
    }
    registerDoSEQ()
  } else {
    # *nix systems
    sims <- pbmclapply(X = 1:n_trials,
                       FUN = single_trial_wrapper,
                       mc.cores = ncores)

    sims <- do.call("rbind", sims)
  }

  sims$trial <- rep(1:n_trials, each = length(n_at_looks))

  args <- list("sens_true"  = sens_true,
               "spec_true"  = spec_true,
               "prev_true"  = prev_true,
               "endpoint"   = endpoint,
               "sens_pg"    = sens_pg,
               "spec_pg"    = spec_pg,
               "prior_sens" = prior_sens,
               "prior_spec" = prior_spec,
               "prior_prev" = prior_prev,
               "succ_sens"  = succ_sens,
               "succ_spec"  = succ_spec,
               "n_at_looks" = n_at_looks,
               "n_mc"       = n_mc,
               "n_trials"   = n_trials)

  out <- list(sims = sims,
              call = Call,
              args = args)

  invisible(out)

}
