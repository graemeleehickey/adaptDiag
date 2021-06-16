##' @title Calculate the minimum number of samples required for a one-sided
##'   exact binomial test
##'
##' @description Calculate the minimum number of samples required for a
##'   one-sided exact binomial test to distinguish between two success
##'   probabilities with specified alpha and power.
##'
##' @param alpha scalar. The desired false positive rate (probability of
##'   incorrectly rejecting the null). Must be be between 0 and 1. Default value
##'   is \code{alpha = 0.05}.
##' @param power scalar. The the minimum probability of correctly rejects the
##'   null when the alternate is true.
##' @param p0 scalar. The expected proportion of successes under the null.
##' @param p1 scalar. The proportion of successes under the alternate
##'   hypothesis.
##'
##' @details This is a one-sided function, such that \eqn{p_0 < p_1}. It
##'   determines the minimum sample size to evaluate the hypothesis test:
##'
##' \deqn{H_0: \, p_1 \le p_0, \, vs.}
##' \deqn{H_1: \, p_1 > p_0}
##'
##' @return A list containing the required sample size and the number of
##'   successful trials required.
##'
##' @references
##' Chow S-C, Shao J, Wang H, Lokhnygina Y. (2017) \emph{Sample Size
##' Calculations in Clinical Research}, Boca Raton, FL: CRC Press.
##'
##' @examples
##'
##' # The minimum number of reference positive cases required to demonstrate
##' # the true sensitivity is >0.7, assuming that the true value is 0.824, with
##' # 90% power is
##'
##' binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.7, p1 = 0.824)
##'
##' # With a sample size of n = 104, if the true prevalence is 0.2, we would
##' # require a sample size of at least n = 520 randomly sampled subjects to
##' # have adequate power to demonstrate the sensitivity of the new test.
##'
##' # The minimum number of reference negative cases required to demonstrate
##' # the true specificity is >0.9, assuming that the true value is 0.963, with
##' # 90% power is
##'
##' binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.9, p1 = 0.963)
##'
##' # The proposed total sample size of n = 520 would be sufficient to
##' # demonstrate both endpoint goals are met.
##'
##' @importFrom stats qnorm qbinom pbinom
##'
##' @export
binom_sample_size <- function(
  alpha = 0.05,
  power = 0.9,
  p0 = 0.9,
  p1 = 0.95) {

  if (p0 >= p1) {
    stop("p0 must be less than p1")
  }

  # Initial estimate: normal approximation
  # Chow et al. (2017), page 85
  N_approx <- (qnorm(power) + qnorm(1 - alpha))^2 * (p1 * (1 - p1)) / (p1 - p0)^2

  # Search range
  N_start <- floor(0.5 * N_approx)
  N_stop <- ceiling(4 * N_approx)
  N <- N_start:N_stop

  # Required number of events at each sample size under null
  # hypothesis (critical value)
  crit.val <- qbinom(p = 1 - alpha, size = N, prob = p0)

  # Calculate beta (type II error) for each N under alternate
  # hypothesis; 1 - beta is the power
  Power <- 1 - pbinom(crit.val, N, p1)

  # Find the smallest sample size yielding at least the required power
  samp.size <- min(which(Power > power))

  # Get the required number of events to reject the null
  # given the sample size required
  return(list(successes = crit.val[samp.size] + 1,
              N = N[samp.size]))

}
