test_that("match broglio", {

  skip_on_cran()
  skip_on_ci()
  skip_on_covr()
  skip()

  library(VGAM)

  design.parameters <- list(sens.goal = 0.7,
                            spec.goal = 0.9,
                            looks = seq(200, 700, by = 50),
                            sens.alpha = 0.1,
                            sens.beta = 0.1,
                            spec.alpha = 0.1,
                            spec.beta = 0.1,
                            ref.alpha = 0.1,
                            ref.beta = 0.1,
                            sens.critv = 0.985,
                            spec.critv = 0.985,
                            min.hist.pos = 35,
                            fut.bound = 0.05)

  needed.to.win.sens <- rep(NA, 700)
  needed.to.win.spec <- rep(NA, 700)

  min.right.sens <- function(N, design.parameters) {
    x <- 0
    post <- 0
    while (post < design.parameters$sens.critv) {
      x <- x + 1
      post <- 1 - pbeta(design.parameters$sens.goal, 0.1 + x, 0.1 + N - x)
    }
    return(x)
  }

  min.right.spec <- function(N, design.parameters) {
    x <- 0
    post <- 0
    while (post < design.parameters$spec.critv) {
      x <- x + 1
      post <- 1 - pbeta(design.parameters$spec.goal, 0.1 + x, 0.1 + N - x)
    }
    return(x)
  }

  for (n in 35:700) {
    needed.to.win.sens[n] <- min.right.sens(n, design.parameters)
    needed.to.win.spec[n] <- min.right.spec(n, design.parameters)
  }

  make.data <- function(N, true.hist.pos, true.sens, true.spec) {
    hist.pos <- rbinom(N, 1, true.hist.pos)
    test.pos <- rbinom(N, 1, ifelse(hist.pos == 1, true.sens, 1 - true.spec))
    data <- cbind(hist.pos, test.pos)
    return(data)
  }

  primary.analysis <- function(N, data, design.parameters, pr) {
    pm <- design.parameters
    tab2x2 <- table(data[1:N, 1], data[1:N, 2])
    PostProbSens <- 1 - pbeta(pm$sens.goal, pm$sens.alpha + tab2x2[2, 2], pm$sens.beta + tab2x2[2, 1])
    PostProbSpec <- 1 - pbeta(pm$spec.goal, pm$spec.alpha + tab2x2[1, 1], pm$spec.beta + tab2x2[1, 2])
    N.left <- max(pm$looks) - N
    if (N.left > 0) {
      x.left.hist.pos <- rbetabinom.ab(10000, N.left, pm$ref.alpha + tab2x2[2, 1] + tab2x2[2, 2], pm$ref.beta + tab2x2[1, 1] + tab2x2[1, 2])
      x.left.hist.neg <- N.left - x.left.hist.pos
      x.left.hist.pos[x.left.hist.pos == 0] <- 1
      x.left.hist.pos.test.pos <- rbetabinom.ab(10000, x.left.hist.pos, pm$sens.alpha + tab2x2[2, 2], pm$sens.beta + tab2x2[2, 1])
      x.left.hist.pos.test.pos[x.left.hist.pos == 0] <- 0
      x.left.hist.neg.test.neg <- rbetabinom.ab(10000, x.left.hist.neg, pm$spec.alpha + tab2x2[1, 1], pm$spec.beta + tab2x2[1, 2])
      hist.pos.at.max <- tab2x2[2, 1] + tab2x2[2, 2] + x.left.hist.pos
      hist.pos.test.pos.at.max <- tab2x2[2, 2] + x.left.hist.pos.test.pos
      win.sens.at.max <- (hist.pos.test.pos.at.max >= needed.to.win.sens[hist.pos.at.max])
      hist.neg.at.max <- tab2x2[1, 1] + tab2x2[1, 2] + x.left.hist.neg
      hist.neg.test.neg.at.max <- tab2x2[1, 1] + x.left.hist.neg.test.neg
      win.spec.at.max <- (hist.neg.test.neg.at.max >= needed.to.win.spec[hist.neg.at.max])
      pred.prob.both.goals.at.max <- mean(win.sens.at.max * win.spec.at.max)
    } else {
      pred.prob.both.goals.at.max <- NA
    }
    if (pr) {
      print(c(N, c(tab2x2), PostProbSens, PostProbSpec, pred.prob.both.goals.at.max))
    }
    return(c(N,
             tab2x2[2, 1] + tab2x2[2, 2],
             PostProbSens,
             PostProbSpec,
             pred.prob.both.goals.at.max))
  }

  stop.check <- function(interim.analysis,
                         design.parameters,
                         pr) {
    go <- 1
    if (interim.analysis[2] >= design.parameters$min.hist.pos &
        interim.analysis[3] >= design.parameters$sens.critv &
        interim.analysis[4] >= design.parameters$spec.critv
    ) {
      go <- 0
      win <- 1
      stop <- 3
      n.trial <- interim.analysis[1]
    } else if (interim.analysis[1] >= max(design.parameters$looks)) {
      go <- 0
      win <- 0
      stop <- 2
      n.trial <- interim.analysis[1]
    } else if (interim.analysis[2] >= design.parameters$min.hist.pos &
               interim.analysis[5] <= design.parameters$fut.bound) {
      go <- 0
      win <- 0
      stop <- 1
      n.trial <- interim.analysis[1]
    } else {
      go <- 1
      win <- NA
      stop <- NA
      n.trial <- NA
    }
    if (pr) {
      print(c(go, win, stop, n.trial))
    }
    return(c(go, win, stop, n.trial))
  }

  sim.trial <- function(design.parameters,
                        true.hist.pos,
                        true.sens,
                        true.spec,
                        pr = FALSE) {
    data <- make.data(max(design.parameters$look),
                      true.hist.pos, true.sens, true.spec)
    go <- 1
    look <- 1
    while (go == 1) {
      n.at.look <- design.parameters$look[look]
      test.statistics <- primary.analysis(n.at.look, data, design.parameters, pr)
      trial.result <- stop.check(test.statistics, design.parameters, pr)
      go <- trial.result[1]
      if (go == 1) {
        look <- look + 1
      }
    }

    N <- trial.result[4]
    return(c(trial.result[2:4],
             mean(data[1:N, 2][data[1:N, 1] == 1]),
             1 - mean(data[1:N, 2][data[1:N, 1] == 0]),
             test.statistics[3:4]))
  }

  sim.N.trials.summarize <- function(
    Nsims,
    design.parameters,
    true.hist.pos,
    true.sens,
    true.spec,
    pr = 10,
    long = rep(TRUE, 5)) {
    results.matrix <- matrix(nrow = Nsims, ncol = 7)
    for (s in 1:Nsims) {
      print(paste("s =", s))
      results.matrix[s, ] <- sim.trial(design.parameters,
                                       true.hist.pos, true.sens,
                                       true.spec, s <= pr)
    }
    t1 <- round(table(factor(results.matrix[, 1],
                             levels = 0:1,
                             labels = c("Fail", "Success"))) / Nsims, 3)
    t2 <- round(table(factor(results.matrix[, 2],
                             levels = 1:3,
                             labels = c("Futility", "MaxN", "EarlySuccess"))) / Nsims, 3)
    t3 <- round(table(factor(results.matrix[, 3], design.parameters$looks),
                      factor(results.matrix[, 1],
                             labels = c("Fail", "Success"))) / Nsims, 3)
    t4 <- round(table(factor(results.matrix[, 3], design.parameters$looks),
                      factor(results.matrix[, 2],
                             levels = 1:3,
                             labels = c("Futility", "MaxN", "EarlySuccess"))) / Nsims, 3)
    t5 <- round(table(factor(as.numeric(results.matrix[, 6] >= design.parameters$sens.critv),
                             levels = 0:1,
                             labels = c("Lose Spec", "Win Spec"))) / Nsims, 3)
    if (long[1]) {
      print(t1)
    }
    if (long[2]) {
      print(t2)
    }
    if (long[3]) {
      print(t3)
    }
    if (long[4]) {
      print(t4)
    }
    if (long[5]) {
      print(t5)
    }
    cat("\n\n")
    return(c(Power = mean(results.matrix[, 1]),
             PrFutStop = mean(results.matrix [, 2] == 1),
             PrMaxN = mean(results.matrix[, 2] == 2),
             PrEarlySuc = mean(results.matrix[, 2] == 3),
             MeanN = mean(results.matrix[, 3]),
             MeanSens = mean(results.matrix[, 4]),
             MeanSpec = mean(results.matrix[, 5])))
  }

  broglio <- sim.N.trials.summarize(1000, design.parameters, 0.1, 0.85, 0.95)

  data <- multi_trial(
    sens_true = 0.85,
    spec_true = 0.95,
    prev_true = 0.1,
    endpoint = "both",
    sens_pg = 0.7,
    spec_pg = 0.9,
    prior_sens = c(0.1, 0.1),
    prior_spec = c(0.1, 0.1),
    prior_prev = c(0.1, 0.1),
    succ_sens = 0.985,
    succ_spec = 0.985,
    n_at_looks = seq(200, 700, by = 50),
    n_mc = 10000,
    n_trials = 1000,
    ncores = 8
  )

  adaptDiag_out <- summarise_trials(data, fut = 0.05, min_pos = 35)

  expect_equal(round(as.numeric(broglio[1]), 2),
               round(adaptDiag_out$power, 2),
               tolerance = 0.01)
  expect_equal(round(as.numeric(broglio[2]), 2),
               round(adaptDiag_out$stop_futility, 2),
               tolerance = 0.01)
  expect_equal(floor(as.numeric(broglio[5])),
               floor(adaptDiag_out$n_avg),
               tolerance = 5)
  expect_equal(round(as.numeric(broglio[6]), 2),
               round(adaptDiag_out$sens, 2),
               tolerance = 0.02)
  expect_equal(round(as.numeric(broglio[7]), 2),
               round(adaptDiag_out$spec, 2),
               tolerance = 0.01)
})


