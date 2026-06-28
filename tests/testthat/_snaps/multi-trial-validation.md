# multi_trial errors on invalid endpoint

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "invalid",
        n_at_looks = c(100, 200), n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! endpoint should be either 'both', 'sens', or 'spec'

# multi_trial errors when sens_pg missing for 'both' endpoint

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "both",
        sens_pg = NULL, spec_pg = 0.8, succ_sens = 0.95, succ_spec = 0.95,
        n_at_looks = c(100, 200), n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! Missing performance goal argument

# multi_trial errors when succ thresholds missing for 'both' endpoint

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "both",
        sens_pg = 0.8, spec_pg = 0.8, succ_sens = NULL, succ_spec = 0.95, n_at_looks = c(
          100, 200), n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! Missing probability threshold argument

# multi_trial errors when sens_pg missing for 'sens' endpoint

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "sens",
        sens_pg = NULL, n_at_looks = c(100, 200), n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! Missing performance goal argument

# multi_trial errors when spec_pg missing for 'spec' endpoint

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "spec",
        spec_pg = NULL, n_at_looks = c(100, 200), n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! Missing performance goal argument

# multi_trial errors when prior arguments are NULL

    Code
      multi_trial(sens_true = 0.9, spec_true = 0.95, prev_true = 0.1, endpoint = "both",
        sens_pg = 0.8, spec_pg = 0.8, succ_sens = 0.95, succ_spec = 0.95, prior_sens = NULL,
        prior_spec = c(1, 1), prior_prev = c(1, 1), n_at_looks = c(100, 200),
        n_trials = 2, ncores = 1)
    Condition
      Error in `multi_trial()`:
      ! Prior distribution parameters must be provided for sensitivity, specificity, and prevalence

