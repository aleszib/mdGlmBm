test_that("the ARI helper is label invariant", {
  truth <- c(1L, 1L, 1L, 2L, 2L, 2L)
  expect_equal(.adjusted_rand_index(truth, truth), 1)
  expect_equal(.adjusted_rand_index(truth, 3L - truth), 1)
  expect_lt(.adjusted_rand_index(truth, c(1L, 2L, 1L, 2L, 1L, 2L)), 1)
})

test_that("strong binomial dynamic structure is substantially recovered", {
  fx <- simulate_binomial_validation_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, k = fx$k, n_starts = 5L, seed = 77L,
    max_iter = 5L, prior = "none"
  )
  estimated <- fit$membership$membership
  overall_ari <- .adjusted_rand_index(estimated, fx$truth_membership$membership)
  per_time_ari <- vapply(seq_along(fx$network$times), function(i) {
    idx <- fx$network$actor_time$time_index == i
    .adjusted_rand_index(estimated[idx], fx$truth_membership$membership[idx])
  }, numeric(1))

  expect_gte(overall_ari, 0.8)
  expect_gte(mean(per_time_ari), 0.8)
  expect_true(all(is.finite(fit$objective)))
})

test_that("entry and exit validation remains finite and aligned", {
  fx <- simulate_entry_exit_fixture()
  dn <- fx$network
  fit <- fit_dynamic_glm_blockmodel(
    dn, k = fx$k, n_starts = 3L, seed = 123L,
    max_iter = 2L, prior = "none"
  )
  expect_equal(nrow(fit$membership), nrow(dn$actor_time))
  expect_true(all(is.finite(c(fit$objective, fit$deviance_total,
                              fit$transition_penalty_total,
                              fit$prior_penalty_total))))

  entry <- which(dn$actor_time$actor_id == "I" & dn$actor_time$time == "t2")
  exit <- which(dn$actor_time$actor_id == "B4" & dn$actor_time$time == "t2")
  scores_entry <- score_actor_time_candidates(
    dn, fit$fits, membership = fit$membership, row_index = entry,
    candidate_clusters = 1:2, transition = fit$transition,
    prior = "none"
  )
  scores_exit <- score_actor_time_candidates(
    dn, fit$fits, membership = fit$membership, row_index = exit,
    candidate_clusters = 1:2, transition = fit$transition,
    prior = "none"
  )
  expect_true(all(scores_entry$previous_transition_penalty == 0))
  expect_true(all(scores_exit$next_transition_penalty == 0))
})

test_that("PPML validation preserves pseudo-likelihood metadata", {
  fx <- simulate_ppml_validation_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, membership = fx$truth_membership, k = fx$k,
    family = "ppml", n_starts = 2L, seed = 123L,
    max_iter = 1L, prior = "none"
  )
  expect_true(fit$pseudo)
  expect_identical(fit$family, "ppml")
  expect_true(is.finite(fit$objective))
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective") %in% names(fit)))
})

test_that("multiple-start selection is no worse than its corresponding first start", {
  fx <- simulate_binomial_validation_fixture()
  fit_one <- fit_dynamic_glm_blockmodel(
    fx$network, k = fx$k, n_starts = 1L, seed = 123L,
    max_iter = 2L, prior = "none"
  )
  fit_many <- fit_dynamic_glm_blockmodel(
    fx$network, k = fx$k, n_starts = 5L, seed = 123L,
    max_iter = 2L, prior = "none"
  )
  expect_equal(
    fit_one$starts[[1L]]$initial_membership,
    fit_many$starts[[1L]]$initial_membership
  )
  expect_lte(fit_many$objective, fit_one$objective + 1e-8)
  expect_equal(fit_many$objective, min(fit_many$start_objectives), tolerance = 1e-8)
  expect_equal(fit_many$best_start, which.min(fit_many$start_objectives))
})

test_that("global cluster relabeling preserves the objective", {
  fx <- simulate_binomial_validation_fixture()
  relabeled <- fx$truth_membership
  relabeled$membership <- 3L - relabeled$membership
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, membership = fx$truth_membership, k = 2L,
    max_iter = 0L, prior = "empirical"
  )
  fit_relabel <- fit_dynamic_glm_blockmodel(
    fx$network, membership = relabeled, k = 2L,
    max_iter = 0L, prior = "empirical"
  )
  expect_equal(fit_relabel$objective, fit$objective, tolerance = 1e-8)
  expect_equal(fit_relabel$deviance_total, fit$deviance_total, tolerance = 1e-8)
  expect_equal(fit_relabel$transition_penalty_total,
               fit$transition_penalty_total, tolerance = 1e-8)
})

test_that("consistent actor reordering preserves the mathematical result", {
  fx <- simulate_binomial_validation_fixture()
  order <- rev(seq_len(nrow(fx$network$matrices[[1L]])))
  Y_reordered <- lapply(fx$network$matrices, function(M) M[order, order, drop = FALSE])
  dn_reordered <- as_dynamic_network(Y_reordered, directed = TRUE, self_ties = FALSE,
                                     times = fx$network$times)
  truth_reordered <- data.frame(
    unit_id = dn_reordered$actor_time$unit_id,
    membership = fx$truth_membership$membership[
      match(dn_reordered$actor_time$actor_id, fx$network$actor_time$actor_id)
    ],
    stringsAsFactors = FALSE
  )
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, membership = fx$truth_membership, k = 2L,
    max_iter = 0L, prior = "empirical"
  )
  fit_reordered <- fit_dynamic_glm_blockmodel(
    dn_reordered, membership = truth_reordered, k = 2L,
    max_iter = 0L, prior = "empirical"
  )
  expect_equal(fit_reordered$objective, fit$objective, tolerance = 1e-8)
  expect_equal(fit_reordered$deviance_total, fit$deviance_total, tolerance = 1e-8)
})

test_that("objective and local score components add on the documented scale", {
  fx <- simulate_binomial_validation_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, membership = fx$truth_membership, k = 2L,
    max_iter = 0L, prior = "empirical"
  )
  recomputed <- fit$deviance_total + fit$transition_penalty_total + fit$prior_penalty_total
  expect_equal(fit$objective, recomputed, tolerance = 1e-8)

  row <- which(fx$network$actor_time$time_index == 2L)[1L]
  scores <- score_actor_time_candidates(
    fx$network, fit$fits, membership = fit$membership, row_index = row,
    candidate_clusters = 1:2, transition = fit$transition, prior = fit$prior
  )
  expect_equal(
    scores$total_score,
    scores$glm_deviance + scores$previous_transition_penalty +
      scores$next_transition_penalty + scores$prior_penalty,
    tolerance = 1e-8
  )
})

test_that("optimizer history exposes finite convergence diagnostics", {
  fx <- simulate_binomial_validation_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, k = 2L, n_starts = 2L, seed = 123L,
    max_iter = 2L, prior = "none"
  )
  expect_true(all(c("iteration", "n_changes", "objective",
                    "deviance_total", "transition_penalty_total",
                    "prior_penalty_total") %in% names(fit$history)))
  expect_true(all(is.finite(fit$history$objective)))
  expect_true(fit$stopping_reason %in% c("no_changes", "max_iter"))
  expect_equal(fit$objective, tail(fit$history$objective, 1L), tolerance = 1e-8)
})

test_that("the known truth is a sensible starting diagnostic", {
  fx <- simulate_binomial_validation_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, membership = fx$truth_membership, k = 2L,
    max_iter = 2L, prior = "none"
  )
  truth_ari <- .adjusted_rand_index(fit$membership$membership,
                                    fx$truth_membership$membership)
  expect_gte(truth_ari, 0.8)
  expect_true(is.finite(fit$objective))
})
