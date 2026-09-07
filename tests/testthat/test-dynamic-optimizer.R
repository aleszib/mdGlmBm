make_dynamic_optimizer_fixture <- function() {
  Y <- list(
    t1 = matrix(c(0, 1, 0,
                  1, 0, 1,
                  0, 1, 0), nrow = 3, byrow = TRUE),
    t2 = matrix(c(0, 1, 1,
                  1, 0, 0,
                  1, 0, 0), nrow = 3, byrow = TRUE),
    t3 = matrix(c(0, 0, 1,
                  0, 0, 1,
                  1, 1, 0), nrow = 3, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dimnames(Y$t2) <- list(c("A", "C", "D"), c("A", "C", "D"))
  dimnames(Y$t3) <- list(c("C", "D", "E"), c("C", "D", "E"))

  dn <- as_dynamic_network(Y)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L),
    stringsAsFactors = FALSE
  )
  list(network = dn, membership = membership)
}

test_that("fit_dynamic_glm_blockmodel runs on a tiny binomial dynamic network", {
  fx <- make_dynamic_optimizer_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network,
    membership = fx$membership,
    k = 2,
    family = "binomial",
    max_iter = 0
  )

  expect_s3_class(fit, "dynamic_glm_blockmodel")
  expect_false(fit$pseudo)
  expect_equal(nrow(fit$membership), nrow(fx$network$actor_time))
  expect_equal(fit$initial_membership$membership, fx$membership$membership)
  expect_equal(fit$membership$membership, fx$membership$membership)
  expect_true(all(c("objective", "objective_history", "history", "transition", "prior") %in% names(fit)))
  expect_equal(length(fit$objective_history), 1L)
  expect_equal(fit$n_iter, 0L)
  expect_false(fit$converged)
  expect_true(is.finite(fit$objective))
  expect_true(is.finite(fit$deviance_total))
  expect_true(is.finite(fit$transition_penalty_total))
  expect_true(is.finite(fit$prior_penalty_total))
})

test_that("fit_dynamic_glm_blockmodel stops when no memberships change", {
  fx <- make_dynamic_optimizer_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network,
    membership = data.frame(
      unit_id = fx$membership$unit_id,
      membership = 1L,
      stringsAsFactors = FALSE
    ),
    k = 1,
    family = "binomial",
    max_iter = 5,
    prior = "none"
  )

  expect_s3_class(fit, "dynamic_glm_blockmodel")
  expect_true(fit$converged)
  expect_equal(fit$n_changes, 0L)
  expect_equal(fit$n_iter, 1L)
  expect_equal(length(fit$objective_history), 2L)
  expect_equal(unique(fit$membership$membership), 1L)
})

test_that("fit_dynamic_glm_blockmodel supports prior modes and entry-exit penalties", {
  fx <- make_dynamic_optimizer_fixture()
  for (prior_mode in c("uniform", "empirical", "none")) {
    fit <- fit_dynamic_glm_blockmodel(
      fx$network,
      membership = data.frame(
        unit_id = fx$membership$unit_id,
        membership = 1L,
        stringsAsFactors = FALSE
      ),
      k = 1,
      family = "binomial",
      max_iter = 1,
      prior = prior_mode
    )

    expect_identical(fit$prior$prior, prior_mode)
    expect_true(is.finite(fit$transition_penalty_total))
    expect_true(is.finite(fit$objective))
    expect_true(all(fit$membership$membership == 1L))
  }
})

test_that("fit_dynamic_glm_blockmodel accepts PPML mode on a stable tiny example", {
  Y <- list(
    t1 = matrix(c(0, 1, 2,
                  1, 0, 1,
                  2, 1, 0), nrow = 3, byrow = TRUE),
    t2 = matrix(c(0, 2, 1,
                  1, 0, 2,
                  0, 1, 0), nrow = 3, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dn <- as_dynamic_network(Y)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = 1L,
    stringsAsFactors = FALSE
  )

  fit <- fit_dynamic_glm_blockmodel(
    dn,
    membership = membership,
    k = 1,
    family = "ppml",
    max_iter = 1,
    prior = "none"
  )

  expect_s3_class(fit, "dynamic_glm_blockmodel")
  expect_true(fit$pseudo)
  expect_identical(fit$family, "ppml")
  expect_true(is.finite(fit$objective))
})

test_that("fit_dynamic_glm_blockmodel rejects missing or misaligned initial memberships", {
  fx <- make_dynamic_optimizer_fixture()

  expect_error(
    fit_dynamic_glm_blockmodel(
      fx$network,
      membership = c(1L, 1L),
      k = 1,
      family = "binomial"
    ),
    "one value per actor-time unit"
  )

  expect_error(
    fit_dynamic_glm_blockmodel(
      fx$network,
      membership = data.frame(
        unit_id = paste0("u", seq_len(nrow(fx$membership))),
        membership = fx$membership$membership,
        stringsAsFactors = FALSE
      ),
      k = 2,
      family = "binomial"
    ),
    "missing rows for one or more actor-time units"
  )
})

test_that("one-based cluster labels are preserved in the dynamic optimizer", {
  fx <- make_dynamic_optimizer_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network,
    membership = fx$membership,
    k = 2,
    family = "binomial",
    max_iter = 0
  )

  expect_true(all(fit$membership$membership >= 1L))
  expect_setequal(unique(fit$membership$membership), unique(fx$membership$membership))
})

test_that("random initialization creates complete one-based partitions", {
  fx <- make_dynamic_optimizer_fixture()
  init <- initialize_dynamic_membership(fx$network, k = 2L)

  expect_equal(nrow(init), nrow(fx$network$actor_time))
  expect_false(anyNA(init$membership))
  expect_true(all(init$membership %in% 1:2))
  expect_setequal(unique(init$membership), 1:2)
  expect_equal(init$unit_id, fx$network$actor_time$unit_id)
})

test_that("missing membership requires fixed k and supports actor entry and exit", {
  fx <- make_dynamic_optimizer_fixture()

  expect_error(
    fit_dynamic_glm_blockmodel(fx$network, max_iter = 0),
    "`k` must be supplied"
  )

  fit <- fit_dynamic_glm_blockmodel(
    fx$network,
    membership = NULL,
    k = 2L,
    n_starts = 2L,
    max_iter = 0,
    seed = 123
  )
  expect_equal(nrow(fit$membership), nrow(fx$network$actor_time))
  expect_true(all(fit$membership$membership %in% 1:2))
  expect_length(fit$starts, 2L)
})

test_that("explicit seeds reproduce starts and preserve the caller RNG state", {
  fx <- make_dynamic_optimizer_fixture()
  set.seed(456)
  state_before <- .Random.seed

  fit1 <- fit_dynamic_glm_blockmodel(
    fx$network, k = 2L, n_starts = 3L, max_iter = 0, seed = 123
  )
  state_after <- .Random.seed
  fit2 <- fit_dynamic_glm_blockmodel(
    fx$network, k = 2L, n_starts = 3L, max_iter = 0, seed = 123
  )

  expect_equal(state_after, state_before)
  expect_identical(fit1$seed, 123)
  expect_equal(
    lapply(fit1$starts, `[[`, "initial_membership"),
    lapply(fit2$starts, `[[`, "initial_membership")
  )
  expect_equal(fit1$objective, fit2$objective)
  expect_equal(fit1$membership, fit2$membership)
})

test_that("seed NULL uses the current RNG stream rather than a fixed seed", {
  fx <- make_dynamic_optimizer_fixture()
  set.seed(789)
  expected <- sample(rep(1:2, length.out = nrow(fx$network$actor_time)))
  set.seed(789)
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, k = 2L, n_starts = 1L, max_iter = 0, seed = NULL
  )

  expect_null(fit$seed)
  expect_equal(fit$starts[[1]]$initial_membership$membership, expected)
})

test_that("supplied membership is retained as the first of multiple starts", {
  fx <- make_dynamic_optimizer_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network,
    membership = fx$membership,
    k = 2L,
    n_starts = 3L,
    max_iter = 0,
    seed = 123
  )

  expect_equal(fit$starts[[1]]$initial_membership$membership, fx$membership$membership)
  expect_equal(fit$n_starts, 3L)
})

test_that("multiple starts select the minimum objective and retain diagnostics", {
  fx <- make_dynamic_optimizer_fixture()
  fit <- fit_dynamic_glm_blockmodel(
    fx$network, k = 2L, n_starts = 3L, max_iter = 0, seed = 123
  )

  expect_equal(fit$objective, min(fit$start_objectives))
  expect_equal(fit$best_start, which.min(fit$start_objectives))
  expect_length(fit$start_objectives, 3L)
  expect_length(fit$start_converged, 3L)
  expect_true(all(vapply(fit$starts, function(x) {
    all(c("initial_membership", "final_membership", "objective",
          "converged", "n_iter") %in% names(x))
  }, logical(1))))
})

test_that("initialization and multiple-start arguments are validated", {
  fx <- make_dynamic_optimizer_fixture()

  expect_error(fit_dynamic_glm_blockmodel(fx$network, k = 2L, n_starts = 0),
               "positive integer")
  expect_error(fit_dynamic_glm_blockmodel(fx$network, k = 2L, n_starts = -1),
               "positive integer")
  expect_error(fit_dynamic_glm_blockmodel(fx$network, k = 2L, n_starts = 1.5),
               "positive integer")
  expect_error(initialize_dynamic_membership(fx$network, k = 20L),
               "every cluster must be represented")
  expect_error(
    fit_dynamic_glm_blockmodel(fx$network, membership = fx$membership, k = 1L, max_iter = 0),
    "cannot be smaller"
  )
})
