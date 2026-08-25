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
