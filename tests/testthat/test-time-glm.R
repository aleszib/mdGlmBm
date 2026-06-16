test_that("fit_time_glm_blockmodel works for one time point with binomial family", {
  Y <- list(
    t1 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dn <- as_dynamic_network(Y)

  fit <- fit_time_glm_blockmodel(dn, membership = c(1L, 1L), family = "binomial")

  expect_s3_class(fit, "time_glm_blockmodel")
  expect_false(fit$pseudo)
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo") %in% names(fit)))
  expect_identical(fit$family, "binomial")
  expect_identical(fit$time, "t1")
  expect_true(is.na(fit$ICL))
})

test_that("fit_time_glm_blockmodels works across multiple times with actor entry and exit", {
  Y <- list(
    t1 = matrix(c(0, 1, 0, 0, 1, 0, 0, 1, 0), nrow = 3, byrow = TRUE),
    t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE),
    t3 = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dimnames(Y$t2) <- list(c("A", "C", "D"), c("A", "C", "D"))
  dimnames(Y$t3) <- list(c("C", "D", "E"), c("C", "D", "E"))
  dn <- as_dynamic_network(Y)

  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 1L, 2L, 1L, 2L, 2L, 2L, 3L, 3L),
    stringsAsFactors = FALSE
  )

  fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")

  expect_s3_class(fit, "time_glm_blockmodels")
  expect_false(fit$pseudo)
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo") %in% names(fit)))
  expect_identical(names(fit$fits), dn$times)
  expect_equal(length(fit$fits), 3L)
  expect_equal(nrow(fit$membership), nrow(dn$actor_time))
  expect_true(is.finite(fit$logLik_total))
  expect_true(is.finite(fit$deviance_total))
  expect_true(all(is.na(fit$ICL)))
})

test_that("PPML time-specific fits set pseudo = TRUE", {
  Y <- list(
    t1 = matrix(c(0, 1,
                  2, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dn <- as_dynamic_network(Y)

  fit <- fit_time_glm_blockmodel(dn, membership = c(1L, 1L), family = "ppml")

  expect_true(fit$pseudo)
  expect_identical(fit$family, "ppml")
  expect_match(fit$criterion_note, "pseudo-likelihood-based")
})

test_that("returned time-specific fits keep common result names", {
  Y <- list(
    t1 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE),
    t2 = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dimnames(Y$t2) <- list(c("A", "B"), c("A", "B"))
  dn <- as_dynamic_network(Y)

  membership <- c(1L, 1L, 1L, 1L)
  fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")

  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo") %in% names(fit)))
  expect_identical(fit$objective, fit$logLik_total)
  expect_identical(fit$family, "binomial")
})

test_that("membership validation fails on missing or misaligned memberships", {
  Y <- list(
    t1 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dn <- as_dynamic_network(Y)

  expect_error(
    fit_time_glm_blockmodel(dn, membership = c(1L), family = "binomial"),
    "one value per actor-time unit"
  )

  expect_error(
    fit_time_glm_blockmodel(
      dn,
      membership = setNames(c(1L, 1L), c("x", "y")),
      family = "binomial"
    ),
    "do not align"
  )
})

test_that("one-based cluster labels are accepted for time-specific fits", {
  Y <- list(
    t1 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dn <- as_dynamic_network(Y)

  fit <- fit_time_glm_blockmodel(dn, membership = c(1L, 2L), family = "binomial")
  expect_identical(fit$membership$membership, c(1L, 2L))
})
