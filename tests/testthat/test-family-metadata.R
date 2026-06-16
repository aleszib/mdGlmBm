test_that("family metadata normalizes binomial and PPML choices", {
  binom <- glm_blockmodel_family("binomial")
  expect_s3_class(binom, "glm_blockmodel_family")
  expect_identical(binom$name, "binomial")
  expect_false(binom$pseudo)
  expect_identical(binom$objective_scale, "logLik")
  expect_identical(binom$local_score_scale, "deviance")
  expect_true(binom$larger_is_better)
  expect_match(binom$criterion_note, "ordinary likelihood-based")

  ppml_family <- glm_blockmodel_family("ppml")
  expect_s3_class(ppml_family, "glm_blockmodel_family")
  expect_identical(ppml_family$name, "ppml")
  expect_true(ppml_family$pseudo)
  expect_identical(ppml_family$objective_scale, "logLik")
  expect_identical(ppml_family$local_score_scale, "deviance")
  expect_true(ppml_family$larger_is_better)
  expect_match(ppml_family$criterion_note, "pseudo-likelihood-based")
})

test_that("family metadata rejects unsupported names", {
  expect_error(
    glm_blockmodel_family("gaussian"),
    "Supported families are 'binomial' and 'ppml'"
  )
})

test_that("result normalization preserves common criterion names", {
  M_bin <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  dimnames(M_bin) <- list(c("a", "b"), c("a", "b"))
  res_bin <- optParGlm(M_bin, clu = c(1L, 1L), maxIter = 0, glmFamily = binomial)
  norm_bin <- dynGLMbm:::normalize_glm_blockmodel_result(
    list(fit = res_bin$fit, clu = res_bin$clu, ICL = res_bin$ICL),
    family = glm_blockmodel_family("binomial")
  )

  expect_false(norm_bin$pseudo)
  expect_identical(norm_bin$family, "binomial")
  expect_match(norm_bin$criterion_note, "ordinary likelihood-based")
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo", "family", "criterion_note") %in% names(norm_bin)))

  M_ppml <- matrix(c(0, 1, 2, 0), nrow = 2, byrow = TRUE)
  dimnames(M_ppml) <- list(c("a", "b"), c("a", "b"))
  res_ppml <- optParGlm(M_ppml, clu = c(1L, 1L), maxIter = 0, glmFamily = ppml())
  norm_ppml <- dynGLMbm:::normalize_glm_blockmodel_result(
    list(fit = res_ppml$fit, clu = res_ppml$clu, ICL = res_ppml$ICL),
    family = glm_blockmodel_family("ppml")
  )

  expect_true(norm_ppml$pseudo)
  expect_identical(norm_ppml$family, "ppml")
  expect_match(norm_ppml$criterion_note, "pseudo-likelihood-based")
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo", "family", "criterion_note") %in% names(norm_ppml)))
  expect_false(any(c("pseudoICL", "pseudoBIC", "pseudoLogLik") %in% names(norm_ppml)))
})

test_that("static results expose pseudo flags and notes consistently", {
  M_bin <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  dimnames(M_bin) <- list(c("a", "b"), c("a", "b"))
  res_bin <- optParGlm(M_bin, clu = c(1L, 1L), maxIter = 0, glmFamily = binomial)
  expect_false(res_bin$pseudo)
  expect_identical(res_bin$family, "binomial")
  expect_match(res_bin$criterion_note, "ordinary likelihood-based")

  M_ppml <- matrix(c(0, 1, 2, 0), nrow = 2, byrow = TRUE)
  dimnames(M_ppml) <- list(c("a", "b"), c("a", "b"))
  res_ppml <- optParGlm(M_ppml, clu = c(1L, 1L), maxIter = 0, glmFamily = ppml())
  expect_true(res_ppml$pseudo)
  expect_identical(res_ppml$family, "ppml")
  expect_match(res_ppml$criterion_note, "pseudo-likelihood-based")
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective") %in% names(res_ppml)))
})
