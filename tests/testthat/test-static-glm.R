test_that("balassaNorm works on a tiny matrix", {
  m <- matrix(c(10, 5, 4, 8), nrow = 2, byrow = TRUE)
  out <- balassaNorm(m)
  expect_equal(dim(out), c(2L, 2L))
  expect_true(all(is.finite(out)))
})

test_that("optParGlm returns consistent metadata for binomial and PPML", {
  M_bin <- matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
  dimnames(M_bin) <- list(c("a", "b"), c("a", "b"))
  res_bin <- optParGlm(M_bin, clu = c(1L, 1L), maxIter = 0, glmFamily = binomial)

  expect_false(res_bin$pseudo)
  expect_identical(res_bin$family, "binomial")
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo", "membership") %in% names(res_bin)))
  expect_identical(res_bin$objective, res_bin$ICL)
  expect_true(is.finite(res_bin$logLik))

  M_ppml <- matrix(c(0, 1, 2, 0), nrow = 2, byrow = TRUE)
  dimnames(M_ppml) <- list(c("a", "b"), c("a", "b"))
  res_ppml <- optParGlm(M_ppml, clu = c(1L, 1L), maxIter = 0, glmFamily = ppml())

  expect_true(res_ppml$pseudo)
  expect_identical(res_ppml$family, "ppml")
  expect_true(all(c("logLik", "BIC", "ICL", "deviance", "objective", "pseudo", "membership") %in% names(res_ppml)))
  expect_identical(res_ppml$objective, res_ppml$ICL)
  expect_true(is.finite(res_ppml$logLik))
})
