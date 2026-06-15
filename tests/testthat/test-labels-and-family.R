test_that("label conversion helpers round-trip", {
  expect_identical(dynGLMbm:::bm_zero_index(c(1L, 3L, 5L)), c(0L, 2L, 4L))
  expect_identical(dynGLMbm:::bm_one_index(c(0L, 2L, 4L)), c(1L, 3L, 5L))
  expect_identical(
    dynGLMbm:::bm_one_index(dynGLMbm:::bm_zero_index(c(2L, 4L))),
    c(2L, 4L)
  )
})

test_that("ppml returns a PPML family object", {
  fam <- ppml()
  expect_s3_class(fam, "family")
  expect_identical(fam$family, "ppml")
  expect_true(isTRUE(attr(fam, "pseudo")))

  fit <- glm(
    y ~ x,
    family = fam,
    data = data.frame(y = c(0, 1, 2, 3), x = c(0, 1, 0, 1))
  )
  expect_s3_class(fit, "glm")
})
