test_that("a tiny legacy mdsbm ICL example runs", {
  x <- matrix(c(0L, 1L, 1L, 0L), nrow = 2, byrow = TRUE)
  dimnames(x) <- list(c("a", "b"), c("a", "b"))
  res <- mdsbm_icl_one_partition(
    x = x,
    sets = 2L,
    k = 1L,
    clu = list(c(0L, 0L))
  )

  expect_type(res, "list")
  expect_true("ICL" %in% names(res))
  expect_true(is.finite(res$ICL))
})
