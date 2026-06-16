test_that("as_dynamic_network builds actor-time and lineage tables for constant actor sets", {
  Y <- list(
    t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
    t2 = matrix(c(0, 0, 1, 1), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dimnames(Y$t2) <- list(c("A", "B"), c("A", "B"))

  dn <- as_dynamic_network(Y)

  expect_s3_class(dn, "dynamic_network")
  expect_equal(dn$times, c("t1", "t2"))
  expect_equal(nrow(dn$actor_time), 4L)
  expect_equal(nrow(dn$lineage), 2L)
  expect_equal(dn$lineage$actor_id, c("A", "B"))
  expect_true(all(dn$lineage$relation == "identity"))
  expect_true(all(dn$lineage$weight == 1))
})

test_that("as_dynamic_network accepts actors entering and leaving over time", {
  Y <- list(
    matrix(c(0, 1, 0, 0, 1, 0, 0, 1, 0), nrow = 3, byrow = TRUE),
    matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE),
    matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, byrow = TRUE)
  )
  times <- c("t1", "t2", "t3")
  dimnames(Y[[1]]) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dimnames(Y[[2]]) <- list(c("A", "C", "D"), c("A", "C", "D"))
  dimnames(Y[[3]]) <- list(c("C", "D", "E"), c("C", "D", "E"))

  dn <- as_dynamic_network(Y, times = times)

  expect_equal(dn$times, times)
  expect_equal(nrow(dn$actor_time), 9L)
  expect_equal(dn$n_actors_by_time, c(3L, 3L, 3L))
  expect_equal(dn$lineage$actor_id, c("A", "C", "C", "D"))
  expect_false("B" %in% dn$lineage$actor_id)
  expect_false("E" %in% dn$lineage$actor_id)
})

test_that("as_dynamic_network rejects matrices without dimnames or with mismatched dimnames", {
  Y_missing <- list(matrix(0, nrow = 2, ncol = 2))
  expect_error(
    as_dynamic_network(Y_missing),
    "row names and column names"
  )

  Y_mismatch <- list(matrix(0, nrow = 2, ncol = 2))
  dimnames(Y_mismatch[[1]]) <- list(c("A", "B"), c("A", "C"))
  expect_error(
    as_dynamic_network(Y_mismatch),
    "match exactly"
  )
})

test_that("self ties are handled explicitly when disabled", {
  Y <- list(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))
  dimnames(Y[[1]]) <- list(c("A", "B"), c("A", "B"))

  dn <- as_dynamic_network(Y, self_ties = FALSE)

  diag_rows <- dn$dyads$sender == dn$dyads$receiver
  expect_true(any(diag_rows))
  expect_true(all(is.na(dn$dyads$value[diag_rows])))
  expect_false(any(dn$dyads$observed[diag_rows]))
})

test_that("print.dynamic_network runs without error", {
  Y <- list(matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE))
  dimnames(Y[[1]]) <- list(c("A", "B"), c("A", "B"))

  dn <- as_dynamic_network(Y)
  expect_silent(capture.output(print(dn)))
  expect_silent(capture.output(summary(dn)))
})
