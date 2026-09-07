make_transition_fixture <- function() {
  Y <- list(
    t1 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE),
    t2 = matrix(c(0, 1,
                  1, 0), nrow = 2, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
  dimnames(Y$t2) <- list(c("A", "B"), c("A", "B"))
  as_dynamic_network(Y)
}

make_scoring_fixture <- function() {
  Y <- list(
    t1 = matrix(c(0, 1, 0,
                  0, 0, 1,
                  1, 0, 0), nrow = 3, byrow = TRUE),
    t2 = matrix(c(0, 1, 1,
                  1, 0, 0,
                  0, 1, 0), nrow = 3, byrow = TRUE),
    t3 = matrix(c(0, 0, 1,
                  1, 0, 1,
                  0, 1, 0), nrow = 3, byrow = TRUE)
  )
  dimnames(Y$t1) <- list(c("A", "B", "C"), c("A", "B", "C"))
  dimnames(Y$t2) <- list(c("A", "C", "D"), c("A", "C", "D"))
  dimnames(Y$t3) <- list(c("A", "C", "D"), c("A", "C", "D"))
  as_dynamic_network(Y)
}

make_scoring_membership <- function(dn) {
  data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L),
    stringsAsFactors = FALSE
  )
}

test_that("estimate_markov_transitions returns expected counts on a tiny fixed-actor network", {
  dn <- make_transition_fixture()
  membership <- c(1L, 2L, 2L, 1L)

  transitions <- estimate_markov_transitions(dn, membership = membership, smoothing = 0)

  expect_s3_class(transitions, "markov_transitions")
  expect_equal(transitions$labels, 1:2)
  expect_equal(transitions$n_lineage_edges_used, 2L)
  expect_equal(
    transitions$counts,
    matrix(c(0, 1, 0, 0,
             1, 0, 0, 0,
             0, 0, 0, 0,
             0, 0, 0, 0), nrow = 4, byrow = TRUE,
           dimnames = list(c("1", "2", "E", "V"), c("1", "2", "E", "V")))
  )
  expect_equal(transitions$counts_by_boundary[[1]]["1", "2"], 1)
  expect_equal(transitions$counts_by_boundary[[1]]["2", "1"], 1)
})

test_that("additive smoothing keeps transition probabilities positive", {
  dn <- make_transition_fixture()
  membership <- c(1L, 2L, 2L, 1L)

  transitions <- estimate_markov_transitions(dn, membership = membership, k = 3, smoothing = 0.5)

  expect_true(all(transitions$probabilities[transitions$probabilities > 0] > 0))
  expect_equal(transitions$probabilities["1", "E"], 0)
  expect_equal(transitions$probabilities["E", "V"], 0)
  expect_equal(nrow(transitions$probabilities), 5L)
  expect_equal(ncol(transitions$probabilities), 5L)
})

test_that("transition penalties equal -2 log probability", {
  dn <- make_scoring_fixture()
  membership <- make_scoring_membership(dn)
  transitions <- estimate_markov_transitions(dn, membership = membership, smoothing = 0.5)
  prior <- estimate_membership_prior(membership$membership, prior = "none")

  target_row <- which(dn$actor_time$time == "t2" & dn$actor_time$actor_id == "A")
  scores <- score_actor_time_candidates(
    dn,
    fit_time_glm_blockmodels(dn, membership = membership, family = "binomial"),
    membership = membership,
    row_index = target_row,
    candidate_clusters = c(1L, 2L),
    transition = transitions,
    prior = prior
  )

  row_two <- scores[scores$candidate_cluster == 2L, , drop = FALSE]
  prev_prob <- transitions$probabilities_by_boundary[[1]]["1", "2"]
  next_prob <- transitions$probabilities_by_boundary[[2]]["2", "1"]
  expect_equal(row_two$previous_transition_penalty, -2 * log(prev_prob), tolerance = 1e-8)
  expect_equal(row_two$next_transition_penalty, -2 * log(next_prob), tolerance = 1e-8)
})

test_that("entry and exit actor-time units use augmented transition penalties", {
  dn <- make_scoring_fixture()
  membership <- make_scoring_membership(dn)
  transitions <- estimate_markov_transitions(dn, membership = membership, smoothing = 0.5)
  prior <- estimate_membership_prior(membership$membership, prior = "uniform")
  fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")

  entry_row <- which(dn$actor_time$time == "t2" & dn$actor_time$actor_id == "D")
  entry_scores <- score_actor_time_candidates(
    dn, fit, membership = membership, row_index = entry_row,
    candidate_clusters = c(1L, 2L), transition = transitions, prior = prior
  )
  expect_equal(
    entry_scores$previous_transition_penalty,
    vapply(entry_scores$candidate_cluster, function(candidate) {
      -2 * log(transitions$probabilities_by_boundary[[1]]["E", as.character(candidate)])
    }, numeric(1)), tolerance = 1e-8
  )
  expect_true(all(is.finite(entry_scores$total_score)))

  exit_row <- which(dn$actor_time$time == "t1" & dn$actor_time$actor_id == "B")
  exit_scores <- score_actor_time_candidates(
    dn, fit, membership = membership, row_index = exit_row,
    candidate_clusters = c(1L, 2L), transition = transitions, prior = prior
  )
  expect_equal(
    exit_scores$next_transition_penalty,
    vapply(exit_scores$candidate_cluster, function(candidate) {
      -2 * log(transitions$probabilities_by_boundary[[1]][as.character(candidate), "V"])
    }, numeric(1)), tolerance = 1e-8
  )
  expect_true(all(is.finite(exit_scores$total_score)))
})

test_that("membership prior helper supports uniform and empirical modes", {
  uniform <- estimate_membership_prior(c(1L, 1L, 2L, 3L), prior = "uniform")
  empirical <- estimate_membership_prior(c(1L, 1L, 2L, 3L), prior = "empirical")

  expect_identical(uniform$prior, "uniform")
  expect_equal(unname(uniform$penalties), rep(-2 * log(1 / 3), 3), tolerance = 1e-8)
  expect_identical(empirical$prior, "empirical")
  expect_equal(sum(empirical$probabilities), 1)
  expect_true(all(empirical$penalties >= 0))
})

test_that("score_actor_time_candidates returns one row per candidate cluster and the required columns", {
  dn <- make_scoring_fixture()
  membership <- make_scoring_membership(dn)
  fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")
  transitions <- estimate_markov_transitions(dn, membership = membership, smoothing = 0.5)
  prior <- estimate_membership_prior(membership$membership, prior = "empirical")

  scores <- score_actor_time_candidates(
    dn,
    fit,
    membership = membership,
    row_index = which(dn$actor_time$time == "t2" & dn$actor_time$actor_id == "A"),
    candidate_clusters = 1:2,
    transition = transitions,
    prior = prior
  )

  expect_equal(nrow(scores), 2L)
  expect_true(all(c(
    "unit_id",
    "actor_id",
    "time",
    "time_index",
    "candidate_cluster",
    "glm_deviance",
    "previous_transition_penalty",
    "next_transition_penalty",
    "prior_penalty",
    "total_score"
  ) %in% names(scores)))
})

test_that("score_actor_time_candidates accepts both binomial and PPML time-GLM fits", {
  dn <- make_scoring_fixture()
  membership <- make_scoring_membership(dn)
  row_index <- which(dn$actor_time$time == "t2" & dn$actor_time$actor_id == "A")

  bin_fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")
  ppml_fit <- fit_time_glm_blockmodels(dn, membership = membership, family = "ppml")
  transitions <- estimate_markov_transitions(dn, membership = membership, smoothing = 0.5)

  bin_scores <- score_actor_time_candidates(
    dn, bin_fit, membership = membership, row_index = row_index,
    candidate_clusters = 1:2, transition = transitions, prior = "none"
  )
  ppml_scores <- score_actor_time_candidates(
    dn, ppml_fit, membership = membership, row_index = row_index,
    candidate_clusters = 1:2, transition = transitions, prior = "none"
  )

  expect_true(is.data.frame(bin_scores))
  expect_true(is.data.frame(ppml_scores))
  expect_true(all(is.finite(bin_scores$total_score)))
  expect_true(all(is.finite(ppml_scores$total_score)))
})

test_that("augmented transition events have hand-checkable entry and exit counts", {
  Y <- lapply(list(c("A", "B", "C"), c("A", "B", "C", "D"),
                   c("A", "C", "D")), function(actors) {
    M <- matrix(0, length(actors), length(actors), dimnames = list(actors, actors))
    M
  })
  names(Y) <- paste0("t", seq_along(Y))
  dn <- as_dynamic_network(Y, self_ties = FALSE)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 2L),
    stringsAsFactors = FALSE
  )
  transitions <- estimate_markov_transitions(dn, membership, k = 2, smoothing = 0.5)

  expect_equal(transitions$counts_by_boundary[[1]]["1", "1"], 2)
  expect_equal(transitions$counts_by_boundary[[1]]["2", "2"], 1)
  expect_equal(transitions$counts_by_boundary[[1]]["E", "2"], 1)
  expect_equal(transitions$counts_by_boundary[[2]]["1", "1"], 2)
  expect_equal(transitions$counts_by_boundary[[2]]["2", "V"], 1)
  expect_equal(transitions$n_entries, 1L)
  expect_equal(transitions$n_exits, 1L)
  expect_equal(transitions$n_persistent, 6L)
  expect_false(any(transitions$events$actor_id == "D" & transitions$events$transition_type == "exit"))
  expect_false(isTRUE(all.equal(
    transitions$probabilities_by_boundary[[1]]["2", c("1", "2")],
    transitions$probabilities_by_boundary[[2]]["2", c("1", "2")]
  )))
  expect_equal(unname(rowSums(transitions$probabilities_by_boundary[[1]][c("1", "2"),
                                                                  c("1", "2")])), c(1, 1), tolerance = 1e-12)
  expect_true(is.finite(transitions$probabilities_by_boundary[[1]]["E", "2"]))
  expect_true(is.finite(transitions$probabilities_by_boundary[[2]]["2", "V"]))
})

test_that("initial probabilities apply only at time one and intermediate single-spell actors use E and V", {
  Y <- lapply(list(c("A", "B"), c("A", "B", "C"), c("A", "B")), function(actors) {
    matrix(0, length(actors), length(actors), dimnames = list(actors, actors))
  })
  names(Y) <- paste0("t", seq_along(Y))
  dn <- as_dynamic_network(Y, self_ties = FALSE)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 1L, 2L, 1L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  transitions <- estimate_markov_transitions(dn, membership, k = 2, smoothing = 0.5)
  expect_equal(transitions$counts_by_boundary[[1]]["E", "1"], 1)
  expect_equal(transitions$counts_by_boundary[[1]]["E", "2"], 0)
  expect_equal(transitions$counts_by_boundary[[2]]["1", "V"], 1)

  prior <- estimate_membership_prior(c(1L, 2L), prior = "empirical", k = 2, smoothing = 0.5)
  expect_equal(sum(prior$probabilities), 1)
  expect_true(is.finite(prior$penalties["1"]))
  prior_with_emergence <- estimate_membership_prior(
    c(1L, 1L, 2L), prior = "empirical", k = 2,
    smoothing = 0, emerging_count = 1
  )
  expect_equal(unname(prior_with_emergence$probabilities), c(2 / 4, 1 / 4))
  expect_lt(sum(prior_with_emergence$probabilities), 1)
})

test_that("an actor present only at an intermediate time has both boundary events", {
  Y <- lapply(list(c("A", "B"), c("A", "B", "C"), c("A", "B", "D")), function(actors) {
    n <- length(actors)
    M <- matrix(0, n, n, dimnames = list(actors, actors))
    if (n > 1L) M[row(M) != col(M)] <- 1
    M
  })
  dn <- as_dynamic_network(Y, self_ties = FALSE)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
    stringsAsFactors = FALSE
  )
  transitions <- estimate_markov_transitions(dn, membership, k = 2, smoothing = 0.5)
  events <- transitions$events[as.character(transitions$events$actor_id) == "C", , drop = FALSE]
  expect_equal(events$transition_type, c("entry", "exit"))
  expect_equal(events$from_state, c("E", "1"))
  expect_equal(events$to_state, c("1", "V"))
  expect_true(all(is.finite(c(
    transitions$probabilities_by_boundary[[1]]["E", "1"],
    transitions$probabilities_by_boundary[[2]]["1", "V"]
  ))))
})

test_that("re-entry is represented as separate exit and entry events", {
  actors_by_time <- list(c("A", "B"), c("B", "C"), c("A", "B", "C"))
  Y <- lapply(actors_by_time, function(actors) {
    matrix(0, length(actors), length(actors), dimnames = list(actors, actors))
  })
  dn <- as_dynamic_network(Y, self_ties = FALSE)
  membership <- data.frame(
    unit_id = dn$actor_time$unit_id,
    membership = c(1L, 2L, 2L, 1L, 1L, 2L, 1L),
    stringsAsFactors = FALSE
  )
  events <- estimate_markov_transitions(dn, membership, k = 2, smoothing = 0.5)$events
  expect_equal(events$transition_type[events$actor_id == "A"], c("exit", "entry"))
  expect_equal(events$from_state[events$actor_id == "A"], c("1", "E"))
  expect_equal(events$to_state[events$actor_id == "A"], c("V", "1"))
})
