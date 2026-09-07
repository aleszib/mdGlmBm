# Deterministic, optimizer-independent fixtures for integrated validation.

.with_test_seed <- function(seed, code) {
  had_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_state) old_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_state) assign(".Random.seed", old_state, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

.adjusted_rand_index <- function(x, y) {
  if (length(x) != length(y)) stop("Partitions must have equal length.", call. = FALSE)
  n <- length(x)
  if (n < 2L) return(1)
  tab <- table(x, y)
  choose_two <- function(z) sum(z * (z - 1) / 2)
  total <- choose_two(n)
  a <- choose_two(rowSums(tab))
  b <- choose_two(colSums(tab))
  c <- choose_two(tab)
  expected <- a * b / total
  maximum <- (a + b) / 2
  if (maximum == expected) return(if (c == maximum) 1 else 0)
  (c - expected) / (maximum - expected)
}

.make_truth_membership <- function(network, truth_by_time) {
  data.frame(
    unit_id = network$actor_time$unit_id,
    membership = as.vector(t(truth_by_time)),
    stringsAsFactors = FALSE
  )
}

simulate_binomial_validation_fixture <- function(seed = 5L) {
  .with_test_seed(seed, {
    n <- 16L
    n_time <- 4L
    actors <- sprintf("a%02d", seq_len(n))
    truth <- matrix(rep(rep(1:2, each = n / 2L), n_time), nrow = n_time, byrow = TRUE)
    Y <- lapply(seq_len(n_time), function(i) {
      z <- truth[i, ]
      probability <- ifelse(outer(z, z, "=="), 0.85, 0.08)
      diag(probability) <- 0
      M <- matrix(rbinom(n * n, size = 1L, prob = probability), nrow = n)
      diag(M) <- 0
      dimnames(M) <- list(actors, actors)
      M
    })
    names(Y) <- paste0("t", seq_len(n_time))
    network <- as_dynamic_network(Y, directed = TRUE, self_ties = FALSE)
    list(
      network = network,
      truth = truth,
      truth_membership = .make_truth_membership(network, truth),
      within_probability = 0.85,
      between_probability = 0.08,
      k = 2L
    )
  })
}

simulate_entry_exit_fixture <- function(seed = 900L) {
  .with_test_seed(seed, {
    actors_by_time <- list(
      c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"),
      c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "I"),
      c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "I", "J")
    )
    group <- c(setNames(rep(1L, 4L), actors_by_time[[1L]][1:4]),
               setNames(rep(2L, 4L), actors_by_time[[1L]][5:8]), I = 2L, J = 1L)
    Y <- lapply(actors_by_time, function(actors) {
      z <- unname(group[actors])
      probability <- ifelse(outer(z, z, "=="), 0.82, 0.10)
      diag(probability) <- 0
      M <- matrix(rbinom(length(actors)^2, size = 1L, prob = probability),
                  nrow = length(actors))
      diag(M) <- 0
      dimnames(M) <- list(actors, actors)
      M
    })
    names(Y) <- paste0("t", seq_along(Y))
    network <- as_dynamic_network(Y, directed = TRUE, self_ties = FALSE)
    truth <- data.frame(
      unit_id = network$actor_time$unit_id,
      membership = unname(group[network$actor_time$actor_id]),
      stringsAsFactors = FALSE
    )
    list(network = network, truth_membership = truth, k = 2L)
  })
}

simulate_ppml_validation_fixture <- function(seed = 88L) {
  .with_test_seed(seed, {
    n <- 8L
    actors <- paste0("p", seq_len(n))
    truth <- rep(1:2, each = n / 2L)
    means <- ifelse(outer(truth, truth, "=="), 3, 0.4)
    diag(means) <- 0
    Y <- lapply(seq_len(2L), function(i) {
      M <- matrix(rpois(n * n, means), nrow = n)
      diag(M) <- 0
      dimnames(M) <- list(actors, actors)
      M
    })
    names(Y) <- c("t1", "t2")
    network <- as_dynamic_network(Y, directed = TRUE, self_ties = FALSE)
    list(
      network = network,
      truth_membership = .make_truth_membership(network, matrix(truth, nrow = 2L, byrow = TRUE)),
      k = 2L,
      within_mean = 3,
      between_mean = 0.4
    )
  })
}
