#' Estimate Markov transition probabilities from a dynamic network
#'
#' Estimate cluster-to-cluster transition counts and smoothed transition
#' probabilities from the identity lineage edges in a `dynamic_network` object.
#'
#' @param x A `dynamic_network` object from [as_dynamic_network()].
#' @param membership Current actor-time memberships aligned to `x$actor_time`.
#'   Supported inputs follow the same conventions as
#'   [fit_time_glm_blockmodels()].
#' @param k Optional number of clusters. If omitted, the number is inferred
#'   from the supplied memberships.
#' @param smoothing Additive smoothing constant. The default `0.5` avoids zero
#'   transition probabilities when a row has sparse or empty counts.
#'
#' @return A list of class `markov_transitions` with elements:
#' \describe{
#'   \item{counts}{Transition count matrix.}
#'   \item{probabilities}{Row-normalized transition probability matrix.}
#'   \item{penalties}{Transition penalties on the deviance scale.}
#'   \item{labels}{Cluster labels used for the matrix dimensions.}
#'   \item{smoothing}{Additive smoothing constant.}
#'   \item{n_lineage_edges_used}{Number of lineage edges included.}
#'   \item{membership}{Normalized actor-time membership table.}
#'   \item{criterion_note}{Short note describing the estimate.}
#' }
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1,
#'                 1, 0), nrow = 2, byrow = TRUE),
#'   t2 = matrix(c(0, 1,
#'                 1, 0), nrow = 2, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dimnames(Y$t2) <- list(c("A", "B"), c("A", "B"))
#' dn <- as_dynamic_network(Y)
#' estimate_markov_transitions(dn, membership = c(1L, 2L, 2L, 1L))
#'
#' @export
estimate_markov_transitions <- function(x, membership, k = NULL, smoothing = 0.5) {
  if (!inherits(x, "dynamic_network")) {
    stop("`x` must be a `dynamic_network` object.", call. = FALSE)
  }
  if (!is.numeric(smoothing) || length(smoothing) != 1L || is.na(smoothing) || smoothing < 0) {
    stop("`smoothing` must be a single non-negative numeric value.", call. = FALSE)
  }

  membership_table <- .time_glm_normalize_membership(
    membership = membership,
    actor_time = x$actor_time,
    time_labels = x$times
  )
  labels <- .transition_cluster_labels(membership_table$membership, k = k)
  labels_chr <- as.character(labels)

  counts <- matrix(0, nrow = length(labels), ncol = length(labels),
                   dimnames = list(labels_chr, labels_chr))

  lineage <- x$lineage
  if (nrow(lineage) > 0L) {
    lineage <- lineage[lineage$relation == "identity", , drop = FALSE]
  }
  if (nrow(lineage) > 0L) {
    mem_lookup <- membership_table$membership
    names(mem_lookup) <- as.character(membership_table$unit_id)
    from_cluster <- unname(mem_lookup[as.character(lineage$from_unit)])
    to_cluster <- unname(mem_lookup[as.character(lineage$to_unit)])

    if (anyNA(from_cluster) || anyNA(to_cluster)) {
      stop(
        "Memberships must be supplied for all actor-time units referenced by the lineage table.",
        call. = FALSE
      )
    }

    valid <- lineage$weight > 0 & is.finite(lineage$weight)
    if (any(valid)) {
      from_factor <- factor(from_cluster[valid], levels = labels)
      to_factor <- factor(to_cluster[valid], levels = labels)
      transition_counts <- xtabs(
        lineage$weight[valid] ~ from_factor + to_factor,
        drop.unused.levels = FALSE
      )
      counts[rownames(transition_counts), colnames(transition_counts)] <- transition_counts
    }
  }

  probabilities <- counts + smoothing
  row_sums <- rowSums(probabilities)
  probabilities <- sweep(probabilities, 1, row_sums, "/")
  penalties <- .transition_penalty(probabilities)

  out <- list(
    counts = counts,
    probabilities = probabilities,
    penalties = penalties,
    labels = labels,
    smoothing = smoothing,
    n_lineage_edges_used = nrow(lineage),
    membership = membership_table,
    criterion_note = "Transition probabilities were estimated from identity lineage edges with additive smoothing.",
    call = match.call()
  )
  class(out) <- "markov_transitions"
  out
}

#' @keywords internal
.transition_cluster_labels <- function(membership, k = NULL) {
  membership <- as.integer(membership)
  if (anyNA(membership)) {
    stop("Membership labels must be finite integers.", call. = FALSE)
  }
  if (any(membership < 1L)) {
    stop("Membership labels must be one-based positive integers.", call. = FALSE)
  }

  inferred <- max(membership)
  if (is.null(k)) {
    k <- inferred
  } else {
    if (length(k) != 1L || is.na(k) || k < 1L) {
      stop("`k` must be a single positive integer.", call. = FALSE)
    }
    k <- as.integer(k)
    if (k < inferred) {
      stop("`k` cannot be smaller than the largest membership label.", call. = FALSE)
    }
  }

  seq_len(k)
}

#' @keywords internal
.transition_penalty <- function(probability) {
  probability <- as.numeric(probability)
  out <- rep(Inf, length(probability))
  valid <- is.finite(probability) & probability > 0
  out[valid] <- -2 * log(probability[valid])
  out
}

#' Estimate a membership prior for dynamic scoring
#'
#' Estimate a prior distribution over cluster memberships for use in deviance-
#' scale local scoring.
#'
#' @param membership Membership labels or a membership table. Character and
#'   factor inputs are treated as factor-like labels and converted to one-based
#'   integers.
#' @param prior Prior mode: `"empirical"`, `"uniform"`, or `"none"`.
#' @param k Optional number of clusters. If omitted, the value is inferred from
#'   the supplied memberships.
#' @param smoothing Additive smoothing used for the empirical prior.
#' @param labels Optional cluster labels to carry into the result.
#'
#' @return A list of class `membership_prior` with elements:
#' \describe{
#'   \item{prior}{Prior mode.}
#'   \item{labels}{Cluster labels.}
#'   \item{counts}{Cluster counts or `NA` when the prior is disabled.}
#'   \item{probabilities}{Prior probabilities.}
#'   \item{penalties}{Prior penalties on the deviance scale.}
#'   \item{smoothing}{Additive smoothing constant.}
#'   \item{n}{Number of memberships used.}
#'   \item{criterion_note}{Short note describing the prior.}
#' }
#'
#' @examples
#' estimate_membership_prior(c(1L, 1L, 2L, 3L), prior = "empirical")
#'
#' @export
estimate_membership_prior <- function(membership,
                                      prior = c("empirical", "uniform", "none"),
                                      k = NULL,
                                      smoothing = 0.5,
                                      labels = NULL) {
  prior <- match.arg(prior)
  membership <- .membership_to_integer_vector(membership)

  if (!is.numeric(smoothing) || length(smoothing) != 1L || is.na(smoothing) || smoothing < 0) {
    stop("`smoothing` must be a single non-negative numeric value.", call. = FALSE)
  }

  if (is.null(labels)) {
    labels <- .transition_cluster_labels(membership, k = k)
  } else {
    labels <- as.integer(labels)
    if (length(labels) == 0L) {
      stop("`labels` must contain at least one cluster label.", call. = FALSE)
    }
    if (anyNA(labels) || any(labels < 1L)) {
      stop("`labels` must be one-based positive integers.", call. = FALSE)
    }
    if (is.null(k)) {
      k <- length(labels)
    } else if (length(k) != 1L || is.na(k) || as.integer(k) != length(labels)) {
      stop("`k` must match the length of `labels` when both are supplied.", call. = FALSE)
    }
    labels <- as.integer(labels)
  }

  k <- length(labels)
  labels_chr <- as.character(labels)
  counts <- rep(NA_integer_, k)
  names(counts) <- labels_chr

  if (prior == "none") {
    probabilities <- rep(1 / k, k)
    names(probabilities) <- labels_chr
    penalties <- rep(0, k)
    names(penalties) <- labels_chr
    criterion_note <- "Membership prior disabled; prior penalties are zero."
  } else {
    counts <- tabulate(factor(membership, levels = labels), nbins = k)
    names(counts) <- labels_chr
    if (prior == "uniform") {
      probabilities <- rep(1 / k, k)
    } else {
      probabilities <- (counts + smoothing) / (sum(counts) + k * smoothing)
    }
    names(probabilities) <- labels_chr
    penalties <- .transition_penalty(probabilities)
    names(penalties) <- labels_chr
    criterion_note <- if (prior == "uniform") {
      "Uniform membership prior."
    } else {
      "Empirical membership prior with additive smoothing."
    }
  }

  out <- list(
    prior = prior,
    labels = labels,
    counts = counts,
    probabilities = probabilities,
    penalties = penalties,
    smoothing = smoothing,
    n = length(membership),
    criterion_note = criterion_note,
    call = match.call()
  )
  class(out) <- "membership_prior"
  out
}

#' @keywords internal
.membership_to_integer_vector <- function(membership) {
  if (is.data.frame(membership)) {
    col_membership <- intersect(c("membership", "cluster", "clu", "label"), names(membership))
    col_membership <- if (length(col_membership)) col_membership[[1]] else NULL
    if (is.null(col_membership)) {
      stop(
        "Membership data frames must contain a membership column named one of: membership, cluster, clu, label.",
        call. = FALSE
      )
    }
    membership <- membership[[col_membership]]
  } else if (is.list(membership) && !is.null(membership$membership)) {
    membership <- membership$membership
  }

  if (is.factor(membership)) {
    membership <- as.integer(membership)
  } else if (is.character(membership)) {
    membership <- as.integer(factor(membership))
  } else {
    membership <- suppressWarnings(as.numeric(membership))
  }

  if (anyNA(membership)) {
    stop("Membership labels must be finite, positive integers or factor-like values.", call. = FALSE)
  }
  membership <- as.integer(round(membership))
  if (any(membership < 1L)) {
    stop("Membership labels must be one-based positive integers.", call. = FALSE)
  }
  membership
}

#' @keywords internal
.markov_transition_matrix <- function(x) {
  if (inherits(x, "markov_transitions")) {
    return(x$probabilities)
  }
  if (!is.null(x$probabilities)) {
    return(x$probabilities)
  }
  if (!is.null(x$transition_probabilities)) {
    return(x$transition_probabilities)
  }
  stop("`transition` must be a `markov_transitions` object or contain transition probabilities.", call. = FALSE)
}

#' @keywords internal
.markov_transition_labels <- function(x) {
  if (inherits(x, "markov_transitions") && !is.null(x$labels)) {
    return(as.integer(x$labels))
  }
  if (!is.null(x$labels)) {
    return(as.integer(x$labels))
  }
  probs <- .markov_transition_matrix(x)
  if (is.null(dim(probs))) {
    stop("Transition probabilities must be a matrix-like object.", call. = FALSE)
  }
  seq_len(nrow(probs))
}

#' @keywords internal
.membership_prior_probs <- function(x) {
  if (inherits(x, "membership_prior") && !is.null(x$probabilities)) {
    return(x$probabilities)
  }
  if (!is.null(x$probabilities)) {
    return(x$probabilities)
  }
  stop("`prior` must be a `membership_prior` object or contain prior probabilities.", call. = FALSE)
}

#' @keywords internal
.membership_prior_labels <- function(x) {
  if (inherits(x, "membership_prior") && !is.null(x$labels)) {
    return(as.integer(x$labels))
  }
  if (!is.null(x$labels)) {
    return(as.integer(x$labels))
  }
  probs <- .membership_prior_probs(x)
  seq_len(length(probs))
}

#' @keywords internal
.resolve_time_glm_fit <- function(fit, time_index, time_label) {
  if (inherits(fit, "time_glm_blockmodels")) {
    if (time_index < 1L || time_index > length(fit$fits)) {
      stop("`time_index` is out of range for the supplied fit object.", call. = FALSE)
    }
    selected <- fit$fits[[time_index]]
    if (!identical(as.character(selected$time), as.character(time_label))) {
      stop("The supplied fit object does not match the target time point.", call. = FALSE)
    }
    return(selected)
  }

  if (inherits(fit, "time_glm_blockmodel")) {
    if (!identical(as.character(fit$time), as.character(time_label)) && !identical(fit$time_index, time_index)) {
      stop("The supplied fit object does not match the target time point.", call. = FALSE)
    }
    return(fit)
  }

  stop("`fit` must be a `time_glm_blockmodel` or `time_glm_blockmodels` object.", call. = FALSE)
}

#' @keywords internal
.time_glm_candidate_deviance <- function(fit, candidate_data) {
  if (nrow(candidate_data) == 0L) {
    return(0)
  }

  family_obj <- stats::family(fit$fit)
  has_block <- !is.null(fit$fit$xlevels) && !is.null(fit$fit$xlevels$block)
  block_levels <- if (has_block) fit$fit$xlevels$block else character(0)
  known <- rep(TRUE, nrow(candidate_data))
  if (has_block) {
    known <- candidate_data$block %in% block_levels
  }

  dev <- numeric(nrow(candidate_data))
  if (any(known)) {
    known_data <- candidate_data[known, , drop = FALSE]
    if (has_block) {
      known_data$block <- factor(known_data$block, levels = block_levels)
    }
    mu_known <- stats::predict(fit$fit, newdata = known_data, type = "response")
    dev[known] <- family_obj$dev.resids(known_data$value, mu_known, wt = known_data$weight)
  }

  if (any(!known)) {
    fallback_fit <- if (has_block) {
      fallback_formula <- stats::update(stats::formula(fit$fit), . ~ . - block)
      stats::glm(fallback_formula, data = fit$data, family = stats::family(fit$fit))
    } else {
      fit$fit
    }
    unknown_data <- candidate_data[!known, , drop = FALSE]
    mu_unknown <- stats::predict(fallback_fit, newdata = unknown_data, type = "response")
    dev[!known] <- stats::family(fallback_fit)$dev.resids(unknown_data$value, mu_unknown, wt = unknown_data$weight)
  }

  sum(dev)
}

#' Score candidate memberships for one actor-time unit
#'
#' Compute deviance-scale local scores for assigning a single actor-time unit
#' to candidate clusters. The score combines the time-specific GLM deviance
#' contribution with Markov transition penalties and an optional membership
#' prior penalty.
#'
#' @param x A `dynamic_network` object from [as_dynamic_network()].
#' @param fit A fitted time-GLM object from [fit_time_glm_blockmodel()] or
#'   [fit_time_glm_blockmodels()].
#' @param membership Current actor-time memberships. When omitted, the
#'   memberships stored in `fit` are used if they align with `x$actor_time`.
#' @param unit_id Optional actor-time `unit_id` identifying the target row.
#' @param row_index Optional one-based row index into `x$actor_time`.
#' @param candidate_clusters Candidate cluster labels. Defaults to all labels
#'   implied by the transition or prior object, or by the current memberships.
#' @param transition Optional pre-estimated Markov transition object. If not
#'   supplied, `estimate_markov_transitions()` is called.
#' @param prior Optional membership prior specification or
#'   `membership_prior` object.
#' @param prior_smoothing Additive smoothing used when estimating an empirical
#'   membership prior.
#' @param smoothing Additive smoothing used when estimating transitions if no
#'   transition object is supplied.
#' @param return_components Logical; kept for future compatibility. The
#'   reference implementation always returns the component columns.
#'
#' @return A data frame with one row per candidate cluster and at least:
#' `unit_id`, `actor_id`, `time`, `time_index`, `candidate_cluster`,
#' `glm_deviance`, `previous_transition_penalty`, `next_transition_penalty`,
#' `prior_penalty`, and `total_score`.
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1,
#'                 1, 0), nrow = 2, byrow = TRUE),
#'   t2 = matrix(c(0, 1,
#'                 1, 0), nrow = 2, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dimnames(Y$t2) <- list(c("A", "B"), c("A", "B"))
#' dn <- as_dynamic_network(Y)
#' mem <- c(1L, 2L, 2L, 1L)
#' fit <- fit_time_glm_blockmodels(dn, membership = mem, family = "binomial")
#' score_actor_time_candidates(dn, fit, membership = mem, row_index = 1L)
#'
#' @export
score_actor_time_candidates <- function(x, fit, membership = NULL,
                                        unit_id = NULL, row_index = NULL,
                                        candidate_clusters = NULL,
                                        transition = NULL,
                                        prior = c("none", "uniform", "empirical"),
                                        prior_smoothing = 0.5,
                                        smoothing = 0.5,
                                        return_components = TRUE) {
  if (!inherits(x, "dynamic_network")) {
    stop("`x` must be a `dynamic_network` object.", call. = FALSE)
  }

  if (!is.null(unit_id) && !is.null(row_index)) {
    stop("Specify only one of `unit_id` or `row_index`.", call. = FALSE)
  }

  if (!is.null(row_index)) {
    if (length(row_index) != 1L || is.na(row_index) || row_index < 1L || row_index > nrow(x$actor_time)) {
      stop("`row_index` is out of range.", call. = FALSE)
    }
    target <- x$actor_time[as.integer(row_index), , drop = FALSE]
  } else if (!is.null(unit_id)) {
    idx <- match(unit_id, x$actor_time$unit_id)
    if (is.na(idx)) {
      stop("`unit_id` was not found in the actor-time table.", call. = FALSE)
    }
    target <- x$actor_time[idx, , drop = FALSE]
  } else {
    stop("Either `unit_id` or `row_index` must be supplied.", call. = FALSE)
  }

  if (is.null(membership)) {
    if (!is.null(fit$membership) && is.data.frame(fit$membership) && nrow(fit$membership) == nrow(x$actor_time)) {
      membership <- fit$membership
    } else if (!is.null(fit$membership) && is.vector(fit$membership) && length(fit$membership) == nrow(x$actor_time)) {
      membership <- fit$membership
    } else {
      stop("`membership` must be supplied when the fit object does not store full actor-time memberships.", call. = FALSE)
    }
  }

  membership_table <- .time_glm_normalize_membership(
    membership = membership,
    actor_time = x$actor_time,
    time_labels = x$times
  )
  mem_lookup <- membership_table$membership
  names(mem_lookup) <- as.character(membership_table$unit_id)

  time_fit <- .resolve_time_glm_fit(fit, target$time_index, target$time)

  if (is.null(transition)) {
    transition <- estimate_markov_transitions(
      x = x,
      membership = membership_table,
      k = NULL,
      smoothing = smoothing
    )
  }
  transition_probs <- .markov_transition_matrix(transition)
  transition_labels <- .markov_transition_labels(transition)

  prior_mode <- prior
  if (inherits(prior, "membership_prior")) {
    prior_obj <- prior
  } else {
    prior_mode <- match.arg(prior)
    prior_obj <- estimate_membership_prior(
      membership = membership_table$membership,
      prior = prior_mode,
      k = length(transition_labels),
      smoothing = prior_smoothing,
      labels = transition_labels
    )
  }
  prior_probs <- .membership_prior_probs(prior_obj)
  prior_penalties <- .membership_prior_penalty(prior_obj)

  if (is.null(candidate_clusters)) {
    candidate_clusters <- transition_labels
  } else {
    candidate_clusters <- as.integer(candidate_clusters)
    if (anyNA(candidate_clusters) || any(candidate_clusters < 1L)) {
      stop("`candidate_clusters` must contain one-based positive integers.", call. = FALSE)
    }
  }

  if (any(!candidate_clusters %in% transition_labels)) {
    stop("`candidate_clusters` must be among the labels used by the transition object.", call. = FALSE)
  }

  candidate_clusters <- sort(unique(candidate_clusters))
  candidate_clusters_chr <- as.character(candidate_clusters)

  lineages_prev <- x$lineage[x$lineage$to_unit == target$unit_id & x$lineage$relation == "identity", , drop = FALSE]
  lineages_next <- x$lineage[x$lineage$from_unit == target$unit_id & x$lineage$relation == "identity", , drop = FALSE]

  candidate_blocks <- character(0)
  if (nrow(time_fit$data) > 0L) {
    target_rows <- time_fit$data$sender_unit == target$unit_id | time_fit$data$receiver_unit == target$unit_id
    candidate_data_base <- time_fit$data[target_rows, , drop = FALSE]
    if (nrow(candidate_data_base) > 0L) {
      for (candidate in candidate_clusters) {
        candidate_lookup <- mem_lookup
        candidate_lookup[as.character(target$unit_id)] <- candidate
        candidate_data <- candidate_data_base
        candidate_data$sender_membership <- unname(candidate_lookup[as.character(candidate_data$sender_unit)])
        candidate_data$receiver_membership <- unname(candidate_lookup[as.character(candidate_data$receiver_unit)])
        candidate_data$block <- paste(candidate_data$sender_membership, candidate_data$receiver_membership, sep = "#")
        candidate_blocks <- c(candidate_blocks, as.character(unique(candidate_data$block)))
      }
    }
  }
  candidate_blocks <- unique(candidate_blocks)

  score_rows <- vector("list", length(candidate_clusters))
  for (i in seq_along(candidate_clusters)) {
    candidate <- candidate_clusters[i]
    candidate_lookup <- mem_lookup
    candidate_lookup[as.character(target$unit_id)] <- candidate

    candidate_data <- time_fit$data[
      time_fit$data$sender_unit == target$unit_id | time_fit$data$receiver_unit == target$unit_id,
      ,
      drop = FALSE
    ]
    if (nrow(candidate_data) > 0L) {
      candidate_data$sender_membership <- unname(candidate_lookup[as.character(candidate_data$sender_unit)])
      candidate_data$receiver_membership <- unname(candidate_lookup[as.character(candidate_data$receiver_unit)])
      candidate_data$block <- paste(candidate_data$sender_membership, candidate_data$receiver_membership, sep = "#")
    }

    glm_deviance <- .time_glm_candidate_deviance(time_fit, candidate_data)

    prev_penalty <- 0
    if (nrow(lineages_prev) > 0L) {
      prev_clusters <- unname(mem_lookup[as.character(lineages_prev$from_unit)])
      prev_probs <- transition_probs[cbind(as.character(prev_clusters), as.character(candidate))]
      prev_penalty <- sum(.transition_penalty(prev_probs))
    }

    next_penalty <- 0
    if (nrow(lineages_next) > 0L) {
      next_clusters <- unname(mem_lookup[as.character(lineages_next$to_unit)])
      next_probs <- transition_probs[cbind(as.character(candidate), as.character(next_clusters))]
      next_penalty <- sum(.transition_penalty(next_probs))
    }

    prior_penalty <- 0
    if (identical(prior_obj$prior, "none")) {
      prior_penalty <- 0
    } else {
      prior_penalty <- .membership_prior_penalty(prior_obj)[as.character(candidate)]
      prior_penalty <- as.numeric(prior_penalty)
    }

    total_score <- glm_deviance + prev_penalty + next_penalty + prior_penalty
    score_rows[[i]] <- data.frame(
      unit_id = target$unit_id,
      actor_id = target$actor_id,
      time = target$time,
      time_index = target$time_index,
      candidate_cluster = candidate,
      glm_deviance = glm_deviance,
      previous_transition_penalty = prev_penalty,
      next_transition_penalty = next_penalty,
      prior_penalty = prior_penalty,
      total_score = total_score,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, score_rows)
  rownames(out) <- NULL
  attr(out, "transition") <- transition
  attr(out, "prior") <- prior_obj
  attr(out, "fit_time_index") <- time_fit$time_index
  attr(out, "fit_time_label") <- time_fit$time
  out
}

#' @keywords internal
.membership_prior_penalty <- function(x) {
  probs <- .membership_prior_probs(x)
  penalties <- .transition_penalty(probs)
  names(penalties) <- names(probs)
  penalties
}

#' Print a Markov transition estimate
#'
#' @param x A `markov_transitions` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.markov_transitions <- function(x, ...) {
  cat("markov_transitions object\n")
  cat(sprintf("  labels: %s\n", paste(x$labels, collapse = ", ")))
  cat(sprintf("  smoothing: %s\n", format(x$smoothing)))
  cat(sprintf("  lineage edges used: %s\n", format(x$n_lineage_edges_used)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a Markov transition estimate
#'
#' @param object A `markov_transitions` object.
#' @param ... Ignored.
#'
#' @return A compact summary list.
#' @export
summary.markov_transitions <- function(object, ...) {
  summary <- list(
    labels = object$labels,
    smoothing = object$smoothing,
    n_lineage_edges_used = object$n_lineage_edges_used,
    criterion_note = object$criterion_note
  )
  class(summary) <- "summary.markov_transitions"
  summary
}

#' @export
print.summary.markov_transitions <- function(x, ...) {
  cat("Summary of markov_transitions\n")
  cat(sprintf("  labels: %s\n", paste(x$labels, collapse = ", ")))
  cat(sprintf("  smoothing: %s\n", format(x$smoothing)))
  cat(sprintf("  lineage edges used: %s\n", format(x$n_lineage_edges_used)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Print a membership prior estimate
#'
#' @param x A `membership_prior` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.membership_prior <- function(x, ...) {
  cat("membership_prior object\n")
  cat(sprintf("  prior: %s\n", x$prior))
  cat(sprintf("  labels: %s\n", paste(x$labels, collapse = ", ")))
  cat(sprintf("  smoothing: %s\n", format(x$smoothing)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a membership prior estimate
#'
#' @param object A `membership_prior` object.
#' @param ... Ignored.
#'
#' @return A compact summary list.
#' @export
summary.membership_prior <- function(object, ...) {
  summary <- list(
    prior = object$prior,
    labels = object$labels,
    smoothing = object$smoothing,
    criterion_note = object$criterion_note
  )
  class(summary) <- "summary.membership_prior"
  summary
}

#' @export
print.summary.membership_prior <- function(x, ...) {
  cat("Summary of membership_prior\n")
  cat(sprintf("  prior: %s\n", x$prior))
  cat(sprintf("  labels: %s\n", paste(x$labels, collapse = ", ")))
  cat(sprintf("  smoothing: %s\n", format(x$smoothing)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}
