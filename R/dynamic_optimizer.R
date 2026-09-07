#' Compute dynamic GLM-Markov objective components
#'
#' Internal helper used by the reference optimizer to keep the reporting scale
#' explicit and stable across iterations.
#'
#' @param network A `dynamic_network` object.
#' @param fits A `time_glm_blockmodels` object.
#' @param transition A `markov_transitions` object.
#' @param prior A `membership_prior` object.
#' @param membership Normalized actor-time membership table.
#'
#' @return A list with the fitted observation-model and dynamic objective
#'   components.
#' @keywords internal
.dynamic_glm_objective_components <- function(network, fits, transition, prior, membership) {
  if (!inherits(network, "dynamic_network")) {
    stop("`network` must be a `dynamic_network` object.", call. = FALSE)
  }
  if (!inherits(fits, "time_glm_blockmodels")) {
    stop("`fits` must be a `time_glm_blockmodels` object.", call. = FALSE)
  }

  mem_lookup <- membership$membership
  names(mem_lookup) <- as.character(membership$unit_id)

  transition_penalty_total <- 0
  transition_logLik_total <- 0
  lineage <- network$lineage
  if (nrow(lineage) > 0L) {
    lineage <- lineage[lineage$relation == "identity", , drop = FALSE]
    if (nrow(lineage) > 0L) {
      from_cluster <- unname(mem_lookup[as.character(lineage$from_unit)])
      to_cluster <- unname(mem_lookup[as.character(lineage$to_unit)])
      valid <- !is.na(from_cluster) & !is.na(to_cluster) & is.finite(lineage$weight) & lineage$weight > 0
      if (any(valid)) {
        probs <- transition$probabilities[
          cbind(as.character(from_cluster[valid]), as.character(to_cluster[valid]))
        ]
        transition_penalty_total <- sum(.transition_penalty(probs) * lineage$weight[valid])
        transition_logLik_total <- sum(log(probs) * lineage$weight[valid])
      }
    }
  }

  prior_penalty_total <- 0
  prior_logLik_total <- 0
  if (!is.null(prior) && !identical(prior$prior, "none")) {
    prior_probs <- .membership_prior_probs(prior)
    mem_prob <- prior_probs[as.character(membership$membership)]
    prior_penalty_total <- sum(.transition_penalty(mem_prob))
    prior_logLik_total <- sum(log(mem_prob))
  }

  deviance_total <- fits$deviance_total
  objective <- deviance_total + transition_penalty_total + prior_penalty_total
  logLik_total <- fits$logLik_total + transition_logLik_total + prior_logLik_total

  list(
    logLik = logLik_total,
    BIC = fits$BIC_total,
    ICL = fits$ICL_total,
    deviance = objective,
    deviance_total = deviance_total,
    transition_penalty_total = transition_penalty_total,
    prior_penalty_total = prior_penalty_total,
    objective = objective,
    objective_scale = "deviance",
    fit_logLik_total = fits$logLik_total,
    fit_deviance_total = deviance_total,
    transition_logLik_total = transition_logLik_total,
    prior_logLik_total = prior_logLik_total
  )
}

#' Initialize dynamic actor-time memberships
#'
#' Generate a one-based random initial partition for every active actor-time
#' unit. The helper is intentionally separate from optimization so that
#' initialization methods can be extended later without changing the
#' optimizer interface.
#'
#' @param network A `dynamic_network` object from [as_dynamic_network()].
#' @param k Positive fixed number of clusters.
#' @param method Initialization method. Currently only `"random"` is
#'   supported.
#'
#' @return A normalized actor-time membership data frame.
#' @export
initialize_dynamic_membership <- function(network, k, method = "random") {
  if (!inherits(network, "dynamic_network")) {
    stop("`network` must be a `dynamic_network` object.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k < 1 || k != as.integer(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)
  method <- match.arg(method, c("random"))
  n_units <- nrow(network$actor_time)
  if (n_units < k) {
    stop(
      "Cannot initialize ", k, " clusters with only ", n_units,
      " actor-time units; every cluster must be represented.",
      call. = FALSE
    )
  }

  membership <- if (identical(method, "random")) {
    sample(rep(seq_len(k), length.out = n_units))
  }

  out <- network$actor_time[, c(
    "unit_id", "actor_id", "time", "time_index", "row_index"
  ), drop = FALSE]
  out$membership <- as.integer(membership)
  out
}

#' Fit a dynamic GLM-Markov blockmodel with fixed K and multiple starts
#'
#' The optimizer uses an existing membership as the first start when supplied;
#' otherwise it generates random one-based initial memberships. Additional
#' starts are random. `seed = NULL` is the default and does not reset the RNG;
#' supply a numeric seed when reproducible starts are required.
#'
#' Deterministic R reference optimizer for the first dynamic GLM-Markov model.
#' The optimizer alternates between:
#' - fitting independent time-specific GLMs conditional on memberships;
#' - estimating Markov transitions and membership priors;
#' - scoring actor-time candidate clusters on the deviance scale;
#' - applying a conservative strict-improvement reassignment sweep.
#'
#' @param network A `dynamic_network` object from [as_dynamic_network()].
#' @param membership Initial actor-time memberships. Supported inputs follow
#'   the same conventions as [fit_time_glm_blockmodels()]. If `NULL`, random
#'   initialization requires a supplied fixed `k`.
#' @param k Optional fixed number of clusters. If omitted, it is inferred from
#'   the supplied initial memberships.
#' @param family GLM family specification. Supported values are `"binomial"`
#'   and `"ppml"`.
#' @param max_iter Maximum number of reassignment sweeps.
#' @param smoothing Additive smoothing used for Markov transitions and the
#'   empirical prior.
#' @param prior Membership prior mode: `"empirical"`, `"uniform"`, or `"none"`.
#' @param sweep_order Sweep order. This reference implementation supports only
#'   `"actor_time"`.
#' @param refit Refit strategy. This reference implementation supports only
#'   `"sweep"`.
#' @param verbose Logical; if `TRUE`, print progress messages.
#' @param tol Strict-improvement tolerance on the deviance scale.
#' @param n_starts Positive number of independent starts at fixed `k`.
#' @param seed Optional finite numeric seed. The default `NULL` uses the
#'   current RNG stream without resetting it. Explicit seeds are locally
#'   scoped and restore the caller's RNG state on return.
#'
#' @return A `dynamic_glm_blockmodel` object with final memberships, fitted
#'   time-specific GLM models, Markov transition estimates, objective history,
#'   and convergence diagnostics.
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
#'   t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))
#' dn <- as_dynamic_network(Y)
#' init <- data.frame(unit_id = dn$actor_time$unit_id, membership = 1L)
#' fit_dynamic_glm_blockmodel(dn, membership = init, k = 1, max_iter = 1)
#'
#' @export
fit_dynamic_glm_blockmodel <- function(network, membership = NULL, k = NULL,
                                       family = c("binomial", "ppml"),
                                       n_starts = 1L, seed = NULL,
                                       max_iter = 10,
                                       smoothing = 0.5,
                                       prior = c("empirical", "uniform", "none"),
                                       sweep_order = c("actor_time"),
                                       refit = c("sweep"),
                                       verbose = FALSE,
                                       tol = 1e-8) {
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
      stop("`seed` must be NULL or a single finite numeric value.", call. = FALSE)
    }
    had_state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_state) {
      old_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_state) {
        assign(".Random.seed", old_state, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  .fit_dynamic_glm_blockmodel_multiple(
    network = network,
    membership = membership,
    k = k,
    family = family,
    n_starts = n_starts,
    seed = seed,
    max_iter = max_iter,
    smoothing = smoothing,
    prior = prior,
    sweep_order = sweep_order,
    refit = refit,
    verbose = verbose,
    tol = tol
  )
}

#' @keywords internal
.fit_dynamic_glm_blockmodel_multiple <- function(network, membership, k,
                                                  family, n_starts, seed,
                                                  max_iter, smoothing, prior,
                                                  sweep_order, refit, verbose,
                                                  tol) {
  family <- force(family)
  prior <- force(prior)
  sweep_order <- force(sweep_order)
  refit <- force(refit)
  if (!inherits(network, "dynamic_network")) {
    stop("`network` must be a `dynamic_network` object.", call. = FALSE)
  }
  if (!is.numeric(n_starts) || length(n_starts) != 1L || is.na(n_starts) ||
      !is.finite(n_starts) || n_starts < 1 || n_starts != as.integer(n_starts)) {
    stop("`n_starts` must be a single positive integer.", call. = FALSE)
  }
  n_starts <- as.integer(n_starts)
  if (!is.numeric(max_iter) || length(max_iter) != 1L || is.na(max_iter) || max_iter < 0) {
    stop("`max_iter` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("`tol` must be a single non-negative numeric value.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single logical value.", call. = FALSE)
  }

  family_info <- glm_blockmodel_family(family)
  prior <- match.arg(prior, c("empirical", "uniform", "none"))
  sweep_order <- match.arg(sweep_order, c("actor_time"))
  refit <- match.arg(refit, c("sweep"))

  if (is.null(membership) && is.null(k)) {
    stop("`k` must be supplied when `membership` is NULL.", call. = FALSE)
  }
  if (!is.null(k) && (!is.numeric(k) || length(k) != 1L || is.na(k) ||
                      !is.finite(k) || k < 1 || k != as.integer(k))) {
    stop("`k` must be a single positive integer when supplied.", call. = FALSE)
  }
  if (!is.null(k)) {
    k <- as.integer(k)
  }

  membership_table <- if (is.null(membership)) {
    initialize_dynamic_membership(network, k = k)
  } else {
    .time_glm_normalize_membership(
      membership = membership,
      actor_time = network$actor_time,
      time_labels = network$times
    )
  }
  if (is.null(k)) {
    k <- max(membership_table$membership)
  } else if (k < max(membership_table$membership)) {
    stop("`k` cannot be smaller than the largest supplied membership label.", call. = FALSE)
  }

  start_memberships <- vector("list", n_starts)
  start_memberships[[1L]] <- membership_table
  if (n_starts > 1L) {
    for (start_idx in 2:n_starts) {
      start_memberships[[start_idx]] <- initialize_dynamic_membership(network, k = k)
    }
  }

  runs <- lapply(start_memberships, function(start_membership) {
    .fit_dynamic_glm_blockmodel_single(
      network = network,
      membership = start_membership,
      k = k,
      family = family,
      max_iter = max_iter,
      smoothing = smoothing,
      prior = prior,
      sweep_order = sweep_order,
      refit = refit,
      verbose = verbose,
      tol = tol
    )
  })
  start_objectives <- vapply(runs, function(x) x$objective, numeric(1))
  best_start <- which.min(start_objectives)
  best <- runs[[best_start]]
  compact_starts <- lapply(runs, function(x) {
    list(
      initial_membership = x$initial_membership,
      final_membership = x$membership,
      objective = x$objective,
      converged = x$converged,
      n_iter = x$n_iter,
      n_changes = x$n_changes,
      objective_history = x$objective_history,
      history = x$history
    )
  })

  best$n_starts <- n_starts
  best$best_start <- best_start
  best$start_objectives <- start_objectives
  best$start_converged <- vapply(runs, function(x) x$converged, logical(1))
  best$starts <- compact_starts
  best$seed <- seed
  best$control$n_starts <- n_starts
  best$control$seed <- seed
  best$call <- match.call()
  best
}

#' @keywords internal
.fit_dynamic_glm_blockmodel_single <- function(network, membership, k,
                                                family, max_iter, smoothing,
                                                prior, sweep_order, refit,
                                                verbose, tol) {
  family_info <- glm_blockmodel_family(family)
  membership_table <- .time_glm_normalize_membership(
    membership = membership,
    actor_time = network$actor_time,
    time_labels = network$times
  )
  initial_membership <- membership_table
  if (is.null(k)) {
    k <- max(membership_table$membership)
  } else {
    if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
        k < 1 || k != as.integer(k)) {
      stop("`k` must be a single positive integer when supplied.", call. = FALSE)
    }
    k <- as.integer(k)
    if (k < max(membership_table$membership)) {
      stop("`k` cannot be smaller than the largest supplied membership label.", call. = FALSE)
    }
  }
  candidate_clusters <- seq_len(k)

  if (verbose) {
    cat("Initial refit\n")
  }

  fits <- fit_time_glm_blockmodels(
    network,
    membership = membership_table,
    family = family_info
  )
  transition <- estimate_markov_transitions(
    network,
    membership = membership_table,
    k = k,
    smoothing = smoothing
  )
  prior_obj <- estimate_membership_prior(
    membership = membership_table$membership,
    prior = prior,
    k = k,
    smoothing = smoothing,
    labels = candidate_clusters
  )

  comps <- .dynamic_glm_objective_components(
    network = network,
    fits = fits,
    transition = transition,
    prior = prior_obj,
    membership = membership_table
  )

  history <- list(
    data.frame(
      iteration = 0L,
      n_changes = 0L,
      logLik = comps$logLik,
      BIC = comps$BIC,
      ICL = comps$ICL,
      deviance = comps$deviance,
      deviance_total = comps$deviance_total,
      transition_penalty_total = comps$transition_penalty_total,
      prior_penalty_total = comps$prior_penalty_total,
      objective = comps$objective,
      stringsAsFactors = FALSE
    )
  )

  total_changes <- 0L
  n_iter <- 0L
  converged <- FALSE

  if (max_iter > 0L) {
    for (iter in seq_len(as.integer(max_iter))) {
      if (verbose) {
        cat("Sweep ", iter, "\n", sep = "")
      }

      sweep_changes <- 0L
      for (row_idx in seq_len(nrow(network$actor_time))) {
        scores <- score_actor_time_candidates(
          x = network,
          fit = fits,
          membership = membership_table,
          row_index = row_idx,
          candidate_clusters = candidate_clusters,
          transition = transition,
          prior = prior_obj
        )

        current_cluster <- membership_table$membership[row_idx]
        current_score <- scores$total_score[scores$candidate_cluster == current_cluster][1L]
        best_idx <- which.min(scores$total_score)
        best_row <- scores[best_idx, , drop = FALSE]

        if (best_row$candidate_cluster != current_cluster &&
            is.finite(best_row$total_score) &&
            is.finite(current_score) &&
            best_row$total_score + tol < current_score) {
          membership_table$membership[row_idx] <- best_row$candidate_cluster
          sweep_changes <- sweep_changes + 1L
        }
      }

      total_changes <- total_changes + sweep_changes
      n_iter <- iter

      if (sweep_changes > 0L && identical(refit, "sweep")) {
        fits <- fit_time_glm_blockmodels(
          network,
          membership = membership_table,
          family = family_info
        )
        transition <- estimate_markov_transitions(
          network,
          membership = membership_table,
          k = k,
          smoothing = smoothing
        )
        prior_obj <- estimate_membership_prior(
          membership = membership_table$membership,
          prior = prior,
          k = k,
          smoothing = smoothing,
          labels = candidate_clusters
        )
      }

      comps <- .dynamic_glm_objective_components(
        network = network,
        fits = fits,
        transition = transition,
        prior = prior_obj,
        membership = membership_table
      )

      history[[length(history) + 1L]] <- data.frame(
        iteration = iter,
        n_changes = sweep_changes,
        logLik = comps$logLik,
        BIC = comps$BIC,
        ICL = comps$ICL,
        deviance = comps$deviance,
        deviance_total = comps$deviance_total,
        transition_penalty_total = comps$transition_penalty_total,
        prior_penalty_total = comps$prior_penalty_total,
        objective = comps$objective,
        stringsAsFactors = FALSE
      )

      if (verbose) {
        cat("  changes: ", sweep_changes, "\n", sep = "")
        cat("  objective: ", format(comps$objective), "\n", sep = "")
      }

      if (sweep_changes == 0L) {
        converged <- TRUE
        break
      }
    }
  }

  history <- do.call(rbind, history)
  rownames(history) <- NULL

  out <- list(
    membership = membership_table,
    initial_membership = initial_membership,
    k = k,
    family = family_info$name,
    pseudo = family_info$pseudo,
    fits = fits,
    transition = transition,
    prior = prior_obj,
    history = history,
    n_iter = n_iter,
    converged = converged,
    n_changes = total_changes,
    logLik = comps$logLik,
    BIC = comps$BIC,
    ICL = comps$ICL,
    deviance = comps$deviance,
    deviance_total = comps$deviance_total,
    transition_penalty_total = comps$transition_penalty_total,
    prior_penalty_total = comps$prior_penalty_total,
    objective = comps$objective,
    objective_history = history$objective,
    criterion_note = paste(
      family_info$criterion_note,
      "Objective is a deviance-scale reference criterion for fixed-K deterministic sweeps.",
      "Smaller values are better."
    ),
    control = list(
      max_iter = as.integer(max_iter),
      smoothing = smoothing,
      prior = prior,
      sweep_order = sweep_order,
      refit = refit,
      tol = tol
    ),
    call = match.call()
  )
  class(out) <- "dynamic_glm_blockmodel"
  out
}

#' Print a dynamic GLM-Markov optimizer fit
#'
#' @param x A `dynamic_glm_blockmodel` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.dynamic_glm_blockmodel <- function(x, ...) {
  cat("dynamic_glm_blockmodel object\n")
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  k: %s\n", format(x$k)))
  cat(sprintf("  n_iter: %s\n", format(x$n_iter)))
  cat(sprintf("  converged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  cat(sprintf("  logLik: %s\n", format(x$logLik)))
  cat(sprintf("  BIC: %s\n", format(x$BIC)))
  cat(sprintf("  ICL: %s\n", format(x$ICL)))
  cat(sprintf("  deviance: %s\n", format(x$deviance)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  cat(sprintf("  deviance_total: %s\n", format(x$deviance_total)))
  cat(sprintf("  transition_penalty_total: %s\n", format(x$transition_penalty_total)))
  cat(sprintf("  prior_penalty_total: %s\n", format(x$prior_penalty_total)))
  cat(sprintf("  changes: %s\n", format(x$n_changes)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a dynamic GLM-Markov optimizer fit
#'
#' @param object A `dynamic_glm_blockmodel` object.
#' @param ... Ignored.
#'
#' @return A compact summary list.
#' @export
summary.dynamic_glm_blockmodel <- function(object, ...) {
  summary <- list(
    k = object$k,
    family = object$family,
    pseudo = object$pseudo,
    n_iter = object$n_iter,
    converged = object$converged,
    n_changes = object$n_changes,
    logLik = object$logLik,
    BIC = object$BIC,
    ICL = object$ICL,
    deviance = object$deviance,
    objective = object$objective,
    objective_history = object$objective_history,
    deviance_total = object$deviance_total,
    transition_penalty_total = object$transition_penalty_total,
    prior_penalty_total = object$prior_penalty_total,
    criterion_note = object$criterion_note
  )
  class(summary) <- "summary.dynamic_glm_blockmodel"
  summary
}

#' @export
print.summary.dynamic_glm_blockmodel <- function(x, ...) {
  cat("Summary of dynamic_glm_blockmodel\n")
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  k: %s\n", format(x$k)))
  cat(sprintf("  n_iter: %s\n", format(x$n_iter)))
  cat(sprintf("  converged: %s\n", if (isTRUE(x$converged)) "TRUE" else "FALSE"))
  cat(sprintf("  changes: %s\n", format(x$n_changes)))
  cat(sprintf("  logLik: %s\n", format(x$logLik)))
  cat(sprintf("  BIC: %s\n", format(x$BIC)))
  cat(sprintf("  ICL: %s\n", format(x$ICL)))
  cat(sprintf("  deviance: %s\n", format(x$deviance)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  cat(sprintf("  deviance_total: %s\n", format(x$deviance_total)))
  cat(sprintf("  transition_penalty_total: %s\n", format(x$transition_penalty_total)))
  cat(sprintf("  prior_penalty_total: %s\n", format(x$prior_penalty_total)))
  cat(sprintf("  objective history length: %s\n", length(x$objective_history)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}
