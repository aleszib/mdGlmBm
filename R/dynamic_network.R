#' Convert time-indexed matrices to a dynamic network object
#'
#' Build the package's internal dynamic data representation from a list of
#' square matrices with actor IDs in row and column names. Actors may enter or
#' leave across time points; only adjacent same-ID observations create lineage
#' links.
#'
#' @param Y A non-empty named or unnamed list of square matrix-like objects.
#' @param directed Logical scalar. If `TRUE` (default), dyads are treated as
#'   ordered pairs. If `FALSE`, dyads are stored once per unordered pair using
#'   the upper triangle of each matrix.
#' @param self_ties Logical scalar. If `FALSE` (default), diagonal entries are
#'   retained in the dyad table as unobserved rows with `value = NA` so diagonal
#'   handling is explicit. If `TRUE`, diagonal ties remain observed when present.
#' @param times Optional time labels. If `NULL`, names from `Y` are used when
#'   available; otherwise sequential labels `t1`, `t2`, ... are generated.
#'
#' @return An object of class `dynamic_network` with at least:
#' \describe{
#'   \item{times}{Time labels.}
#'   \item{actor_time}{Actor-time table.}
#'   \item{dyads}{Dyad table.}
#'   \item{lineage}{Identity lineage links between adjacent time points.}
#'   \item{directed}{Logical directedness flag.}
#'   \item{self_ties}{Logical self-tie handling flag.}
#' }
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
#'   t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))
#'
#' dn <- as_dynamic_network(Y)
#' dn$actor_time
#' dn$lineage
#'
#' @export
as_dynamic_network <- function(Y, directed = TRUE, self_ties = FALSE, times = NULL) {
  if (!is.list(Y) || length(Y) == 0L) {
    stop("`Y` must be a non-empty list of matrix-like objects.", call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("`directed` must be a single logical value.", call. = FALSE)
  }
  if (!is.logical(self_ties) || length(self_ties) != 1L || is.na(self_ties)) {
    stop("`self_ties` must be a single logical value.", call. = FALSE)
  }
  directed <- isTRUE(directed)
  self_ties <- isTRUE(self_ties)

  n_time <- length(Y)
  time_labels <- dynnet_validate_times(Y, times, n_time)
  matrices <- vector("list", n_time)
  actor_ids_by_time <- vector("list", n_time)
  dimensions <- vector("list", n_time)

  for (i in seq_len(n_time)) {
    validated <- dynnet_validate_matrix(Y[[i]], time_labels[i], i)
    matrices[[i]] <- validated$matrix
    actor_ids_by_time[[i]] <- validated$actor_ids
    dimensions[[i]] <- validated$dimension
  }

  actor_time <- dynnet_build_actor_time(actor_ids_by_time, time_labels)
  actor_lookup <- vector("list", n_time)
  for (i in seq_len(n_time)) {
    idx <- actor_time$time_index == i
    actor_lookup[[i]] <- setNames(actor_time$unit_id[idx], actor_time$actor_id[idx])
  }
  names(actor_lookup) <- as.character(seq_len(n_time))
  dyads <- dynnet_build_dyads(
    matrices = matrices,
    actor_ids_by_time = actor_ids_by_time,
    actor_lookup = actor_lookup,
    time_labels = time_labels,
    directed = directed,
    self_ties = self_ties
  )
  lineage <- dynnet_build_lineage(actor_lookup, actor_ids_by_time, time_labels)

  n_actors_by_time <- vapply(actor_ids_by_time, length, integer(1))
  n_dyads_by_time <- vapply(dyads, nrow, integer(1))
  dyads <- do.call(rbind, dyads)
  rownames(dyads) <- NULL

  structure(
    list(
      times = time_labels,
      matrices = matrices,
      actor_time = actor_time,
      dyads = dyads,
      lineage = lineage,
      directed = directed,
      self_ties = self_ties,
      dimensions = dimensions,
      actors_by_time = actor_ids_by_time,
      actors = sort(unique(unlist(actor_ids_by_time, use.names = FALSE))),
      n_actors_by_time = n_actors_by_time,
      n_dyads_by_time = n_dyads_by_time
    ),
    class = "dynamic_network"
  )
}

dynnet_validate_times <- function(Y, times, n_time) {
  if (is.null(times)) {
    nm <- names(Y)
    if (!is.null(nm) && all(nzchar(nm))) {
      time_labels <- as.character(nm)
    } else if (is.null(nm) || all(!nzchar(nm))) {
      time_labels <- paste0("t", seq_len(n_time))
    } else {
      stop(
        "If `times` is not supplied, list names must be either all present or all absent.",
        call. = FALSE
      )
    }
  } else {
    if (length(times) != n_time) {
      stop("`times` must have the same length as `Y`.", call. = FALSE)
    }
    time_labels <- as.character(times)
  }

  if (anyNA(time_labels) || any(!nzchar(time_labels))) {
    stop("Time labels must be non-missing and non-empty.", call. = FALSE)
  }
  if (anyDuplicated(time_labels)) {
    stop("Time labels must be unique.", call. = FALSE)
  }
  time_labels
}

dynnet_validate_matrix <- function(x, time_label, time_index) {
  if (is.null(dim(x)) || length(dim(x)) != 2L) {
    stop("Each element of `Y` must be a square matrix-like object.", call. = FALSE)
  }

  mat <- as.matrix(x)
  if (nrow(mat) != ncol(mat)) {
    stop("Each matrix in `Y` must be square.", call. = FALSE)
  }

  actor_ids <- rownames(mat)
  col_ids <- colnames(mat)
  if (is.null(actor_ids) || is.null(col_ids)) {
    stop("Each matrix in `Y` must have both row names and column names.", call. = FALSE)
  }
  if (anyNA(actor_ids) || anyNA(col_ids) || any(!nzchar(actor_ids)) || any(!nzchar(col_ids))) {
    stop("Row and column names must be non-missing and non-empty.", call. = FALSE)
  }
  if (!identical(actor_ids, col_ids)) {
    stop(
      "Row and column names must match exactly within each time point.",
      call. = FALSE
    )
  }
  if (anyDuplicated(actor_ids)) {
    stop("Row and column names must be unique within each time point.", call. = FALSE)
  }

  list(
    matrix = mat,
    actor_ids = as.character(actor_ids),
    dimension = c(nrow(mat), ncol(mat))
  )
}

dynnet_build_actor_time <- function(actor_ids_by_time, time_labels) {
  rows <- vector("list", length(actor_ids_by_time))
  unit_id <- 0L

  for (i in seq_along(actor_ids_by_time)) {
    actor_ids <- actor_ids_by_time[[i]]
    n_i <- length(actor_ids)
    if (n_i == 0L) {
      rows[[i]] <- data.frame(
        unit_id = integer(),
        actor_id = character(),
        time = character(),
        time_index = integer(),
        row_index = integer(),
        active = logical(),
        stringsAsFactors = FALSE
      )
      next
    }

    unit_ids <- seq.int(unit_id + 1L, length.out = n_i)
    unit_id <- unit_id + n_i
    rows[[i]] <- data.frame(
      unit_id = unit_ids,
      actor_id = actor_ids,
      time = rep.int(time_labels[i], n_i),
      time_index = rep.int(i, n_i),
      row_index = seq_len(n_i),
      active = rep.int(TRUE, n_i),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

dynnet_build_dyads <- function(matrices, actor_ids_by_time, actor_lookup, time_labels,
                               directed, self_ties) {
  rows <- vector("list", length(matrices))

  for (i in seq_along(matrices)) {
    mat <- matrices[[i]]
    actor_ids <- actor_ids_by_time[[i]]
    n_i <- length(actor_ids)

    if (n_i == 0L) {
      rows[[i]] <- data.frame(
        time = character(),
        time_index = integer(),
        sender = character(),
        receiver = character(),
        sender_unit = integer(),
        receiver_unit = integer(),
        value = numeric(),
        observed = logical(),
        dyad_id = character(),
        weight = numeric(),
        stringsAsFactors = FALSE
      )
      next
    }

    if (directed) {
      pair_idx <- expand.grid(
        sender = seq_len(n_i),
        receiver = seq_len(n_i),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
    } else {
      pair_idx <- which(upper.tri(mat, diag = TRUE), arr.ind = TRUE)
      pair_idx <- data.frame(sender = pair_idx[, 1], receiver = pair_idx[, 2])
    }

    sender_idx <- pair_idx$sender
    receiver_idx <- pair_idx$receiver
    values <- mat[cbind(sender_idx, receiver_idx)]
    observed <- !is.na(values)
    if (!self_ties) {
      diag_idx <- sender_idx == receiver_idx
      if (any(diag_idx)) {
        values[diag_idx] <- NA
        observed[diag_idx] <- FALSE
      }
    }

    sender_actor <- actor_ids[sender_idx]
    receiver_actor <- actor_ids[receiver_idx]
    sender_unit <- unname(actor_lookup[[as.character(i)]][sender_actor])
    receiver_unit <- unname(actor_lookup[[as.character(i)]][receiver_actor])

    rows[[i]] <- data.frame(
      time = rep.int(time_labels[i], length(sender_idx)),
      time_index = rep.int(i, length(sender_idx)),
      sender = sender_actor,
      receiver = receiver_actor,
      sender_unit = sender_unit,
      receiver_unit = receiver_unit,
      value = values,
      observed = observed,
      dyad_id = paste0(time_labels[i], ":", sender_actor, "->", receiver_actor),
      weight = rep.int(1, length(sender_idx)),
      stringsAsFactors = FALSE
    )
  }

  rows
}

dynnet_build_lineage <- function(actor_lookup, actor_ids_by_time, time_labels) {
  if (length(actor_lookup) <= 1L) {
    return(data.frame(
      from_unit = integer(),
      to_unit = integer(),
      from_time = character(),
      to_time = character(),
      from_time_index = integer(),
      to_time_index = integer(),
      actor_id = character(),
      relation = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  rows <- vector("list", length(actor_lookup) - 1L)
  edge_count <- 0L
  for (i in seq_len(length(actor_lookup) - 1L)) {
    current_lookup <- actor_lookup[[as.character(i)]]
    next_lookup <- actor_lookup[[as.character(i + 1L)]]
    common <- intersect(names(current_lookup), names(next_lookup))

    if (length(common) == 0L) {
      next
    }

    edge_count <- edge_count + 1L
    rows[[edge_count]] <- data.frame(
      from_unit = unname(current_lookup[common]),
      to_unit = unname(next_lookup[common]),
      from_time = rep.int(time_labels[i], length(common)),
      to_time = rep.int(time_labels[i + 1L], length(common)),
      from_time_index = rep.int(i, length(common)),
      to_time_index = rep.int(i + 1L, length(common)),
      actor_id = common,
      relation = rep.int("identity", length(common)),
      weight = rep.int(1, length(common)),
      stringsAsFactors = FALSE
    )
  }

  if (edge_count == 0L) {
    return(data.frame(
      from_unit = integer(),
      to_unit = integer(),
      from_time = character(),
      to_time = character(),
      from_time_index = integer(),
      to_time_index = integer(),
      actor_id = character(),
      relation = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows[seq_len(edge_count)])
  rownames(out) <- NULL
  out
}

#' Print a dynamic network object
#'
#' @param x A `dynamic_network` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.dynamic_network <- function(x, ...) {
  summary <- dynamic_network_summary(x)
  cat("dynamic_network object\n")
  cat(sprintf("  time points: %d\n", summary$n_timepoints))
  cat(sprintf("  unique actors: %d\n", summary$n_actors))
  cat(sprintf("  actor-time rows: %d\n", summary$n_actor_time))
  cat(sprintf("  dyad rows: %d\n", summary$n_dyads))
  cat(sprintf("  lineage edges: %d\n", summary$n_lineage))
  cat(sprintf("  directed: %s\n", if (summary$directed) "TRUE" else "FALSE"))
  cat(sprintf("  self ties: %s\n", if (summary$self_ties) "TRUE" else "FALSE"))
  cat("  times: ", paste(summary$times, collapse = ", "), "\n", sep = "")
  invisible(x)
}

#' Summarize a dynamic network object
#'
#' @param object A `dynamic_network` object.
#' @param ... Ignored.
#'
#' @return A list with compact summary information.
#' @export
summary.dynamic_network <- function(object, ...) {
  summary <- dynamic_network_summary(object)
  class(summary) <- "summary.dynamic_network"
  summary
}

#' @export
print.summary.dynamic_network <- function(x, ...) {
  cat("Summary of dynamic_network\n")
  cat(sprintf("  time points: %d\n", x$n_timepoints))
  cat(sprintf("  unique actors: %d\n", x$n_actors))
  cat(sprintf("  actor-time rows: %d\n", x$n_actor_time))
  cat(sprintf("  dyad rows: %d\n", x$n_dyads))
  cat(sprintf("  lineage edges: %d\n", x$n_lineage))
  cat(sprintf("  directed: %s\n", if (x$directed) "TRUE" else "FALSE"))
  cat(sprintf("  self ties: %s\n", if (x$self_ties) "TRUE" else "FALSE"))
  cat("  times: ", paste(x$times, collapse = ", "), "\n", sep = "")
  invisible(x)
}

dynamic_network_summary <- function(x) {
  if (!inherits(x, "dynamic_network")) {
    stop("`x` must be a `dynamic_network` object.", call. = FALSE)
  }

  list(
    times = x$times,
    n_timepoints = length(x$times),
    n_actors = length(x$actors),
    n_actor_time = nrow(x$actor_time),
    n_dyads = nrow(x$dyads),
    n_lineage = nrow(x$lineage),
    directed = isTRUE(x$directed),
    self_ties = isTRUE(x$self_ties),
    actors_by_time = setNames(x$n_actors_by_time, x$times),
    dyads_by_time = setNames(x$n_dyads_by_time, x$times)
  )
}
