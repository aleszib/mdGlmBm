#' Normalize actor-time memberships for time-specific GLM fitting
#'
#' @param membership Membership specification for actor-time units. Supported
#'   inputs are a vector aligned with `actor_time`, a data frame with `unit_id`
#'   and a membership column, or a list of per-time memberships.
#' @param actor_time The `actor_time` table from [as_dynamic_network()].
#' @param time_labels Character vector of time labels.
#'
#' @return A data frame with `unit_id`, `actor_id`, `time`, `time_index`,
#'   `row_index`, and `membership`.
#' @keywords internal
.time_glm_normalize_membership <- function(membership, actor_time, time_labels) {
  if (!is.data.frame(actor_time)) {
    stop("`actor_time` must be a data frame.", call. = FALSE)
  }

  if (is.list(membership) && !is.data.frame(membership)) {
    if (length(membership) != length(time_labels)) {
      stop("Per-time membership lists must match the number of time points.", call. = FALSE)
    }

    pieces <- vector("list", length(time_labels))
    for (i in seq_along(time_labels)) {
      idx <- actor_time$time_index == i
      pieces[[i]] <- .time_glm_normalize_single_membership(
        membership = membership[[i]],
        actor_time = actor_time[idx, , drop = FALSE],
        time_label = time_labels[i]
      )
    }
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    return(out)
  }

  .time_glm_normalize_single_membership(
    membership = membership,
    actor_time = actor_time,
    time_label = paste(time_labels, collapse = ", ")
  )
}

#' @keywords internal
.time_glm_normalize_single_membership <- function(membership, actor_time, time_label) {
  n <- nrow(actor_time)

  if (is.data.frame(membership)) {
    col_unit <- if ("unit_id" %in% names(membership)) "unit_id" else NULL
    col_membership <- intersect(c("membership", "cluster", "clu", "label"), names(membership))
    col_membership <- if (length(col_membership)) col_membership[[1]] else NULL

    if (is.null(col_membership)) {
      stop(
        "Membership data frames must contain a membership column named one of: membership, cluster, clu, label.",
        call. = FALSE
      )
    }

    if (is.null(col_unit)) {
      if (nrow(membership) != n) {
        stop("Membership data frames must either include `unit_id` or match the actor-time row count.", call. = FALSE)
      }
      values <- membership[[col_membership]]
    } else {
      unit_match <- match(actor_time$unit_id, membership[[col_unit]])
      if (anyNA(unit_match)) {
        stop(
          "Membership data are missing rows for one or more actor-time units.",
          call. = FALSE
        )
      }
      values <- membership[[col_membership]][unit_match]
    }
  } else {
    values <- membership
    if (is.null(values)) {
      stop("Memberships must be supplied for all active actor-time units.", call. = FALSE)
    }

    if (length(values) != n) {
      stop(
        "Membership vectors must have one value per actor-time unit.",
        call. = FALSE
      )
    }

    nm <- names(values)
    if (!is.null(nm)) {
      nm <- as.character(nm)
      if (anyDuplicated(nm)) {
        stop("Membership names must be unique when supplied.", call. = FALSE)
      }
      if (setequal(nm, as.character(actor_time$unit_id))) {
        values <- values[match(as.character(actor_time$unit_id), nm)]
      } else if (nrow(actor_time) == length(unique(actor_time$actor_id)) && setequal(nm, actor_time$actor_id)) {
        values <- values[match(actor_time$actor_id, nm)]
      } else {
        stop(
          "Membership names do not align with the actor-time units for time point '",
          time_label, "'.",
          call. = FALSE
        )
      }
    }
  }

  values <- .time_glm_coerce_memberships(values, time_label)

  data.frame(
    unit_id = actor_time$unit_id,
    actor_id = actor_time$actor_id,
    time = actor_time$time,
    time_index = actor_time$time_index,
    row_index = actor_time$row_index,
    membership = values,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.time_glm_coerce_memberships <- function(x, time_label) {
  if (is.factor(x)) {
    x <- as.integer(x)
  } else if (is.character(x)) {
    x <- as.integer(factor(x))
  } else {
    x <- suppressWarnings(as.numeric(x))
  }

  if (anyNA(x)) {
    stop(
      "Memberships must be positive integer or factor-like labels for time point '",
      time_label, "'.",
      call. = FALSE
    )
  }

  if (any(abs(x - round(x)) > sqrt(.Machine$double.eps))) {
    stop(
      "Memberships must be positive integer or factor-like labels for time point '",
      time_label, "'.",
      call. = FALSE
    )
  }

  x <- as.integer(round(x))
  if (any(x < 1L)) {
    stop(
      "Membership labels must be one-based positive integers for time point '",
      time_label, "'.",
      call. = FALSE
    )
  }

  x
}

#' @keywords internal
.time_glm_time_index <- function(x, time = NULL) {
  if (!inherits(x, "dynamic_network")) {
    stop("`x` must be a `dynamic_network` object.", call. = FALSE)
  }

  if (is.null(time)) {
    if (length(x$times) != 1L) {
      stop("`time` must be supplied when fitting one time point from a multi-time network.", call. = FALSE)
    }
    return(1L)
  }

  if (is.numeric(time) && length(time) == 1L) {
    time_index <- as.integer(time)
    if (is.na(time_index) || time_index < 1L || time_index > length(x$times)) {
      stop("`time` is out of range.", call. = FALSE)
    }
    return(time_index)
  }

  time <- as.character(time)[1]
  idx <- match(time, x$times)
  if (is.na(idx)) {
    stop("`time` was not found in the dynamic network object.", call. = FALSE)
  }
  idx
}

#' @keywords internal
.time_glm_build_data <- function(x, membership_table, time_index) {
  dyads_t <- x$dyads[x$dyads$time_index == time_index & x$dyads$observed, , drop = FALSE]
  if (nrow(dyads_t) == 0L) {
    stop("No observed dyads are available for the selected time point.", call. = FALSE)
  }

  mem_lookup <- membership_table$membership
  names(mem_lookup) <- as.character(membership_table$unit_id)
  sender_membership <- unname(mem_lookup[as.character(dyads_t$sender_unit)])
  receiver_membership <- unname(mem_lookup[as.character(dyads_t$receiver_unit)])

  if (anyNA(sender_membership) || anyNA(receiver_membership)) {
    stop("Memberships must be supplied for all active actor-time units referenced by the dyads.", call. = FALSE)
  }

  dyads_t$sender_membership <- sender_membership
  dyads_t$receiver_membership <- receiver_membership
  dyads_t$block <- paste(sender_membership, receiver_membership, sep = "#")
  dyads_t$block <- factor(dyads_t$block)
  dyads_t
}

#' Fit one time-specific GLM blockmodel
#'
#' Fit the observation model for a single time point of a `dynamic_network`
#' object conditional on supplied actor-time memberships.
#'
#' @param x A `dynamic_network` object from [as_dynamic_network()].
#' @param membership Membership labels for actor-time units. Supported inputs
#'   are a full actor-time vector, a matching data frame, or a list of
#'   per-time memberships.
#' @param time The time point to fit. For multi-time networks this may be a
#'   time label or one-based index. If omitted, `x` must contain exactly one
#'   time point.
#' @param family GLM family specification. Supported values are `"binomial"`
#'   and `"ppml"`.
#' @param formula Optional model formula. The default is `value ~ block`.
#'
#' @return A `time_glm_blockmodel` object with the fitted `glm` model, common
#'   criterion fields, and normalized membership data for the selected time.
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1, 1, 0), nrow = 2, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dn <- as_dynamic_network(Y)
#' fit_time_glm_blockmodel(dn, membership = c(1L, 1L), family = "binomial")
#'
#' @export
fit_time_glm_blockmodel <- function(x, membership, time = NULL,
                                    family = c("binomial", "ppml"),
                                    formula = value ~ block) {
  family_info <- glm_blockmodel_family(family)
  time_index <- .time_glm_time_index(x, time)
  time_label <- x$times[time_index]
  membership_table <- .time_glm_normalize_membership(
    membership = membership,
    actor_time = x$actor_time,
    time_labels = x$times
  )
  membership_time <- membership_table[membership_table$time_index == time_index, , drop = FALSE]
  if (nrow(membership_time) == 0L) {
    stop("No memberships were found for the selected time point.", call. = FALSE)
  }

  fit_data <- .time_glm_build_data(x, membership_table, time_index)
  if (length(unique(fit_data$block)) <= 1L) {
    fit_formula <- stats::as.formula("value ~ 1")
  } else if (is.null(formula)) {
    fit_formula <- stats::as.formula("value ~ block")
  } else if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.", call. = FALSE)
  } else {
    fit_formula <- formula
    if (!grepl("\\bblock\\b", paste(deparse(fit_formula), collapse = " "))) {
      fit_formula <- stats::update(fit_formula, . ~ . + block)
    }
  }

  fit <- stats::glm(fit_formula, data = fit_data, family = family_info$family)
  logLik_value <- as.numeric(stats::logLik(fit))
  bic_value <- as.numeric(stats::BIC(fit))
  icl_value <- NA_real_
  objective_value <- logLik_value
  criterion_note <- paste(
    family_info$criterion_note,
    "Objective is the fitted log-likelihood for the time-specific observation model.",
    "ICL is not computed for fixed-membership time-specific fits."
  )

  out <- normalize_glm_blockmodel_result(
    list(
      time = time_label,
      time_index = time_index,
      fit = fit,
      data = fit_data,
      membership = membership_time,
      time_membership = setNames(list(setNames(membership_time$membership, membership_time$actor_id)), time_label),
      formula = fit_formula,
      n_dyads = nrow(fit_data),
      n_memberships = nrow(membership_time)
    ),
    fit = fit,
    family = family_info,
    ICL = icl_value,
    objective = objective_value,
    pseudo = family_info$pseudo,
    criterion_note = criterion_note
  )

  out$family_object <- family_info$family
  out$objective_scale <- family_info$objective_scale
  out$local_score_scale <- family_info$local_score_scale
  out$larger_is_better <- family_info$larger_is_better
  out$time <- time_label
  out$time_index <- time_index
  out$time_membership <- setNames(list(setNames(membership_time$membership, membership_time$actor_id)), time_label)
  class(out) <- c("time_glm_blockmodel", "list")
  out
}

#' Fit time-specific GLM blockmodels across all time points
#'
#' Fit the observation model independently at each time point of a
#' `dynamic_network` object using supplied actor-time memberships.
#'
#' @param x A `dynamic_network` object from [as_dynamic_network()].
#' @param membership Membership labels for actor-time units. See
#'   `fit_time_glm_blockmodel()` for accepted formats.
#' @param family GLM family specification. Supported values are `"binomial"`
#'   and `"ppml"`.
#' @param formula Optional model formula. The default is `value ~ block`.
#'
#' @return A `time_glm_blockmodels` object with per-time fits and aggregated
#'   criterion summaries.
#'
#' @examples
#' Y <- list(
#'   t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
#'   t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
#' )
#' dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
#' dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))
#' dn <- as_dynamic_network(Y)
#' mem <- c(1L, 1L, 1L, 1L, 2L)
#' fit_time_glm_blockmodels(dn, membership = mem, family = "binomial")
#'
#' @export
fit_time_glm_blockmodels <- function(x, membership,
                                     family = c("binomial", "ppml"),
                                     formula = value ~ block) {
  family_info <- glm_blockmodel_family(family)
  membership_table <- .time_glm_normalize_membership(
    membership = membership,
    actor_time = x$actor_time,
    time_labels = x$times
  )

  fits <- lapply(seq_along(x$times), function(i) {
    fit_time_glm_blockmodel(
      x = x,
      membership = membership_table,
      time = i,
      family = family_info,
      formula = formula
    )
  })
  names(fits) <- x$times

  logLik_vec <- vapply(fits, function(z) z$logLik, numeric(1))
  deviance_vec <- vapply(fits, function(z) z$deviance, numeric(1))
  bic_vec <- vapply(fits, function(z) z$BIC, numeric(1))
  icl_vec <- vapply(fits, function(z) z$ICL, numeric(1))
  objective_vec <- vapply(fits, function(z) z$objective, numeric(1))
  total_logLik <- sum(logLik_vec)
  total_deviance <- sum(deviance_vec)
  total_bic <- if (all(is.finite(bic_vec))) sum(bic_vec) else NA_real_
  total_icl <- if (all(is.finite(icl_vec))) sum(icl_vec) else NA_real_

  criterion_note <- paste(
    family_info$criterion_note,
    "Per-time fits are independent observation models.",
    "ICL is not computed for the fixed-membership time-specific layer."
  )

  out <- list(
    fits = fits,
    family = family_info$name,
    family_object = family_info$family,
    pseudo = family_info$pseudo,
    criterion_note = criterion_note,
    membership = membership_table,
    time_membership = setNames(
      lapply(seq_along(x$times), function(i) {
        idx <- membership_table$time_index == i
        setNames(membership_table$membership[idx], membership_table$actor_id[idx])
      }),
      x$times
    ),
    logLik = logLik_vec,
    logLik_total = total_logLik,
    deviance = deviance_vec,
    deviance_total = total_deviance,
    BIC = bic_vec,
    BIC_total = total_bic,
    ICL = icl_vec,
    ICL_total = total_icl,
    objective = total_logLik,
    objective_by_time = objective_vec,
    objective_scale = family_info$objective_scale,
    local_score_scale = family_info$local_score_scale,
    larger_is_better = family_info$larger_is_better,
    formula = formula,
    n_timepoints = length(x$times),
    time_labels = x$times,
    call = match.call()
  )
  class(out) <- "time_glm_blockmodels"
  out
}

#' Print a time-specific GLM blockmodel fit
#'
#' @param x A `time_glm_blockmodel` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.time_glm_blockmodel <- function(x, ...) {
  cat("time_glm_blockmodel object\n")
  cat(sprintf("  time: %s\n", x$time))
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  logLik: %s\n", format(x$logLik)))
  cat(sprintf("  deviance: %s\n", format(x$deviance)))
  cat(sprintf("  BIC: %s\n", format(x$BIC)))
  cat(sprintf("  ICL: %s\n", format(x$ICL)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a time-specific GLM blockmodel fit
#'
#' @param object A `time_glm_blockmodel` object.
#' @param ... Ignored.
#'
#' @return A compact summary list.
#' @export
summary.time_glm_blockmodel <- function(object, ...) {
  summary <- list(
    time = object$time,
    family = object$family,
    pseudo = object$pseudo,
    logLik = object$logLik,
    deviance = object$deviance,
    BIC = object$BIC,
    ICL = object$ICL,
    objective = object$objective,
    n_dyads = object$n_dyads,
    n_memberships = object$n_memberships,
    criterion_note = object$criterion_note
  )
  class(summary) <- "summary.time_glm_blockmodel"
  summary
}

#' @export
print.summary.time_glm_blockmodel <- function(x, ...) {
  cat("Summary of time_glm_blockmodel\n")
  cat(sprintf("  time: %s\n", x$time))
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  logLik: %s\n", format(x$logLik)))
  cat(sprintf("  deviance: %s\n", format(x$deviance)))
  cat(sprintf("  BIC: %s\n", format(x$BIC)))
  cat(sprintf("  ICL: %s\n", format(x$ICL)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  cat(sprintf("  observed dyads: %s\n", format(x$n_dyads)))
  cat(sprintf("  memberships: %s\n", format(x$n_memberships)))
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Print a multi-time GLM blockmodel fit
#'
#' @param x A `time_glm_blockmodels` object.
#' @param ... Ignored.
#'
#' @return Invisibly returns `x`.
#' @export
print.time_glm_blockmodels <- function(x, ...) {
  cat("time_glm_blockmodels object\n")
  cat(sprintf("  time points: %d\n", length(x$time_labels)))
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  total logLik: %s\n", format(x$logLik_total)))
  cat(sprintf("  total deviance: %s\n", format(x$deviance_total)))
  cat(sprintf("  total BIC: %s\n", format(x$BIC_total)))
  cat(sprintf("  total ICL: %s\n", format(x$ICL_total)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  cat("  times: ", paste(x$time_labels, collapse = ", "), "\n", sep = "")
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}

#' Summarize a multi-time GLM blockmodel fit
#'
#' @param object A `time_glm_blockmodels` object.
#' @param ... Ignored.
#'
#' @return A compact summary list.
#' @export
summary.time_glm_blockmodels <- function(object, ...) {
  summary <- list(
    times = object$time_labels,
    n_timepoints = length(object$time_labels),
    family = object$family,
    pseudo = object$pseudo,
    logLik_total = object$logLik_total,
    deviance_total = object$deviance_total,
    BIC_total = object$BIC_total,
    ICL_total = object$ICL_total,
    objective = object$objective,
    criterion_note = object$criterion_note
  )
  class(summary) <- "summary.time_glm_blockmodels"
  summary
}

#' @export
print.summary.time_glm_blockmodels <- function(x, ...) {
  cat("Summary of time_glm_blockmodels\n")
  cat(sprintf("  time points: %d\n", x$n_timepoints))
  cat(sprintf("  family: %s\n", x$family))
  cat(sprintf("  pseudo: %s\n", if (isTRUE(x$pseudo)) "TRUE" else "FALSE"))
  cat(sprintf("  total logLik: %s\n", format(x$logLik_total)))
  cat(sprintf("  total deviance: %s\n", format(x$deviance_total)))
  cat(sprintf("  total BIC: %s\n", format(x$BIC_total)))
  cat(sprintf("  total ICL: %s\n", format(x$ICL_total)))
  cat(sprintf("  objective: %s\n", format(x$objective)))
  cat("  times: ", paste(x$times, collapse = ", "), "\n", sep = "")
  if (!is.null(x$criterion_note)) {
    cat("  note: ", x$criterion_note, "\n", sep = "")
  }
  invisible(x)
}
