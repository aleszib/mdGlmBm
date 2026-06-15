#' Convert cluster labels to zero-based indexing
#'
#' @param x An integer vector or list of integer vectors.
#'
#' @return `x` with labels shifted down by one.
#' @keywords internal
bm_zero_index <- function(x) {
  if (is.list(x)) {
    return(lapply(x, bm_zero_index))
  }
  as.integer(x) - 1L
}

#' Convert cluster labels to one-based indexing
#'
#' @param x An integer vector or list of integer vectors.
#'
#' @return `x` with labels shifted up by one.
#' @keywords internal
bm_one_index <- function(x) {
  if (is.list(x)) {
    return(lapply(x, bm_one_index))
  }
  as.integer(x) + 1L
}

#' Extract a family name
#'
#' @param family A family object or family constructor.
#'
#' @return A single family name string when available.
#' @keywords internal
bm_family_name <- function(family) {
  if (inherits(family, "family") && !is.null(family$family)) {
    return(family$family)
  }
  if (is.function(family)) {
    family <- tryCatch(family(), error = function(...) NULL)
  }
  if (is.list(family) && !is.null(family$family)) {
    return(family$family)
  }
  as.character(family)[1]
}

#' Detect PPML-style families
#'
#' @param family A family object or family constructor.
#'
#' @return `TRUE` for `ppml()`-style families, otherwise `FALSE`.
#' @keywords internal
bm_is_ppml_family <- function(family) {
  if (identical(tolower(bm_family_name(family)), "ppml")) {
    return(TRUE)
  }
  FALSE
}

#' Build common result metadata
#'
#' @param fit A fitted `glm` object or `NULL`.
#' @param clu Optional hard partition.
#' @param ICL Optional ICL value to preserve.
#' @param pseudo Logical pseudo-likelihood flag.
#' @param objective Optional optimization objective.
#'
#' @return A list of common result fields.
#' @keywords internal
bm_result_metadata <- function(fit = NULL, clu = NULL, ICL = NULL,
                               pseudo = FALSE, objective = NULL) {
  logLik_value <- if (!is.null(fit)) as.numeric(stats::logLik(fit)) else NA_real_
  deviance_value <- if (!is.null(fit)) fit$deviance else NA_real_
  BIC_value <- if (!is.null(fit)) as.numeric(stats::BIC(fit)) else NA_real_
  if (is.null(objective)) {
    objective <- if (!is.null(ICL)) ICL else logLik_value
  }
  list(
    membership = clu,
    logLik = logLik_value,
    BIC = BIC_value,
    ICL = ICL,
    deviance = deviance_value,
    objective = objective,
    pseudo = isTRUE(pseudo)
  )
}

#' Compare partitions while ignoring label permutations
#'
#' @param a,clu2 Vector or list partitions to compare.
#'
#' @return `TRUE` when the partitions are equivalent up to relabeling.
#' @keywords internal
bm_partition_equal <- function(a, clu2) {
  if (is.list(a) || is.list(clu2)) {
    if (!is.list(a) || !is.list(clu2) || length(a) != length(clu2)) {
      return(FALSE)
    }
    return(isTRUE(all(vapply(seq_along(a), function(i) bm_partition_equal(a[[i]], clu2[[i]]), logical(1)))))
  }

  a <- as.integer(a)
  clu2 <- as.integer(clu2)
  if (length(a) != length(clu2)) {
    return(FALSE)
  }
  if (length(a) == 0L) {
    return(TRUE)
  }

  a_key <- match(a, unique(a))
  clu2_key <- match(clu2, unique(clu2))
  identical(a_key, clu2_key)
}

#' Generate a simple random partition
#'
#' @param n Number of units.
#' @param k Number of groups.
#' @param mingr,maxgr Ignored in the local fallback.
#' @param addParam Additional parameters ignored by the fallback.
#'
#' @return An integer vector of 1-based cluster labels.
#' @keywords internal
bm_random_partition <- function(n, k, mingr = NULL, maxgr = NULL, addParam = list()) {
  if (length(n) != 1L || length(k) != 1L) {
    stop("n and k must be scalar values.", call. = FALSE)
  }
  if (n <= 0L || k <= 0L) {
    stop("n and k must be positive.", call. = FALSE)
  }
  sample(rep(seq_len(k), length.out = n), size = n, replace = FALSE)
}

clu <- function(x) {
  if (is.list(x)) {
    if (!is.null(x$clu)) {
      return(x$clu)
    }
    if (!is.null(x$part)) {
      return(x$part)
    }
  }
  x
}

encodeToZeroIndexed <- function(x) {
  if (is.list(x)) {
    return(lapply(x, encodeToZeroIndexed))
  }
  ux <- sort(unique(x))
  match(x, ux) - 1L
}

numClust <- function(x) {
  if (is.list(x)) {
    return(vapply(x, function(y) length(unique(y)), integer(1)))
  }
  length(unique(x))
}

ss <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(0)
  }
  sum((x - mean(x))^2)
}

stochBlock <- function(M, clu, ...) {
  stop("stochBlock is not implemented in dynGLMbm yet.")
}
