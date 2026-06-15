#' Dynamic blockmodeling optimizer
#'
#' Runs the compiled dynamic blockmodeling optimizer. Cluster labels are
#' zero-based, matching the original Python implementation.
#'
#' @param x Square adjacency matrix.
#' @param sets Integer vector of time-set sizes.
#' @param k Integer vector with the number of ordinary clusters per set.
#' @param runs Number of random starts.
#' @param parallel Accepted for compatibility. The compiled implementation
#'   currently runs starts serially to keep R's random number generator safe.
#' @param symmetric Whether to use the symmetric likelihood.
#' @param verbose Print progress from compiled code.
#' @param epsilon Lower/upper clipping value for intra-time densities.
#' @param epsilonTrans Lower/upper clipping value for inter-time transition
#'   probabilities.
#' @param seed Optional R random seed.
#' @param ... Optimizer controls: mincluster, maxiter, chng_ratio, maxRuns,
#'   maxNoImp.
#' @return A list compatible with the old Python bridge: clu, sets, k, cluSize,
#'   ICL, and params. Additional diagnostics are also returned.
#' @export
mdsbm_opr <- function(x, sets, k, runs, parallel = FALSE, symmetric = TRUE,
                      verbose = FALSE, epsilon = 0.001,
                      epsilonTrans = 0.001, seed = NULL, ...) {
  controls <- bm_controls(...)
  x <- bm_integer_matrix(x)
  sets <- bm_integer_vector(sets, "sets")
  k <- bm_integer_vector(k, "k")
  runs <- as.integer(runs)
  bm_validate_dimensions(x, sets, k)
  if (!is.null(seed)) set.seed(seed)
  if (isTRUE(parallel)) {
    warning("parallel = TRUE is accepted for compatibility; mdsbm runs starts serially in C++.")
  }

  params <- as.list(environment())
  params$x <- NULL
  params$controls <- controls

  out <- bm_cpp_multiple_optimize(
    x, sets, k, runs, isTRUE(symmetric), isTRUE(verbose),
    epsilon, epsilonTrans, controls$mincluster, controls$maxiter,
    controls$chng_ratio, controls$maxRuns, controls$maxNoImp
  )
  out$params <- params
  out
}

#' Optimize one supplied dynamic blockmodeling partition
#'
#' @inheritParams mdsbm_opr
#' @param clu List of zero-based cluster membership vectors, one per set.
#' @export
mdsbm_one_partition <- function(x, sets, k, clu, symmetric = TRUE,
                                verbose = FALSE, epsilon = 0.001,
                                epsilonTrans = 0.001, seed = NULL, ...) {
  controls <- bm_controls(...)
  x <- bm_integer_matrix(x)
  sets <- bm_integer_vector(sets, "sets")
  k <- bm_integer_vector(k, "k")
  clu <- bm_partition_list(clu, sets, k)
  bm_validate_dimensions(x, sets, k)
  if (!is.null(seed)) set.seed(seed)

  params <- as.list(environment())
  params$x <- NULL
  params$controls <- controls

  out <- bm_cpp_optimize_partition(
    x, sets, k, clu, isTRUE(symmetric), isTRUE(verbose),
    epsilon, epsilonTrans, controls$mincluster, controls$maxiter,
    controls$chng_ratio, controls$maxRuns, controls$maxNoImp
  )
  out$params <- params
  out
}

#' Compute ICL for one supplied dynamic blockmodeling partition
#'
#' @inheritParams mdsbm_one_partition
#' @export
mdsbm_icl_one_partition <- function(x, sets, k, clu, symmetric = TRUE,
                                    epsilon = 0.001,
                                    epsilonTrans = 0.001) {
  x <- bm_integer_matrix(x)
  sets <- bm_integer_vector(sets, "sets")
  k <- bm_integer_vector(k, "k")
  clu <- bm_partition_list(clu, sets, k)
  bm_validate_dimensions(x, sets, k)

  params <- list(
    sets = sets, k = k, symmetric = symmetric,
    epsilon = epsilon, epsilonTrans = epsilonTrans
  )
  out <- bm_cpp_icl_partition(x, sets, k, clu, isTRUE(symmetric), epsilon, epsilonTrans)
  out$params <- params
  out
}

bm_controls <- function(...,
                        mincluster = 3L,
                        maxiter = 100L,
                        chng_ratio = 0.4,
                        maxRuns = 10L,
                        maxNoImp = 10L) {
  dots <- list(...)
  if (length(dots)) {
    unknown <- setdiff(names(dots), character())
    if (any(!nzchar(names(dots)))) {
      stop("All optimizer controls passed through ... must be named.", call. = FALSE)
    }
    formals_here <- names(formals(sys.function()))
    unknown <- setdiff(names(dots), formals_here)
    if (length(unknown)) {
      stop("Unknown optimizer control(s): ", paste(unknown, collapse = ", "), call. = FALSE)
    }
  }
  args <- utils::modifyList(
    list(
      mincluster = mincluster,
      maxiter = maxiter,
      chng_ratio = chng_ratio,
      maxRuns = maxRuns,
      maxNoImp = maxNoImp
    ),
    dots
  )
  args$mincluster <- as.integer(args$mincluster)
  args$maxiter <- as.integer(args$maxiter)
  args$maxRuns <- as.integer(args$maxRuns)
  args$maxNoImp <- as.integer(args$maxNoImp)
  args$chng_ratio <- as.numeric(args$chng_ratio)
  args
}

bm_integer_matrix <- function(x) {
  x <- as.matrix(x)
  if (length(dim(x)) != 2L || nrow(x) != ncol(x)) {
    stop("x must be a square matrix.", call. = FALSE)
  }
  storage.mode(x) <- "integer"
  x
}

bm_integer_vector <- function(x, name) {
  if (!is.numeric(x) && !is.integer(x)) {
    stop(name, " must be numeric/integer.", call. = FALSE)
  }
  x <- as.integer(x)
  if (anyNA(x) || any(x <= 0L)) {
    stop(name, " must contain positive integers.", call. = FALSE)
  }
  x
}

bm_partition_list <- function(clu, sets, k) {
  if (!is.list(clu) || length(clu) != length(sets)) {
    stop("clu must be a list with one vector per set.", call. = FALSE)
  }
  out <- vector("list", length(clu))
  for (i in seq_along(clu)) {
    zi <- as.integer(clu[[i]])
    if (length(zi) != sets[[i]]) {
      stop("clu[[", i, "]] must have length sets[[", i, "]].", call. = FALSE)
    }
    if (anyNA(zi) || any(zi < 0L) || any(zi >= k[[i]])) {
      stop("clu[[", i, "]] must use zero-based labels in 0:(k[[", i, "]] - 1).", call. = FALSE)
    }
    out[[i]] <- zi
  }
  out
}

bm_validate_dimensions <- function(x, sets, k) {
  if (length(sets) != length(k)) {
    stop("sets and k must have the same length.", call. = FALSE)
  }
  if (sum(sets) != nrow(x)) {
    stop("sum(sets) must equal nrow(x) and ncol(x).", call. = FALSE)
  }
  if (any(k > sets)) {
    stop("Each k value must be no larger than the corresponding set size.", call. = FALSE)
  }
  invisible(TRUE)
}
