#' Balassa Normalization for a One-Mode Network
#'
#' Applies Balassa normalization to a matrix representing a one-mode network.
#' The input is first coerced to a matrix using [as.matrix()]. The normalization
#' divides each observed entry by its expected value under independence of row
#' and column totals.
#'
#' @param x A numeric object coercible to a matrix.
#'
#' @return A numeric matrix of the same dimensions as `x`, containing Balassa
#'   normalized values.
#'
#' @details
#' For an input matrix `x`, the expected value for cell `(i, j)` is computed as
#' `rowSums(x)[i] * colSums(x)[j] / sum(x)`. The returned matrix is `x /
#' expected`.
#'
#' @examples
#' m <- matrix(c(10, 5, 4, 8), nrow = 2, byrow = TRUE)
#' balassaNorm(m)
#'
#' @export
balassaNorm <- function(x) {
  x <- as.matrix(x)
  
  if (!is.numeric(x)) {
    stop("`x` must be numeric and coercible to a numeric matrix.")
  }
  
  total_exports <- rowSums(x, na.rm = TRUE)
  total_imports <- colSums(x, na.rm = TRUE)
  total_flow <- sum(x, na.rm = TRUE)
  
  if (total_flow == 0) {
    stop("The total sum of the matrix is 0, so Balassa normalization cannot be computed.")
  }
  
  expected <- total_exports %*% t(total_imports) / total_flow
  result <- x / expected
  
  return(result)
}

#' Compute ICL for Models Returned by `select.dynsbm`
#'
#' Recomputes the Integrated Classification Likelihood (ICL) values for the
#' models returned by the legacy `dynsbm` package's `select.dynsbm()` helper.
#' The implementation follows that package's internal `compute.icl()` function
#' used for plotting the ICL curve.
#'
#' @param x A list returned by the legacy `dynsbm` package's
#'   `select.dynsbm()` helper.
#'
#' @return A `data.frame` with one row per fitted model and the columns:
#'   \describe{
#'     \item{Q}{Number of groups.}
#'     \item{logLik}{Completed-data log-likelihood stored in the fitted model.}
#'     \item{ICL}{Recomputed ICL value.}
#'   }
#'
#' @details
#' The ICL is computed using the formula implemented in the `dynsbm` package:
#'
#' \deqn{
#' \mathrm{ICL} = \mathrm{loglikelihood}
#' - \mathbf{1}_{T > 1}\frac{1}{2}Q(Q - 1)N(T - 1)
#' - \mathrm{pen},
#' }
#'
#' where, in the binary case,
#'
#' \deqn{
#' \mathrm{pen} =
#' \frac{1}{2}Q\frac{N(N - 1)T}{2} +
#' \frac{1}{4}Q(Q - 1)T\frac{N(N - 1)}{2}.
#' }
#'
#' For continuous models, the penalty is doubled, matching the package source.
#' If a model contains a `gamma` component (discrete edge type), the original
#' package does not provide ICL values, so `NA` is returned.
#'
#' @examples
#' \dontrun{
#' library(dynsbm)
#' data(simdataT5Q4N40binary)
#' sel <- select.dynsbm(simdataT5Q4N40binary, Qmin = 1, Qmax = 6, nstart = 1)
#' dynsbmICL(sel)
#' }
#'
#' @references
#' dynsbm source code for `select.dynsbm()`, `draw.icl()` and `compute.icl()`:
#' [rdrr source](https://rdrr.io/cran/dynsbm/src/R/select.dynsbm.R)
#'
#' @export
dynsbmICL <- function(x) {
  if (!is.list(x) || length(x) == 0) {
    stop("`x` must be a non-empty list returned by `select.dynsbm()`.")
  }
  
  oneICL <- function(model) {
    required <- c("membership", "trans", "loglikelihood")
    missing_required <- setdiff(required, names(model))
    if (length(missing_required) > 0) {
      stop(
        "Each element of `x` must contain: ",
        paste(required, collapse = ", "),
        ". Missing: ",
        paste(missing_required, collapse = ", "),
        "."
      )
    }
    
    T <- ncol(model$membership)
    Q <- nrow(model$trans)
    N <- nrow(model$membership)
    logLik <- model$loglikelihood
    
    if ("gamma" %in% names(model)) {
      icl <- NA_real_
    } else {
      pen <- 0.5 * Q * (N * (N - 1) * T / 2) +
        0.25 * Q * (Q - 1) * T * (N * (N - 1) / 2)
      
      if ("sigma" %in% names(model)) {
        pen <- 2 * pen
      }
      
      trans_pen <- if (T > 1) 0.5 * Q * (Q - 1) * N * (T - 1) else 0
      icl <- logLik - trans_pen - pen
    }
    
    data.frame(
      Q = Q,
      logLik = logLik,
      ICL = icl,
      stringsAsFactors = FALSE
    )
  }
  
  res <- do.call(rbind, lapply(x, oneICL))
  rownames(res) <- NULL
  res
}
