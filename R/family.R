#' Normalize GLM blockmodel family choices
#'
#' Convert the package's supported family inputs into a structured metadata
#' object for static and future dynamic GLM blockmodeling.
#'
#' @param family A family specification. Supported values are `"binomial"`
#'   and `"ppml"`, or an equivalent family object.
#'
#' @return A list of class `glm_blockmodel_family` with elements:
#' \describe{
#'   \item{name}{Normalized family name.}
#'   \item{family}{A `stats::family` object ready for `glm()`.}
#'   \item{pseudo}{Logical pseudo-likelihood flag.}
#'   \item{objective_scale}{Scale used for model-selection objectives.}
#'   \item{local_score_scale}{Scale used for local reassignment scores.}
#'   \item{larger_is_better}{Logical flag for the objective direction.}
#'   \item{criterion_note}{Short text explaining how to interpret criteria.}
#' }
#'
#' @examples
#' glm_blockmodel_family("binomial")
#' glm_blockmodel_family("ppml")
#'
#' @export
glm_blockmodel_family <- function(family = c("binomial", "ppml")) {
  spec <- .glm_blockmodel_family_spec(family)
  structure(spec, class = "glm_blockmodel_family")
}

.glm_blockmodel_family <- glm_blockmodel_family

.glm_blockmodel_family_spec <- function(family) {
  if (inherits(family, "glm_blockmodel_family")) {
    return(unclass(family))
  }

  family_obj <- NULL
  if (inherits(family, "family")) {
    family_obj <- family
    family <- family$family
  } else if (is.function(family)) {
    family_obj <- tryCatch(family(), error = function(...) NULL)
    if (inherits(family_obj, "family")) {
      family <- family_obj$family
    } else {
      family <- deparse(substitute(family))
    }
  } else if (is.list(family) && !is.null(family$family)) {
    if (inherits(family$family, "family")) {
      family_obj <- family$family
      family <- family_obj$family
    } else {
      family <- family$family
    }
  }

  family_name <- tolower(as.character(family)[1])

  switch(
    family_name,
    binomial = {
      family_obj <- stats::binomial()
      list(
        name = "binomial",
        family = family_obj,
        pseudo = FALSE,
        objective_scale = "logLik",
        local_score_scale = "deviance",
        larger_is_better = TRUE,
        criterion_note = "Likelihood criteria are ordinary likelihood-based."
      )
    },
    ppml = {
      family_obj <- ppml()
      list(
        name = "ppml",
        family = family_obj,
        pseudo = TRUE,
        objective_scale = "logLik",
        local_score_scale = "deviance",
        larger_is_better = TRUE,
        criterion_note = "Likelihood criteria are pseudo-likelihood-based."
      )
    },
    if (!is.null(family_obj) && inherits(family_obj, "family")) {
      pseudo <- isTRUE(attr(family_obj, "pseudo"))
      list(
        name = family_name,
        family = family_obj,
        pseudo = pseudo,
        objective_scale = "logLik",
        local_score_scale = "deviance",
        larger_is_better = TRUE,
        criterion_note = if (pseudo) {
          "Likelihood criteria are pseudo-likelihood-based."
        } else {
          "Likelihood criteria follow the supplied family."
        }
      )
    } else {
      stop(
        "Unsupported family '", family_name, "'. Supported families are 'binomial' and 'ppml'.",
        call. = FALSE
      )
    }
  )
}

#' Normalize a GLM blockmodel result object
#'
#' Ensure a fitted result carries the common criterion fields used by the
#' package across binomial and PPML modes.
#'
#' @param result A list-like fit result.
#' @param fit Optional fitted `glm` object. If missing, `result$fit` is used
#'   when available.
#' @param family Optional family specification or normalized family object.
#' @param ICL Optional ICL value to preserve.
#' @param objective Optional objective value to preserve.
#' @param pseudo Optional pseudo-likelihood flag. If `NULL`, the flag is taken
#'   from `family` when available.
#' @param criterion_note Optional note about the reported criteria.
#'
#' @return `result` with normalized fields `logLik`, `BIC`, `ICL`,
#'   `deviance`, `objective`, `pseudo`, `family`, and `criterion_note`.
#' @keywords internal
normalize_glm_blockmodel_result <- function(result, fit = NULL, family = NULL,
                                            ICL = NULL, objective = NULL,
                                            pseudo = NULL, criterion_note = NULL) {
  if (is.null(result) || !is.list(result)) {
    stop("`result` must be a list-like object.", call. = FALSE)
  }

  if (is.null(fit) && !is.null(result$fit)) {
    fit <- result$fit
  }

  family_info <- NULL
  if (inherits(family, "glm_blockmodel_family")) {
    family_info <- family
  } else if (!is.null(family)) {
    family_info <- glm_blockmodel_family(family)
  } else if (!is.null(result$family)) {
    family_info <- tryCatch(glm_blockmodel_family(result$family), error = function(...) NULL)
  }

  if (is.null(pseudo)) {
    if (!is.null(family_info)) {
      pseudo <- family_info$pseudo
    } else if (!is.null(result$pseudo)) {
      pseudo <- result$pseudo
    } else {
      pseudo <- FALSE
    }
  }
  pseudo <- isTRUE(pseudo)

  if (is.null(criterion_note)) {
    if (!is.null(family_info)) {
      criterion_note <- family_info$criterion_note
    } else if (!is.null(result$criterion_note)) {
      criterion_note <- result$criterion_note
    }
  }

  if (is.null(fit) && !is.null(result$fit)) {
    fit <- result$fit
  }

  if (is.null(result$membership) && !is.null(result$clu)) {
    result$membership <- result$clu
  }

  if (is.null(result$logLik) && !is.null(fit)) {
    result$logLik <- as.numeric(stats::logLik(fit))
  }
  if (is.null(result$BIC) && !is.null(fit)) {
    result$BIC <- as.numeric(stats::BIC(fit))
  }
  if (is.null(result$deviance) && !is.null(fit)) {
    result$deviance <- fit$deviance
  }

  if (is.null(ICL) && !is.null(result$ICL)) {
    ICL <- result$ICL
  }
  result$ICL <- ICL

  if (is.null(objective)) {
    if (!is.null(result$objective)) {
      objective <- result$objective
    } else if (!is.null(ICL)) {
      objective <- ICL
    } else if (!is.null(result$logLik)) {
      objective <- result$logLik
    } else if (!is.null(fit)) {
      objective <- as.numeric(stats::logLik(fit))
    }
  }
  result$objective <- objective

  result$pseudo <- pseudo
  if (!is.null(family_info)) {
    result$family <- family_info$name
  } else if (!is.null(result$family)) {
    result$family <- bm_family_name(result$family)
  }
  if (is.null(result$criterion_note) && !is.null(criterion_note)) {
    result$criterion_note <- criterion_note
  }

  result
}
