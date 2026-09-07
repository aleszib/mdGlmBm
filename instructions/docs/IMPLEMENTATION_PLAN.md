# Implementation Plan

## Backend plan

### R reference backend

The first implementation should use R for:

- input validation;
- dynamic data object construction;
- time-specific GLM fitting;
- PPML/binomial family handling;
- local candidate scoring;
- objective accounting;
- deterministic tests.

This backend is the correctness baseline.

### C++ performance backend

After the R backend is correct, move these parts to C++/Rcpp:

- candidate scoring loops;
- actor-time membership update sweeps;
- transition count updates;
- sparse dyad traversal where applicable.

Keep GLM coefficient fitting in R initially.

## Phased implementation

1. **Package foundation:** establish the R package structure, namespace, roxygen2, testthat, CI, and preservation checks for existing `mdsbm` behavior.
2. **Dynamic data layer:** validate lists of square matrices, build `dynamic_network` with `actor_time`, `dyads`, and `lineage`, and handle actor entry/exit.
3. **Family/objective metadata:** define binomial and PPML handling, common result fields, and `pseudo = FALSE/TRUE` semantics.
4. **Time-specific GLM fitting:** fit independent per-time models with R GLM machinery and expose stable internal fit objects.
5. **Markov transition and local scoring:** implement smoothed transition/prior penalties on the `-2 * log(probability)` deviance/loss scale.
6. **First R reference optimizer:** implement deterministic mdsbm-style actor-time reassignment and objective accounting.
7. **Initialization/random starts:** add reproducible initialization, controlled random starts, and deterministic seed behavior.
8. **Diagnostics/objective history:** report convergence, stopping criteria, objective history, and pseudo-likelihood notes.
9. **C++ scoring backend:** add Rcpp scoring/update loops only after R reference equivalence tests are established.
10. **Model selection and later families:** add model-selection support and, after review, Poisson and Gaussian/normal families.
11. **Split/merge lineage support:** generalize lineage aggregation without blocking future unit split/merge behavior.

The current initialization layer uses a simple independent random partition for
fixed `K`. Explicit seeds are locally scoped; the default `seed = NULL` uses the
current RNG stream. Automatic `K` selection is not part of this phase.

## Early-task non-goals

Early tasks do not include broad refactoring, pooled/shared GLM coefficients across time, exact refitting for every candidate move as the default, full C++ GLM fitting, split/merge estimation, arbitrary C++ formula support, or CRAN-readiness claims. Preserve existing static GLM and `mdsbm` behavior unless a task explicitly changes it.

### Avoid full C++ GLM fitting initially

Do not reimplement GLM IRLS in C++ at the beginning. It would add risk around family handling, link functions, weights, offsets, convergence, separation, and rank deficiency.

## First public functions to aim for

These names are provisional and may be adjusted before public release:

```r
as_dynamic_network()
fit_static_glm_blockmodel()
fit_markov_glm_blockmodel()
score_actor_time_candidates()
compute_dynamic_objective()
```

## Result object requirements

All model results should include common fields:

```r
list(
  membership = ...,      # user-facing one-based labels
  time_membership = ...,
  family = ...,
  pseudo = TRUE/FALSE,
  logLik = ...,
  BIC = ...,
  ICL = ...,
  deviance = ...,
  objective = ...,
  transition = ...,
  glm = ...,
  call = ...,
  control = ...,
  diagnostics = ...
)
```

`print()` and `summary()` methods must show a PPML pseudo-likelihood note when `pseudo = TRUE`.

## Label conventions

User-facing labels should be one-based R labels. Existing zero-based labels from `mdsbm` internals must be converted at boundaries and tested.
