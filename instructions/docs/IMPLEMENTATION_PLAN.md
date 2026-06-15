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
