# dynGLMbm

`dynGLMbm` is the package foundation for dynamic GLM-based blockmodeling.

Current status:

- the legacy `mdsbm` dynamic optimizer interface is preserved;
- cleaned static GLM blockmodeling helpers are imported and documented;
- a first dynamic data-layer constructor is available for time-indexed network input;
- time-specific GLM observation-model fitting is available for fixed memberships;
- the dynamic GLM-Markov optimizer is not implemented yet.

## Static GLM example

```r
library(dynGLMbm)

M <- matrix(c(0, 1,
              1, 0), nrow = 2, byrow = TRUE)
dimnames(M) <- list(c("a", "b"), c("a", "b"))

res <- optParGlm(M, clu = c(1L, 1L), maxIter = 0, glmFamily = binomial)
res$ICL
res$pseudo
res$criterion_note
```

`binomial` is the true-likelihood mode.

PPML uses the same result field names, but sets `pseudo = TRUE` and reports a
note that the likelihood-style criteria are pseudo-likelihood-based.

PPML is available through the `ppml()` family:

```r
res_ppml <- optParGlm(
  matrix(c(0, 1,
           2, 0), nrow = 2, byrow = TRUE),
  clu = c(1L, 1L),
  maxIter = 0,
  glmFamily = ppml()
)
res_ppml$pseudo
res_ppml$criterion_note
```

## Time-specific GLM example

```r
Y <- list(
  t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
  t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
)
dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))

dn <- as_dynamic_network(Y)

membership <- data.frame(
  unit_id = dn$actor_time$unit_id,
  membership = c(1L, 1L, 1L, 1L, 2L)
)

fit_time_glm_blockmodels(dn, membership = membership, family = "binomial")

fit_time_glm_blockmodels(dn, membership = membership, family = "ppml")
```

## Reference scoring layer

```r
transitions <- estimate_markov_transitions(dn, membership = membership)
prior <- estimate_membership_prior(membership$membership, prior = "empirical")

score_actor_time_candidates(
  dn,
  fit_time_glm_blockmodels(dn, membership = membership, family = "binomial"),
  membership = membership,
  row_index = 1L,
  transition = transitions,
  prior = prior
)
```

`score_actor_time_candidates()` reports deviance-scale local scores for
candidate memberships. The dynamic GLM-Markov optimizer is still not
implemented, so this is a reference scoring layer only.

## Reference optimizer

```r
Y <- list(
  t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
  t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
)
dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))

dn <- as_dynamic_network(Y)
init <- data.frame(unit_id = dn$actor_time$unit_id, membership = 1L)

fit <- fit_dynamic_glm_blockmodel(dn, membership = init, k = 1, max_iter = 1)
fit$membership
fit$objective_history
```

This is the first R reference dynamic optimizer. It requires a supplied initial
membership table and fixed `K`; it does not yet do random starts, C++ scoring,
split/merge, or automatic model selection.

## Legacy dynamic example

```r
x <- matrix(c(0L, 1L,
              1L, 0L), nrow = 2, byrow = TRUE)
dimnames(x) <- list(c("a", "b"), c("a", "b"))

mdsbm_icl_one_partition(
  x = x,
  sets = 2L,
  k = 1L,
  clu = list(c(0L, 0L))
)
```

## Dynamic data layer example

```r
Y <- list(
  t1 = matrix(c(0, 1, 0, 0), nrow = 2, byrow = TRUE),
  t2 = matrix(c(0, 1, 1, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = TRUE)
)

dimnames(Y$t1) <- list(c("A", "B"), c("A", "B"))
dimnames(Y$t2) <- list(c("A", "B", "C"), c("A", "B", "C"))

dn <- as_dynamic_network(Y)
dn$actor_time
dn$lineage
```

## Notes

- The exported `mdsbm_*` functions are preserved as the legacy dynamic API.
- Static GLM helpers are available for binomial and PPML-style use.
- `as_dynamic_network()` is the first internal/public data representation step
  for actor-time and lineage handling.
- `fit_time_glm_blockmodel()` and `fit_time_glm_blockmodels()` fit independent
  observation models at fixed memberships.
- `estimate_markov_transitions()`, `estimate_membership_prior()`, and
  `score_actor_time_candidates()` provide the reference dynamic scoring layer.
- `fit_dynamic_glm_blockmodel()` is the first R reference dynamic optimizer
  for fixed `K` and supplied starting memberships.
- The dynamic GLM-Markov optimizer is still pending and will be implemented in a later task.
