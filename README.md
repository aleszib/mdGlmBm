# dynGLMbm

`dynGLMbm` is the package foundation for dynamic GLM-based blockmodeling.

Current status:

- the legacy `mdsbm` dynamic optimizer interface is preserved;
- cleaned static GLM blockmodeling helpers are imported and documented;
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
```

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
```

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

## Notes

- The exported `mdsbm_*` functions are preserved as the legacy dynamic API.
- Static GLM helpers are available for binomial and PPML-style use.
- The dynamic GLM-Markov optimizer is still pending and will be implemented in a later task.
