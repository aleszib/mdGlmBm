# Testing Strategy

## Principles

Tests are part of the statistical specification. Every behavior that affects the model objective, labels, data representation, or result interpretation should be tested.

Every exported function must have testthat coverage. Every behavior change must add or update focused tests, and existing tests must continue to pass.

Use fast deterministic tests for CI. Put large examples and long comparisons outside routine CI.

## Initial test areas

### Package integrity

- package loads;
- exported functions exist;
- roxygen documentation builds;
- no top-level `library()` calls in package R files.

### Data representation

- list-of-matrices input validates;
- actors can enter;
- actors can leave;
- same actor ID in adjacent time points creates lineage;
- diagonal handling is explicit;
- missing dyads are not confused with observed zeros.

### Families and result metadata

- binomial mode sets `pseudo = FALSE`;
- PPML mode sets `pseudo = TRUE`;
- PPML result still has fields named `logLik`, `BIC`, `ICL`, `deviance`, and `objective`;
- print/summary notes pseudo-likelihood when `pseudo = TRUE`.

### Objective scale

- dynamic transition penalties are computed as `-2 * log(probability)` on the deviance scale;
- raw log-likelihood terms are not added directly to deviance terms;
- local score is smaller-is-better.
- transition smoothing prevents unintended infinite penalties.
- local scoring exposes and tests observation, predecessor, successor, and prior components.

### Optimizer behavior

- convergence and stopping criteria are tested once the optimizer exists;
- initialization and random starts are reproducible;
- actor entry has no predecessor penalty and actor exit has no successor penalty.

### Existing source code preservation

- existing `mdsbm` functions run on tiny examples after package consolidation;
- selected static GLM functions run on tiny examples.

### Future R/C++ parity

When C++ scoring is added, compare R and C++ scoring on fixed tiny inputs.

Equivalence tests should cover individual components, full candidate scores, membership updates, and objective history within documented numerical tolerances.

## CI expectations

GitHub Actions should run:

- `R CMD check` via `r-lib/actions`;
- package tests via `testthat`;
- optional coverage upload only if configured later.
