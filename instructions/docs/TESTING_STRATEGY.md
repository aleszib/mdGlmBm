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
- actor entry uses an emerging-to-substantive penalty and actor exit uses a
  substantive-to-vanishing penalty;
- exact augmented transition counts, smoothing, structural-zero behavior, and
  initial-state probabilities are tested;

### Existing source code preservation

- existing `mdsbm` functions run on tiny examples after package consolidation;
- selected static GLM functions run on tiny examples.

### Future R/C++ parity

When C++ scoring is added, compare R and C++ scoring on fixed tiny inputs.

Equivalence tests should cover individual components, full candidate scores, membership updates, and objective history within documented numerical tolerances.

The corrected reference semantics also require tests for persistent, entry,
and exit events, single-observation intermediate actors, and the distinction
between an observed exit and the end of the observation window.

### Simulation-based validation

Keep simulation validation separate from exact software invariants. Seeded,
optimizer-independent fixtures cover a strong binomial block structure,
multiple time points, actor entry/exit, and a small PPML count-data example.
Use ARI or another label-invariant measure for recovery. Also validate
multiple-start selection, global label permutation and actor-order invariance,
objective-component sums, local score-component sums, and convergence
diagnostics. Do not assume a fully recomputed global objective is monotone
unless the update/refit scheme proves it; record history and report observed
behavior instead. Keep larger simulation studies outside routine `testthat` and
`R CMD check`.

## CI expectations

GitHub Actions should run:

- `R CMD check` via `r-lib/actions`;
- package tests via `testthat`;
- optional coverage upload only if configured later.
