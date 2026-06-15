# AGENTS.md

## Project mission

This repository is being developed into an R package for dynamic GLM-based blockmodeling. The long-term model should preserve the dynamic optimization logic of `mdsbm`, while replacing the static Bernoulli/block-density observation component with GLM-based blockmodeling.

The target package name is provisional: `dynGLMbm`. Rename only if explicitly instructed.

## Human and agent roles

The human lead owns the statistical meaning, research scope, acceptance criteria, and release decisions. The strategic AI owns architecture guidance, work-order design, evidence review, and long-term project memory. The execution agent owns bounded implementation tasks, tests, documentation updates, commits, and pull requests.

Do not expand scope beyond the current task. Do not infer new statistical methodology without documenting it and asking for review.

## Non-negotiable statistical doctrine

1. Keep the `mdsbm`-style dynamic optimization logic as the long-term target.
2. The first dynamic GLM-Markov estimator should use independent time-specific GLM observation models.
3. Temporal coupling enters through the latent membership process, not through shared GLM coefficients in the first implementation.
4. Candidate actor-time cluster assignments are scored on the deviance/loss scale:
   - individual GLM deviance contribution at time `t`;
   - plus `-2 * log(previous transition probability)` when a predecessor exists;
   - plus `-2 * log(next transition probability)` when a successor exists;
   - plus `-2 * log(membership prior)` if a prior term is used.
5. Raw log-likelihood terms must never be added directly to raw deviance terms. Convert dynamic log-probability terms to deviance scale by multiplying by `-2`.
6. Final reporting computes a total log-likelihood-style criterion after full iterations, but local reassignment uses individual deviance-style scores.
7. The first supported dynamic families are `binomial` and `ppml`.
8. PPML uses the same result field names as other families (`logLik`, `BIC`, `ICL`, `deviance`, `objective`) but must set `pseudo = TRUE` and print a clear note that likelihood criteria are pseudo-likelihood-based.
9. True likelihood families set `pseudo = FALSE`.
10. Actors may enter and leave over time. Do not assume a fixed actor set across all time points.
11. Split and merge of units over time are not first-release requirements, but the internal lineage design must not block them.

## Data representation doctrine

The first public data interface should be a list of time-specific square networks with actor IDs in row and column names. Each time point may have a different actor set.

Example:

```r
Y <- list(
  t1 = Y1,
  t2 = Y2,
  t3 = Y3
)
```

where `rownames(Y[[t]])` and `colnames(Y[[t]])` identify active actors at time `t`.

Internally, convert input to a dynamic data object with:

- `actor_time`: actor-time rows with actor ID, time, unit ID, and time-specific index;
- `dyads`: dyadic observation rows or efficient equivalent per time point;
- `lineage`: predecessor/successor links between actor-time units.

For the first version, lineage is inferred from the same actor ID appearing in adjacent time points. Entry means no predecessor. Exit means no successor.

Keep the giant `mdsbm`-style matrix as a compatibility/internal representation only, not as the primary public API.

## Implementation doctrine

Build correctness first, but design for speed.

Initial implementation:

- R reference backend for data construction, time-specific GLM fitting, local scoring, objective accounting, and tests.
- Use R's GLM machinery (`stats::glm`, `stats::glm.fit`, `model.matrix`, family objects) rather than reimplementing GLM fitting in C++ immediately.

Performance path:

- Move local candidate scoring and membership updates to C++ after the R reference backend is correct and tested.
- Keep GLM coefficient fitting in R initially.
- Full C++ GLM fitting is a later optimization only if profiling proves it necessary.

Do not optimize before tests define the intended statistical behavior.

## Package engineering rules

- Use roxygen2 for documentation.
- Use testthat for tests.
- Use explicit namespace imports; do not put top-level `library()` calls inside package R files.
- Avoid `blockmodeling:::` internal functions. If unavoidable, isolate the use in one internal helper and document the risk.
- Keep public functions stable and small.
- Preserve existing `mdsbm` functions during package consolidation unless explicitly instructed otherwise.
- Do not claim CRAN readiness, production readiness, or statistical finality in early PRs.

## Testing requirements

Every implementation PR must include or update tests unless it is documentation-only.

Minimum first-test areas:

- dynamic input validation;
- actor entry/exit lineage construction;
- one-based user-facing partition labels;
- binomial family metadata sets `pseudo = FALSE`;
- PPML metadata sets `pseudo = TRUE` while still using common result names;
- deviance/log-probability scale conversion;
- local candidate score behavior on tiny examples;
- existing `mdsbm` functions still run on tiny examples after consolidation.

Tests should be small, deterministic, and fast. Long-running examples belong in vignettes or optional tests, not CI.

## Workflow rules

- Start each task from current `main` unless instructed otherwise.
- Create a feature branch.
- Commit only related files.
- Open a pull request if remote authentication allows it.
- Do not merge your own pull request.
- If GitHub authentication fails, first try to refresh authentication. If that fails, stop and ask the human to authenticate.
- You are running in a hardened WSL Ubuntu guest and may use passwordless sudo for local system dependencies inside the guest only.
- Do not access production secrets, host credentials, cloud credentials, or files outside the repository except for normal local dependency installation in the guest.

## Required final report format

Every execution task must end with:

- Branch name;
- Commit SHA(s);
- PR URL if opened;
- Files changed;
- Summary of implementation;
- Tests and checks run with exact commands and results;
- Local dependencies installed, if any;
- Skipped or blocked checks;
- Known risks;
- Recommended next task.
