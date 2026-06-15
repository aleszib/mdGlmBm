# Task 001: Package foundation and source consolidation

## Goal

Create the initial R package foundation for `dynGLMbm` by using the existing `mdsbm` package as the package skeleton and importing selected static GLM blockmodeling functions in a clean, tested way.

This task is package consolidation only. It should not implement the full dynamic GLM-Markov optimizer yet.

## Inputs expected in repository

The repository may contain or receive source ZIPs or extracted folders such as:

- `mdsbm-main.zip` or `mdsbm-main/`;
- `GLM_blockmodeling-main.zip` or `GLM_blockmodeling-main/`.

If archives are present, extract them inside a temporary/source-material folder and document the paths used.

## Required implementation

1. Use `mdsbm-main/mdsbm` as the starting package skeleton if available.
2. Preserve existing `mdsbm` exported functions unless they cannot build; if they cannot build, report exactly why.
3. Import selected static GLM functions from `GLM_blockmodeling-main` into package R files:
   - `ppml()`;
   - `optParGlm()`;
   - `optRandomParGlm()`;
   - `optRandomParRangeGlm()`;
   - `upAndDownSearch()`;
   - `balassaNorm()`.
4. Do not import duplicate copy files such as `GLMblockmodelingFunctions - Copy.R`.
5. Remove top-level `library()` calls from package R files.
6. Use DESCRIPTION `Imports` and explicit namespace qualification.
7. Avoid `blockmodeling:::` internal functions. If unavoidable, isolate in a single internal helper with a comment and test coverage.
8. Add roxygen2 documentation stubs sufficient for package checks.
9. Add testthat infrastructure.
10. Add basic result-metadata helpers for family/objective modes:
    - binomial sets `pseudo = FALSE`;
    - PPML sets `pseudo = TRUE`;
    - common result names are preserved.

## Required tests

Add fast tests for:

- package load;
- `ppml()` returns expected family-like behavior, or document the implemented PPML interface;
- `balassaNorm()` on a tiny matrix;
- binomial metadata sets `pseudo = FALSE`;
- PPML metadata sets `pseudo = TRUE` and uses common field names;
- one-based/zero-based label conversion helpers if added;
- one tiny static GLM call if feasible;
- one tiny existing `mdsbm` function call if feasible.

Do not force fragile tests that depend on long random searches. Use tiny deterministic examples.

## CI files

Ensure `.github/workflows/R-CMD-check.yaml` exists. Do not remove CI installed in Task 000.

## Required checks

Run, if package structure exists:

```sh
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("rcmdcheck", quietly=TRUE)) install.packages("rcmdcheck", repos="https://cloud.r-project.org")'
Rscript -e 'if (!requireNamespace("roxygen2", quietly=TRUE)) install.packages("roxygen2", repos="https://cloud.r-project.org")'
Rscript -e 'roxygen2::roxygenise()'
Rscript -e 'rcmdcheck::rcmdcheck(args="--no-manual", error_on="never")'
```

If system dependencies are missing, you may use passwordless sudo inside the hardened WSL Ubuntu guest. Document all packages installed.

## Non-goals

- Do not implement `fit_markov_glm_blockmodel()` yet unless only as a documented stub that errors clearly.
- Do not implement the dynamic GLM optimizer yet.
- Do not move local scoring to C++ yet.
- Do not claim CRAN readiness.
- Do not make PPML criteria look like true likelihood criteria without `pseudo = TRUE`.

## Branch and PR

Create a branch such as:

```text
feature/package-foundation
```

Commit only related files and open a PR if GitHub authentication allows it.

## Final report

Include:

- source material detected;
- files imported;
- functions exported;
- tests added;
- commands run and exact results;
- CI files present;
- known issues;
- recommended next task.
