# Project Memory: dynGLMbm

## Current purpose

Create a professional-grade R package for dynamic GLM-based blockmodeling. The project combines ideas and code from:

- `mdsbm`: dynamic Markov-style stochastic blockmodeling implemented efficiently in C++/Rcpp;
- `GLM_blockmodeling`: static GLM-based blockmodeling using deviance-style local optimization.

The human lead wants speed, tests, professional package structure, and a long-term path toward an efficient dynamic GLM-Markov estimator.

## Settled statistical decisions

### Dynamic structure

The long-term estimator should preserve the `mdsbm` dynamic logic: actor-time memberships are updated using within-time observation fit and previous/next dynamic membership terms.

### Observation model

The static/within-time observation model should be GLM-based blockmodeling instead of a simple Bernoulli block-density model.

### Time handling

For the first dynamic GLM-Markov estimator, the GLM observation models are fitted independently at each time point. The different time points are coupled through latent memberships and Markov transition terms, not through shared GLM coefficients.

### Local optimization scale

The GLM code uses individual deviance contributions for local reassignment. This should be preserved.

The dynamic terms should be converted to the same scale:

```text
candidate loss = individual GLM deviance + (-2 * dynamic log-probability terms)
```

Smaller is better.

### Final reporting scale

After a full iteration, compute final objective components for reporting and diagnostics:

```text
total logLik = sum_t GLM logLik_t + Markov transition logLik + membership/prior logLik
```

The local optimizer does not need to use the global logLik directly at every candidate move.

### Families

First dynamic families:

- `binomial`;
- `ppml`.

Later:

- `poisson`;
- Gaussian/normal;
- additional true likelihood or pseudo-likelihood families if statistically documented.

For PPML, use the same result field names as for other families, but set `pseudo = TRUE` and print a note that likelihood criteria are pseudo-likelihood-based.

### Actors over time

Actors may enter and leave networks. Do not assume a fixed actor set. Use actor IDs and an internal actor-time/lineage representation.

Split/merge of units over time is a long-term possibility, but not a first implementation requirement.

## Settled implementation decisions

### Public data interface

Use a list of time-specific networks with actor IDs in row/column names.

### Internal data object

Convert input to an internal dynamic object containing at least:

- actor-time table;
- dyad representation;
- lineage table.

### Main optimizer strategy

Use Option 1: `mdsbm`-style local actor-time reassignment using individual GLM deviance plus Markov penalties.

Do not use static GLM followed only by post-hoc smoothing as the final model. Do not make exact global-refit greedy optimization the production default, because speed is important.

### R vs C++

Initial implementation should use R's GLM machinery and a clear R reference backend. The performance path is a C++ backend for local scoring and membership updates, with GLM coefficient fitting remaining in R initially.

## Development milestones

### Milestone 0: Instruction and CI bootstrap

Install `AGENTS.md`, project memory, statistical doctrine, task files, PR template, and GitHub Actions CI.

### Milestone 1: Package foundation

Consolidate source material into an R package skeleton. Preserve existing `mdsbm` functions. Import and clean selected static GLM functions. Add basic tests and CI.

### Milestone 2: Dynamic data layer

Implement input validation and conversion from list-of-matrices to actor-time/dyad/lineage representation. Support actor entry and exit.

### Milestone 3: Family/objective metadata

Implement family handling for binomial and PPML. Common result names. `pseudo` flag.

### Milestone 4: R reference dynamic scorer

Implement local deviance-scale candidate scoring for actor-time memberships with Markov penalties.

### Milestone 5: R reference optimizer

Implement a first working dynamic GLM-Markov optimizer using independent time-specific GLMs and mdsbm-style dynamic reassignment.

### Milestone 6: C++ scoring backend

Move local candidate scoring and membership updates to C++ after the R reference behavior is tested.
