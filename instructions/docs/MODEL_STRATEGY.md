# Statistical Model Strategy

## Target model class

The target model is a dynamic latent blockmodel with GLM-based within-time observations and Markov-style temporal coupling of memberships.

For time points `t = 1, ..., T`:

- `Y_t` is the observed network at time `t`;
- `Z_t` is the vector of latent block memberships for actors active at time `t`;
- `beta_t` are time-specific GLM coefficients;
- `P` is the Markov transition matrix for latent block evolution.

The first target objective is:

```text
sum_t logLik_GLM(Y_t | Z_t, beta_t)
+ sum_lineage log P(Z_successor | Z_predecessor)
+ membership/prior terms if used
```

Each `beta_t` is fitted independently conditional on `Z_t`.

## Local update criterion

The practical optimizer uses deviance-scale local scoring because the static GLM code optimizes partitions through individual deviance contributions.

For actor-time unit `(i, t)` and candidate block `k`:

```text
score(i, t, k) =
  individual_deviance_t(i, k | beta_t)
  - 2 * log previous_transition_probability
  - 2 * log next_transition_probability
  - 2 * log membership_prior_probability
```

Smaller score is better.

Boundary cases:

- entry: no predecessor term;
- exit: no successor term;
- isolated actor-time unit: only GLM and prior terms;
- future split/merge: dynamic terms become sums over all lineage predecessors/successors.

## Deviance and log-likelihood compatibility

For true likelihood GLMs, deviance is on a `-2 logLik + constant` scale. Therefore, dynamic log-probability terms must be multiplied by `-2` before being added to individual deviance scores.

Never add raw deviance and raw log-likelihood together.

## Final objective reporting

After a full sweep/iteration, compute and store:

- `logLik`;
- `deviance`;
- `BIC`;
- `ICL`;
- `objective`;
- `pseudo`.

For PPML, use the same field names and set `pseudo = TRUE`. Print a note that likelihood-based criteria are pseudo-likelihood-based.

## First supported families

1. `binomial` with `pseudo = FALSE`.
2. `ppml` with `pseudo = TRUE`.

Planned later:

- `poisson`;
- Gaussian/normal;
- other families only after the objective and reporting interpretation are documented.

## Non-goals for the first dynamic optimizer

- No pooled/shared GLM across all time points.
- No exact refit for every candidate move as the default.
- No split/merge lineage estimation.
- No arbitrary formula support in the first C++ scoring backend.
- No CRAN/readiness claims.
