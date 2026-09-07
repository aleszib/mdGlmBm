# Statistical Model Strategy

## Target model class

The target model is a dynamic latent blockmodel with GLM-based within-time observations and Markov-style temporal coupling of memberships.

For time points `t = 1, ..., T`:

- `Y_t` is the observed network at time `t`;
- `Z_t` is the vector of latent block memberships for actors active at time `t`;
- `beta_t` are time-specific GLM coefficients;
- `P` is the augmented Markov transition model for latent block evolution.
  It has `K` substantive states plus emerging (`E`) and vanishing (`V`)
  transition-only states.

The first observation model is explicitly time-specific: `Y_t | Z_t, beta_t`. Each `beta_t` is fitted independently at time `t`; temporal dependence is in `Z_t`, not shared GLM coefficients.

The first target objective is:

```text
sum_t logLik_GLM(Y_t | Z_t, beta_t)
+ sum_transition_events log P(state_successor | state_predecessor)
+ initial-time substantive-state terms if used
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

- entry at an observed later time: use the `E -> k` transition term;
- observed exit before a later time: use the `k -> V` transition term;
- an actor present at only one intermediate time uses both terms;
- actors present at the first time use the initial substantive-state distribution;
- future split/merge: dynamic terms become sums over all lineage predecessors/successors.

Transition probabilities are smoothed before taking logs so numerical zeros do not create unintended infinite candidate penalties. Auxiliary states are never GLM blocks or ordinary candidate labels. The complete objective uses an initial-state term only at the first observed time and does not apply a generic membership prior independently at every actor-time. The current R reference retains pooled substantive transition estimates for compatibility with the established optimizer, while retaining boundary-specific augmented event estimates for entry/exit auditing.

The `prior` options therefore mean: `"empirical"` estimates the initial
substantive distribution from time-one memberships, `"uniform"` uses equal
initial probabilities, and `"none"` omits the explicit initial-state term.
The latter is a conditional/debug specification rather than the preferred
complete model.

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

The later C++ scoring backend must reproduce the R reference scores and objective components on fixed deterministic inputs.
