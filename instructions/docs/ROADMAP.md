# Roadmap

## 0. Bootstrap instructions and CI

Install project instructions, statistical doctrine, PR template, and GitHub Actions workflow.

## 1. Package foundation

Use the existing `mdsbm` package as the starting package skeleton. Import selected static GLM blockmodeling functions. Add tests and CI.

## 2. Dynamic data layer

Implement the public list-of-matrices input and internal actor-time/dyad/lineage representation. Support actor entry and exit.

## 3. Family metadata and result schema

Implement family handling for binomial and PPML. Standardize result fields and `pseudo` behavior.

## 4. R reference local scorer

Implement deviance-scale actor-time candidate scoring with Markov penalties.

## 5. R reference dynamic optimizer

Implement the first `fit_markov_glm_blockmodel()` using independent time-specific GLMs and mdsbm-style dynamic reassignment.

## 6. C++ scoring backend

Optimize local scoring and membership sweeps in C++/Rcpp. Keep R backend for correctness tests.

## 7. Extended families

Add Poisson and Gaussian/normal after binomial and PPML are tested.

## 8. Long-term extensions

- sparse matrix support;
- long dyad table input;
- optional exact verification mode for small networks;
- split/merge lineage support;
- partially pooled GLM observation models;
- performance profiling and benchmark suite.
