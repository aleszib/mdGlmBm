# Task 002: Dynamic data layer with actor entry and exit

## Goal

Implement the first dynamic data representation: public list-of-matrices input converted to internal actor-time, dyad, and lineage structures.

Do this after Task 001 unless explicitly instructed otherwise.

## Required behavior

Implement a function such as:

```r
as_dynamic_network(Y, directed = TRUE, loops = FALSE, covariates = NULL)
```

where `Y` is a named or unnamed list of square matrices with row/column names.

The function should return an internal object with:

- `actor_time`;
- `dyads` or efficient equivalent;
- `lineage`;
- metadata for directedness, loops, times, and actor IDs.

## Required tests

- actor present at all times creates adjacent lineage links;
- actor entering has no predecessor;
- actor exiting has no successor;
- different actor sets across times are accepted;
- row/column name mismatch is rejected;
- non-square matrix is rejected;
- diagonal handling follows `loops`;
- missing values are preserved as missing/observed status, not converted silently to zero.

## Non-goals

- Do not implement split/merge lineage yet.
- Do not implement full optimizer.
- Do not implement long dyad-table input yet.
