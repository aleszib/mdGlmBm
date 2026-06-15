# Data Representation Strategy

## Public input for the first implementation

The primary public input should be a named list of square matrices, one per time point.

```r
Y <- list(
  t1 = Y1,
  t2 = Y2,
  t3 = Y3
)
```

Each matrix must have row and column names identifying active actors at that time point.

Actors may enter or leave over time. The same actor ID in adjacent time points creates a lineage link. Missing actor ID at the previous time point means entry. Missing actor ID at the next time point means exit.

## Accepted network storage

Initial implementation should support base R dense matrices. Design validation so that sparse Matrix classes can be added without changing the public API.

Do not require the user to construct the giant `mdsbm` matrix.

## Internal representation

Convert public input to an internal object with these components:

### actor_time

A table with one row per actor-time unit:

- `unit_id`;
- `actor_id`;
- `time_id`;
- `time_label`;
- `row_index`;
- `active`.

### dyads

A dyadic representation per time point or equivalent efficient structure:

- `time_id`;
- `sender_unit`;
- `receiver_unit`;
- `value`;
- `observed`;
- `weight`;
- optional dyadic covariates.

### lineage

A table of dynamic links between actor-time units:

- `from_unit`;
- `to_unit`;
- `from_time`;
- `to_time`;
- `relation`;
- `weight`.

For the first implementation, lineage uses only one-to-one same-ID links between adjacent time points.

## Long-term split/merge compatibility

Do not implement split/merge first. However, the lineage table should be general enough that future split/merge relations can be represented as multiple predecessor or successor links.

## Diagonal and directedness

The first implementation should default to directed networks and ignore diagonals unless explicitly configured otherwise. Diagonal handling must be tested and documented.

## Missing values

Missing dyads should be treated differently from observed zero ties. The input validator should preserve this distinction.
