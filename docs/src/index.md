# KidneyExchange.jl

Documentation for KidneyExchange.jl

## Quickstart

```@example 1
using KidneyExchange
```

- Generate an instance with 500 pairs of incompatible donors and receivers and 25 altruist donors

```@example 1
generate_sparse_unos_instance(500, 25, 1)
```
- Solve the instance with branch-and-price

```@example 1
solve_with_BP("sparse/sparse_500_25_1", 3, 4)
```
- Specify another branch-and-price formulation (the second true is to keep verbosity)

```@example 1
solve_with_BP("sparse/sparse_500_25_1", 3, 4, BP_params(true, true))
```

- Solve the instance using one of its compact MIP formulations (HPIEF by default)

```@example 1
solve_with_mip("sparse/sparse_500_25_1", 3, 4)
```

- Specify another MIP formulation (true is to keep verbosity)

```@example 1
solve_with_mip("sparse/sparse_500_25_1", 3, 4, MIP_params(KidneyExchange.EXTENDED_EDGE, true))
```

## Index

```@index
```
