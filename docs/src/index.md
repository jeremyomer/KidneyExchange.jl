# KidneyExchange.jl

This a Julia package to solve the deterministic kidney exchange
problem. It provides five different solution methods, two of which
are based on a branch-and-price algorithm. The other three methods
consist in solving compact integer programming formulations. These
methods are described in [omer:hal-03830810](@Citet).

The package also comes with three instance generators that allow
to reproduce a benchmark similar to that used in the aforementionned
report.

## Input data

The package provides a parser (function [`read_kep_file`](@ref)) for the
instances of the [PrefLib](https://www.preflib.org/dataset/00036) library
publicly shared by John P. Dickerson and described in [dickersonOptimizingKidneyExchange2012](@Citet).
Those instances must be downloaded from the PrefLib
[website](https://www.preflib.org/dataset/00036) and stored in the
`data/preflib` folder before solving them with the algorithms of the
package.

Otherwise, three generators functions are provided with the package :

- [`generate_saidman_instance`](@ref)
- [`generate_heterogeneous_instance`](@ref)
- [`generate_sparse_unos_instance`](@ref)). 

!!! note 
   [`generate_saidman_instance`](@ref)) will generate instances similar to those of the [PrefLib library](https://www.preflib.org/dataset/00036).

Users who wish to run our code on other existing instances may input
them as `.wmd` file and read them with [`read_wmd_file`](@ref).

## Integer programming solvers

The package was fully tested with the two commercial MIP solvers
CPLEX and Gurobi. Those can both be downloaded and used under an
Academic licence at <https://www.ibm.com/academic/home> and
<https://www.gurobi.com>.

If you prefer running a fully open version of the package, it is
possible to do so by using a chosen mixture of HiGHS, Clp, Cbc, and
GLPK. Be aware that the execution of the algorithms may take longer
than communicated in our article if you do so. In particular, for
large instances, this is even true for the branch-and-price algorithms
which rely on the capacity of the solver to solve the relaxed master
problem with integer variables. The corresponding packages are
documented at [HiGHS.jl](https://github.com/jump-dev/HiGHS.jl), 
[Cbc.jl](https://github.com/jump-dev/Cbc.jl), [Clp.jl](https://github.com/jump-dev/Clp.jl)
and [Glpk.jl](https://github.com/jump-dev/Glpk.jl).

To choose the solver, you need to set the field `optimizer` of the [`BP_params`](@ref) or [`MIP_params`](@ref) structure with one of the following options:
- `HiGHS`: use exclusively HiGHS, i.e., both for integer programs (IPs) and linear programs (LPs)
- `GLPK`: use exclusively GLPK
- `Clp`: use Cbc for every integer program and solve the linear relaxations with Clp
- `Cbc`: exactly the same as `Clp` (***default for MIP approaches***)
- `GLPK-Cbc`: use Cbc for every integer program and solve the linear relaxations with GLPK (***default for branch-and-price***)
- `CPLEX`: use exclusively CPLEX (_requires a licensed installation of CPLEX_)
- `Gurobi`: use exclusively Gurobi (_requires a licensed installation of Gurobi_)

Some parameters of the solvers (e.g. the number of threads used by
the solver) can be set when calling the constructors of the parameters
used by the desired algorithm. Additional details are given in the
documentation of [`MIP_params`](@ref) and [`BP_params`](@ref).

Packages Clp, Cbc, Gurobi and CPLEX are not included among the
dependencies so they must be added and loaded (with using) if you
wish to use one of them. This was necessary due to some compatibility
issues with Clp and Cbc, and because CPLEX and Gurobi require a
licensed installation.

## Quickstart

```@example 1
using KidneyExchange
```

Generate an instance with 500 pairs of incompatible donors and
receivers and 25 non-directed donors. The corresponding input files
will be created in the "data/sparse/" folder of the package.

```@example 1
generate_sparse_unos_instance(500, 25, 1)
```
Solve the instance with branch-and-price. The basic usage of this
function requires as input, the instance (below "sparse/sparse_500_25_1"),
and the parameters K and L (below 3 and 4, respectively).

```@example 1
solve_with_BP("sparse/sparse_500_25_1", 3, 4)
```

As always with Julia, the first time the solve function is called, it takes extra time. 

In a more advanced use case additional parameters can be changed
through the object [`BP_params`](@ref). For instance, specify another
branch-and-price formulation (the second true is to keep verbosity)

```@example 1
solve_with_BP("sparse/sparse_500_25_1", 3, 4, BP_params(true, true))
```

Solve the instance using one of its compact MIP formulations (HPIEF
by default). The basic usage of this function requires as input,
the instance (below "sparse/sparse_500_25_1"), and the parameters
K and L (below 3 and 4, respectively).

```@example 1
solve_with_mip("sparse/sparse_500_25_1", 3, 4)
```

In a more advanced use case additional parameters can be changed
through the object [`MIP_params`](@ref). For instance, specify another MIP
formulation (true is to keep verbosity)

```@example 1
solve_with_mip("sparse/sparse_500_25_1", 3, 4, MIP_params(KidneyExchange.EXTENDED_EDGE, true))
```

## Index

```@index
```
