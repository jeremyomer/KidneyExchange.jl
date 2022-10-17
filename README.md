# KidneyExchange


| **Documentation**                                                 | **Build Status**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] |


This a Julia package to solve the deterministic kidney exchange problem. It provides five different solving methods, two of which are based on a branch-and-price algorithms. The other three methods consist in solving integer programming compact formulations. These methods are described in the following article:
***give the citation of the HAL report once submitted***

The package also comes with three generators that allow to reproduce a benchmark similar to that used in the aforementionned report.  

## General informations

* Language : `Julia v1.6`
* Dependencies : `Cbc`, `Clp`, `CPLEX`, `DelimitedFiles`, `Distributions`, `GLPK`, `Graphs`, `Gurobi`, `JuMP`, `Random`, `TimerOutputs`

## Input data

Before running the code, the kidney exchange data publicly shared by John Dickerson must be downloaded at https://www.preflib.org/data/MD/00001 and stored in the data/preflib folder. The data shared by John Dickerson is described in:
*Optimizing Kidney Exchange with Transplant Chains: Theory and Reality.* John P. Dickerson, Ariel D. Procaccia, Tuomas Sandholm; Proceedings of AAMAS; 2012

Otherwise, three generators are provided with the package. One of them will generate instances similar to the PrefLib.

## Integer programming solvers

The package was fully tested with the two commercial MIP solvers CPLEX and Gurobi. Those can both be downloaded and used under Academic licence at https://www.ibm.com/academic/home and https://www.gurobi.com.

If you prefer running a fully open version of the package, it is possible by using a chosen mixture of Clp, Cbc, and GLPK. Be aware that the execution of the algorithms may take longer than communicated in our article if you do so. In particular, for large instances, this is even true for the branch-and-price algorithms which rely on the capacity of the solver to solve the relaxed master problem with integer variables. The corresponding packages are documented at https://github.com/jump-dev/Cbc.jl, https://github.com/jump-dev/Clp.jl and https://github.com/jump-dev/Glpk.jl. 

To choose the solver, you need set the field `optimizer` of the `BP_params` or `MIP_params` structure with one of the following options:
- `GLPK`: use exclusively GLPK
- `Clp`: use Cbc for every integer program and solve the linear relaxations with Clp
- `Cbc`: exactly the same as `Clp` (***default for MIP approaches***)
- `GLPK-Cbc`: use Cbc for every integer program and solve the linear relaxations with GLPK (***default for branch-and-price***)
- `CPLEX`: use exclusively CPLEX (_requires a licensed installation of CPLEX_)
- `Gurobi`: use exclusively Gurobi (_requires a licensed installation of Gurobi_)


## Basic usages

Load the module
```
julia> using Pkg
pkg> add "https://github.com/jeremyomer/KidneyExchange.jl.git"
julia> using KidneyExchange
```

The simpler option is to download the code folder and navigate to the folder before opening Julia. The package may then be loaded as:
 ```
julia> using Revise
pkg> activate .
julia> using KidneyExchange
```
Using Revise will allow to modify the code and see directly the effect of these modifications.

If you wish to use a commercial integer programming solver such as Gurobi or CPLEX, you need to have an active license and import the corresponding package. For instance: 
 ```
julia> using Gurobi
```

Generate an instance with 500 pairs of incompatible donors and receivers and 25 altruist donors
`generate_sparse_unos_instance(500, 25, 1)`

Solve the instance with branch-and-price
`solve_with_BP("sparse/sparse_500_25_1", 3, 4);`

Specify another branch-and-price formulation (the second true is to keep verbosity)
`solve_with_BP("sparse/sparse_500_25_1", 3, 4, BP_params(true, true));`

Solve the instance using one of its compact MIP formulations (HPIEF by default)
`solve_with_mip("sparse/sparse_500_25_1", 3, 4);`

Specify another MIP formulation (true is to keep verbosity)
`solve_with_mip("sparse/sparse_500_25_1", 3, 4, MIP_params(KidneyExchange.EXTENDED_EDGE, true));`

The parameters of branch-and-price approaches can be listed with
 ```
julia> ?
help?> BP_params
```
And the same can be done with MIP_params for the solution of compact formulations. For instance, if Gurobi has been loaded and you wish to solve all linear and integer programs with Gurobi, this can be done with:
 ```
julia> bp_params = BP_params()
julia> bp_params.optimizer = "Gurobi"
julia> solve_with_BP("sparse/sparse_500_25_1", 3, 4, bp_params);
```

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jeremyomer.github.io/KidneyExchange.jl/dev/

[ci-img]: https://github.com/jeremyomer/KidneyExchange.jl/workflows/CI/badge.svg?branch=master
[ci-url]: https://github.com/jeremyomer/KidneyExchange.jl/actions?query=workflow%3A%22CI%22

[codecov-img]: https://codecov.io/gh/jeremyomer/KidneyExchange.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jeremyomer/KidneyExchange.jl
