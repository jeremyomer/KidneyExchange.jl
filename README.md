# KidneyExchange


| **Documentation**                                                 | **Build Status**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-dev-img]][docs-dev-url] | [![][ci-img]][ci-url] [![][codecov-img]][codecov-url] |


This a Julia package to solve the deterministic kidney exchange problem. It provides five different solution methods, two of which are based on a branch-and-price algorithm. The other three methods consist in solving compact integer programming formulations. These methods are described in the following article:

Ayse N Arslan, Jérémy Omer, Fulin Yan. KidneyExchange.jl: A Julia package for solving the kidney exchange problem with branch-and-price. 2022. [⟨hal-03830810⟩](https://hal.inria.fr/hal-03830810)

The package also comes with three instance generators that allow to reproduce a benchmark similar to that used in the aforementioned article.  

## General informations

* Language : `Julia v1.6-1.8`
* Dependencies : `DelimitedFiles`, `GLPK`, `Graphs`, `HiGHs`, `JuMP`, `Random`, `Requires`, `TimerOutputs`, `DocStringExtensions`

## Input data

The package provides a parser (function `read_kep_file`) for the instances of the [PrefLib](https://www.preflib.org/) library publicly shared by John P. Dickerson and described in the article 

*Optimizing Kidney Exchange with Transplant Chains: Theory and Reality.* John P. Dickerson, Ariel D. Procaccia, Tuomas Sandholm; Proceedings of AAMAS; 2012. 

These instances must be downloaded through this [link](https://plmbox.math.cnrs.fr/d/4c2ac128352b489aafb5/) and stored in the data/preflib folder before solving them with various algorithms of the package.

Otherwise, three generators are provided with the package (functions `generate_saidman_instance`, `generate_heterogeneous_instance` and `generate_sparse_unos_instance`). One of them (`generate_saidman_instance`) will generate instances similar to those of the PrefLib library.

Users who wish to run our code on other existing instances may input them as a .wmd file and read them with the function `read_kep_file`; a description of the .wmd format is given in the [documentation of our parser](https://jeremyomer.github.io/KidneyExchange.jl/dev/functions/#KidneyExchange.read_kep_file-Tuple{AbstractString,%20AbstractString}).

## Integer programming solvers

The package was fully tested with the two commercial MIP solvers CPLEX and Gurobi. Those can both be downloaded and used under an Academic licence at https://www.ibm.com/academic/home and https://www.gurobi.com, respectively.

If you prefer running a fully open version of the package, it is possible to do so by using a chosen mixture of HiGHS, Clp, Cbc, and GLPK. The corresponding packages are documented at https://github.com/jump-dev/HiGHS.jl, https://github.com/jump-dev/Cbc.jl, https://github.com/jump-dev/Clp.jl and https://github.com/jump-dev/Glpk.jl. Be aware that the execution of the algorithms may take longer than communicated in our article if you do so, especially for large instances. This is true even for the branch-and-price algorithms which rely on the capacity of the solver to solve the relaxed master problem with integer variables. 

To choose the solver, you need to set the field `optimizer` of the `BP_params` or `MIP_params` structure with one of the following options:
- `HiGHS`: use exclusively HiGHS, i.e., both for integer programs (IPs) and linear programs (LPs)
- `GLPK`: use exclusively GLPK
- `Clp`: use Cbc for every integer program and solve the linear relaxations with Clp
- `Cbc`: exactly the same as `Clp` (***default for MIP approaches***)
- `GLPK-Cbc`: use Cbc for every integer program and solve the linear relaxations with GLPK (***default for branch-and-price***)
- `CPLEX`: use exclusively CPLEX (_requires a licensed installation of CPLEX_)
- `Gurobi`: use exclusively Gurobi (_requires a licensed installation of Gurobi_)

Some parameters of the solvers (e.g. the number of threads used by the solver) can be set when calling the constructors of the parameters used by the desired algorithm. Additional details are given in the documentation of [MIP_params](https://jeremyomer.github.io/KidneyExchange.jl/dev/types/#KidneyExchange.MIP_params) for MIP approaches and [BP_params](https://jeremyomer.github.io/KidneyExchange.jl/dev/types/#KidneyExchange.BP_params) for branch-and-price algorithms.

Packages Clp, Cbc, Gurobi and CPLEX are not included among the dependencies so they must be added and loaded (with using) if you wish to use one of them. This was necessary due to some compatibility issues with Clp and Cbc, and because CPLEX and Gurobi require a licensed installation. 

## Basic usage (with Julia REPL)

Below, we give some examples of the basic usage of the package. All the functions that are called below and the other callable functions are also described in the [documentation](https://jeremyomer.github.io/KidneyExchange.jl/dev/) of the package. This documentation includes the description of the input and output parameters of each function.

To use the package, first load the module
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
Using Revise will allow the user to modify the code and see the effect of these modifications directly.

If you wish to use a commercial integer programming solver such as Gurobi or CPLEX, you need to have an active license and import the corresponding package. For instance: 
 ```
julia> using Gurobi
```

Generate a sparse instance with 500 pairs of incompatible donors and receivers and 25 non-directed donors. The corresponding input files will be created in the "data/sparse/" folder of the package. 

`generate_sparse_unos_instance(500, 25, 1)`

Solve the instance with branch-and-price. The basic usage of this function requires as input, the instance (below "sparse/sparse_500_25_1"), and the parameters K and L (below 3 and 4, respectively). 

`solve_with_BP("sparse/sparse_500_25_1", 3, 4);`

As always with Julia, the first time the solve function is called, it takes extra time. 

In a more advanced use case additional parameters can be changed through the object BP_params. For instance, specify another branch-and-price formulation (the second true is to keep verbosity)

`solve_with_BP("sparse/sparse_500_25_1", 3, 4, BP_params(true, true));`

Solve the instance using one of its compact MIP formulations (HPIEF by default). The basic usage of this function requires as input, the instance (below "sparse/sparse_500_25_1"), and the parameters K and L (below 3 and 4, respectively). 

`solve_with_mip("sparse/sparse_500_25_1", 3, 4);`

In a more advanced use case additional parameters can be changed through the object MIP_params. For instance, specify another MIP formulation (true is to keep verbosity)

`solve_with_mip("sparse/sparse_500_25_1", 3, 4, MIP_params(KidneyExchange.EXTENDED_EDGE, true));`

The parameters of branch-and-price algorithms can be listed with
 ```
julia> ?
help?> BP_params
```
And the same can be done with MIP_params for the solution of compact formulations. For instance, if Gurobi has been loaded and you wish to run the branch-and-price algorithm with Gurobi, this can be done with:
 ```
julia> bp_params = BP_params()
julia> bp_params.optimizer = "Gurobi"
julia> solve_with_BP("sparse/sparse_500_25_1", 3, 4, bp_params);
```
after having imported the Gurobi package in your environment. 

## Advanced usage and running tests 

The package can also be used without adding it to the Julia environment and resorting to the Julia REPL. To do so, download the code folder and make sure that all dependencies of the package are present in your Julia environment. 

A set of example script files are provided in the code folder under the name KEPTestxxx.jl for each method tested in our manuscript. 

These files require an instance, and parameters K and L as input. They then solve the provided instance with the chosen algorithm, using default parameter settings (with Gurobi as the MIP solver), and write the output to a dedicated .csv file. 
As an example, the KEPTestBP.jl can be run with:

`julia KEPTestBP.jl preflib/MD-00001-00000011 3 6`

provided that the path to the Julia executable is correctly configured in your environment (see [here](https://julialang.org/downloads/platform/)) and the PrefLib instances are stored in the corresponding data folder. 

We additionally remark that our KEPTestxxx.jl files use the Gurobi solver as the MIP solver, as such their execution requires the package Gurobi.jl to be added to the Julia environment with a licensed installation of the Gurobi solver.    

Interested users can then embed this execution into a shell script in order to test a number of instances and tabulate the results. For instance, on a Unix system 

```
for K in {3..4}
do
    julia KEPTestBP.jl preflib/MD-00001-00000011 $K 6
done
```
will solve the provided PrefLib instance for K=3,4 and create two .csv files as output. 

One can additionally set all the desired attributes of BP_params and MIP_params directly inside KEPTestxxx.jl files. For instance, by adding the lines  
```
bp_params = BP_params()
bp_params.optimizer = "Gurobi"
 ```
or 
```
mip_params = MIP_params()
mip_params.optimizer = "Gurobi"
```
the default solver is changed to Gurobi. Users who do not wish to use the Gurobi solver can start with modifying this line to replace Gurobi with a solver of their choice. 

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://jeremyomer.github.io/KidneyExchange.jl/dev/

[ci-img]: https://github.com/jeremyomer/KidneyExchange.jl/workflows/CI/badge.svg?branch=master
[ci-url]: https://github.com/jeremyomer/KidneyExchange.jl/actions?query=workflow%3A%22CI%22

[codecov-img]: https://codecov.io/gh/jeremyomer/KidneyExchange.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/jeremyomer/KidneyExchange.jl
