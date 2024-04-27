# Creating a sysimage using PackageCompiler


## Create a new julia environment to create your image

```
$ mkdir KidneyExchangeImage
$ cd KidneyExchangeImage/
$ julia -q --startup-file=no
julia> using PackageCompiler
 │ Package PackageCompiler not found, but a package named PackageCompiler is available from a
 │ registry.
 │ Install package?
 │   (@v1.10) pkg> add PackageCompiler
 └ (y/n/o) [y]:
   Resolving package versions...
   Installed PackageCompiler ─ v2.1.17
    Updating `~/.julia/environments/v1.10/Project.toml`
  [9b87118b] + PackageCompiler v2.1.17
    Updating `~/.julia/environments/v1.10/Manifest.toml`
  [9b87118b] + PackageCompiler v2.1.17
  [4af54fe1] + LazyArtifacts
(@v1.10) pkg> activate .
  Activating new project at `~/JuliaProjects/KidneyExchangeImage`
```

## Install the package KidneyExchange.jl

```
(KidneyExchangeImage) pkg> add https://github.com/jeremyomer/KidneyExchange.jl.git
     Cloning git-repo `https://github.com/jeremyomer/KidneyExchange.jl.git`
    Updating git-repo `https://github.com/jeremyomer/KidneyExchange.jl.git`
   Resolving package versions...
    Updating `~/JuliaProjects/KidneyExchangeImage/Project.toml`
  [1df6b590] + KidneyExchange v0.1.0 `https://github.com/jeremyomer/KidneyExchange.jl.git#master`
    Updating `~/JuliaProjects/KidneyExchangeImage/Manifest.toml`
⌅ [ec485272] + ArnoldiMethod v0.2.0
  [6e4b80f9] + BenchmarkTools v1.5.0
  [523fee87] + CodecBzip2 v0.8.2
  [944b1d66] + CodecZlib v0.7.4
  ...
  ...
  [8e850ede] + nghttp2_jll v1.52.0+1
  [3f19e933] + p7zip_jll v17.4.0+2
        Info Packages marked with ⌅ have new versions available but compatibility constraints restrict them from upgrading. To see why use `status --outdated -m`
Precompiling project...
  1 dependency successfully precompiled in 7 seconds. 132 already precompiled.
```

## Create a an example to precompile the package in your image

```
shell> cat precompile.jl

using KidneyExchange
using Printf
using Downloads
using HiGHS

bp_params = BP_params()
bp_params.optimizer = "HiGHS"
timer = TimerOutput()
bp_params.verbose = false
max_time = 7200.0

dataset = @sprintf "%08d" 2  #  (1 - 310) instance index of synthetic kidney donor pools

filename = "00036-" * dataset
dat_file = filename * ".dat"
wmd_file = filename * ".wmd"

Downloads.download("https://www.preflib.org/static/data/kidney/" * dat_file, dat_file)
Downloads.download("https://www.preflib.org/static/data/kidney/" * wmd_file, wmd_file)

cycle_limit, chain_limit =  3, 2
bp_status, graph_info, subgraph_info = solve_with_BP(filename, cycle_limit, chain_limit, bp_params, timer, max_time)
```

## Install your solver

```
(KidneyExchangeImage) pkg> add HiGHS
   Resolving package versions...
    Updating `~/JuliaProjects/KidneyExchangeImage/Project.toml`
  [87dc4568] + HiGHS v1.9.0
  No Changes to `~/JuliaProjects/KidneyExchangeImage/Manifest.toml`
```

## Create your image

```
julia> create_sysimage(["KidneyExchange", "HiGHS"]; sysimage_path="KidneyExchangeSysimage.so", precompile_execution_file="precompile.jl")
Precompiling project...
  123 dependencies successfully precompiled in 82 seconds. 10 already precompiled.
✔ [02m:41s] PackageCompiler: compiling incremental system image

```


## Use it 

```
$ ls
00036-00000001.dat		00036-00000002.wmd		MD-00001-00000002.wmd		precompile.jl
00036-00000001.wmd		KidneyExchangeSysimage.so	Manifest.toml			solve_with_bp.jl
00036-00000002.dat		MD-00001-00000002.dat		Project.toml			solve_with_bppicef.jl
00036-00000002.wmd		Project.toml

$ julia --startup-file=no -J KidneyExchangeSysimage.so solve_with_bp.jl
  0.007080 seconds (11.89 k allocations: 533.320 KiB)
  0.005073 seconds (10.47 k allocations: 441.906 KiB)
```
