## ViscousFlow.jl

_A framework for simulating viscous incompressible flows_

| Documentation | Build Status |
|:---:|:---:|
| [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://jdeldre.github.io/ViscousFlow.jl/latest) | [![Build Status](https://img.shields.io/travis/jdeldre/ViscousFlow.jl/master.svg?label=linux)](https://travis-ci.org/jdeldre/ViscousFlow.jl) [![Build status](https://img.shields.io/appveyor/ci/jdeldre/whirl-jl/master.svg?label=windows)](https://ci.appveyor.com/project/jdeldre/whirl-jl/branch/master) [![codecov](https://codecov.io/gh/jdeldre/ViscousFlow.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jdeldre/ViscousFlow.jl) |

## About the package

The purpose of this package is to enable easy setup and solution of viscous incompressible flows. Documentation can be found at https://jdeldre.github.io/ViscousFlow.jl/latest.

**ViscousFlow.jl** is registered in the general Julia registry. To install in julia `0.6`, type
```julia
julia> Pkg.add("ViscousFlow")
```
in the Julia REPL.

In julia `0.7` or `1.0`, enter the package manager by typing `]` and then type,
e.g.,
```julia
(v1.0) pkg> add ViscousFlow
```

Then, in any version, type
```julia
julia> using ViscousFlow
```
For examples, consult the documentation or see the example Jupyter notebooks in the Examples folder.

![](https://github.com/jdeldre/ViscousFlow.jl/raw/master/cylinderRe400.gif)
