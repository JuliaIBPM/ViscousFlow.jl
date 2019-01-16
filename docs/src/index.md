# Whirl

*a framework for simulating viscous incompressible flows*

The objective of this package is to allow easy setup and fast simulation of incompressible
flows, particularly those past bodies in motion. The package provides
tools for
- constructing grids and body shapes,
- using the operators on those grids,
- specifying the relevant parameters and setting their values,
- solving the problem, and finally,
- visualizing and analyzing the results.

The underlying grids are uniform and Cartesian, allowing the use of the lattice
Green's function (LGF) for inverting the Poisson equation; the diffusion operators are
solved with the integrating factor (Liska and Colonius ref). Many of the core aspects
of the fluid-body interaction are based on the immersed boundary projection method,
developed by Taira and Colonius (ref). The coupled fluid-body interactions are based
on the work of Wang and Eldredge (ref).


## Installation

This package requires Julia `0.6` and above.
To install, simply run
```julia
julia> Pkg.clone("https://github.com/jdeldre/Whirl.jl.git","Whirl")
```
in the Julia REPL.
Since this package is still under heavy development, you should run
```julia
julia> Pkg.test("Whirl") # might take some time
```
to make sure things are working as intended and
```julia
julia> Pkg.update()
```
to get the most recent version of the library and its dependencies.

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that too to follow the examples.

## Basic Usage

Do something here.
