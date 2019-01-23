# ViscousFlow.jl

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
solved with the integrating factor (Liska and Colonius [^1]). Many of the core aspects
of the fluid-body interaction are based on the immersed boundary projection method,
developed by Taira and Colonius [^2]. The coupled fluid-body interactions are based
on the work of Wang and Eldredge [^3].

![](https://github.com/jdeldre/ViscousFlow.jl/raw/master/cylinderRe400.gif)

## Installation

This package works on Julia `0.6`, `0.7` and `1.0` and is registered in the general Julia registry. To install in julia `0.6`, type
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

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that, too, to follow the examples.

## References

[^1]: Liska, S. and Colonius, T. (2017) "A fast immersed boundary method for external incompressible viscous flows using lattice Green's functions," *J. Comput. Phys.*, 331, 257--279.

[^2]: Taira, K. and Colonius, T. (2007) "The immersed boundary method: a projection approach," *J. Comput. Phys.*, 225, 2118--2137.

[^3]: Wang, C. and Eldredge, J. D. (2015) "Strongly coupled dynamics of fluids and rigid-body systems with the immersed boundary projection method," *J. Comput. Phys.*, 295, 87--113. [(DOI)](https://doi.org/10.1016/j.jcp.2015.04.005).
