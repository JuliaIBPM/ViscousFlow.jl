"""
The `ViscousFlow` module is here

"""
module ViscousFlow

using DocStringExtensions
using Reexport

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems
@reexport using GridUtilities

using LinearAlgebra
using SparseArrays

export NavierStokes, PointForce, SpatialGauss, Gaussian!, Gaussian,
       setstepsizes, timestep, timerange,
       vorticity, velocity, streamfunction, nl, force, pressurejump


include("navier_stokes.jl")
include("plot_recipes.jl")



end
