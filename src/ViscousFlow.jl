"""
The `ViscousFlow` module is here

"""
module ViscousFlow

using DocStringExtensions
using Reexport

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

using LinearAlgebra
using SparseArrays

export NavierStokes, PointForce, SpatialGauss, Gaussian!, Gaussian,
       set_stepsizes, timestep,
       vorticity, velocity, streamfunction, nl, force, pressurejump


include("navier_stokes.jl")
include("plot_recipes.jl")



end
