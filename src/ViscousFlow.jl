"""
The `ViscousFlow` module is here

"""
module ViscousFlow

using DocStringExtensions
using Reexport

@reexport using CartesianGrids
@reexport using ImmersedLayers
@reexport using RigidBodyTools
@reexport using ConstrainedSystems
@reexport using GridUtilities


using LinearAlgebra
using SparseArrays

export NavierStokes, PointForce, SpatialGauss, Gaussian!, Gaussian,
       setstepsizes, timestep, timerange,
       vorticity, velocity, streamfunction, nl, force, pressurejump

abstract type PointMotionType end
abstract type StaticPoints <: PointMotionType end
abstract type MovingPoints <: PointMotionType end

abstract type FreestreamType end
abstract type StaticFreestream <: FreestreamType end
abstract type VariableFreestream <: FreestreamType end

include("navier_stokes.jl")
include("plot_recipes.jl")



end
