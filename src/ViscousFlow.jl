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

export NavierStokes, PulseParams, PointForce, SpatialDGaussian,
       ExternalFlow, InternalFlow, ExternalInternalFlow,
       setstepsizes, timestep, timerange, newstate,
       update_immersion_operators!,
       vorticity, velocity, velocity!, streamfunction, streamfunction!,
       scalarpotential, scalarpotential!, convective_derivative, convective_derivative!,
       pressure, traction, force, pressurejump

abstract type PointMotionType end
abstract type StaticPoints <: PointMotionType end
abstract type MovingPoints <: PointMotionType end

abstract type FreestreamType end
abstract type StaticFreestream <: FreestreamType end
abstract type VariableFreestream <: FreestreamType end

abstract type FlowSide end
abstract type ExternalFlow <: FlowSide end
abstract type InternalFlow <: FlowSide end
abstract type ExternalInternalFlow <: FlowSide end


include("utils/pulses.jl")
include("navier_stokes.jl")
include("directmotion.jl")
include("plot_recipes.jl")
include("LidDrivenCavity.jl")


end
