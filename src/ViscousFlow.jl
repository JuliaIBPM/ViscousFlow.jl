module ViscousFlow

#using DocStringExtensions
using Reexport
using UnPack
@reexport using ImmersedLayers
@reexport using GridUtilities


export ViscousIncompressibleFlowProblem
export setup_grid, viscousflow_system, setup_problem, surface_point_spacing,
        surface_velocity_in_translating_frame!

#= Supporting functions =#


include("defaults.jl")
include("reference_frame.jl")
include("fields.jl")
include("ode_operators.jl")
include("api.jl")









end
