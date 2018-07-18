"""
The `Whirl` module is here

"""
module Whirl

using Reexport


#== Imports/Exports ==#

include("Utils.jl")
@reexport using .Utils

include("Fields.jl")

@reexport using .Fields

include("SaddlePointSystems.jl")

@reexport using .SaddlePointSystems

include("TimeMarching.jl")
@reexport using .TimeMarching

include("Systems.jl")
@reexport using .Systems

#== Type Definitions ==#

const ndim = 2

abstract type SolnType end

mutable struct Soln{T} <: SolnType
  # Current time of solution
  t::Float64

  # Solution data
  u::T

end

Soln(u::T)  where {T} = Soln{T}(0.0,u)

mutable struct ConstrainedSoln{T,K} <: SolnType
  # Current time of solution
  t::Float64

  # Solution data
  u::T

  # Constraint data
  f::K

end

ConstrainedSoln(u::T,f::K) where {T,K} = ConstrainedSoln{T,K}(0.0,u,f)

#ConstrainedSoln(u::Vector{Array{T,2}},f::Vector{Array{T,2}}) where {T} = ConstrainedSoln{Vector{Array{T,2}},Vector{Array{T,2}}}(0.0,u,f)

#=
include("Fields.jl")
import .Fields

include("IntFactSystems.jl")
import .IntFactSystems



include("Process.jl")
using .Process

include("Grids.jl")
using .Grids

include("Bodies.jl")
using .Bodies

include("ddf.jl")
using .DDF

include("Systems.jl")
using .Systems

include("NavierStokes.jl")
using .NavierStokes

=#

# include("RigidBodyMotions.jl")
# using .RigidBodyMotions


#== Plot Recipes ==#

include("plot_recipes.jl")

function Base.show(io::IO, s::SolnType)
    println(io, "Solution: t = $(s.t)")
end




end
