"""
The `Whirl2d` module is here

"""
module Whirl2d

#== Imports/Exports ==#

include("Utils.jl")
using .Utils


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


include("TimeMarching.jl")
using .TimeMarching

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

# include("RigidBodyMotions.jl")
# using .RigidBodyMotions


function Base.show(io::IO, s::SolnType)
    println(io, "Solution: t = $(s.t)")
end




end
