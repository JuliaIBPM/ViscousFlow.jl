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

  # Auxiliary solution data
  ψ::T

end

function Soln(u::Array{T,2})  where {T}
  ψ = zeros(u)
  Soln{Array{T,2}}(0.0,u,ψ)
end

function Soln(u::Vector{Array{T,2}}) where {T}
  ψ = zeros(u)
  Soln{Vector{Array{T,2}}}(0.0,u,ψ)
end

mutable struct ConstrainedSoln{T,K} <: SolnType
  # Current time of solution
  t::Float64

  # Solution data
  u::T

  # Constraint data
  f::K

  # Auxiliary solution data
  ψ::T

end

function ConstrainedSoln(u::Array{T,2},f::Array{T,2}) where {T}
  ψ = zeros(u)
  ConstrainedSoln{Array{T,2},Array{T,2}}(0.0,u,f,ψ)
end

function ConstrainedSoln(u::Vector{Array{T,2}},f::Vector{Array{T,2}}) where {T}
  ψ = [zeros(u[i]) for i=1:length(u)]
  ConstrainedSoln{Vector{Array{T,2}},Vector{Array{T,2}}}(0.0,u,f,ψ)
end

include("TimeMarching.jl")
using .TimeMarching

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


function Base.show(io::IO, s::SolnType)
    println(io, "Solution: t = $(s.t)")
end




end
