"""
The `Whirl2d` module is here

"""
module Whirl2d

#== Imports/Exports ==#

include("Utils.jl")
using .Utils


#== Type Definitions ==#

const ndim = 2

abstract type Soln end

mutable struct ConstrainedSoln<:Soln
  # Current time of solution
  t::Float64

  # Solution data
  u

  # Constraint data
  f

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


function Base.show(io::IO, s::ConstrainedSoln)
    println(io, "Solution: t = $(s.t)")
end




end
