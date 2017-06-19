"""
The `Whirl2d` module is here

"""
module Whirl2d

#== Imports/Exports ==#

include("Utils.jl")
using .Utils


#== Type Definitions ==#

const ndim = 2

include("Grids.jl")
using .Grids

include("Bodies.jl")
using .Bodies

include("ddf.jl")
using .DDF

include("Systems.jl")
using .Systems

mutable struct Soln

  # Current time of solution
  t::Float64

  # Domain structure
  dom::Systems.Domain

  # Grid vorticity vector
  w::Array{Float64,Whirl2d.ndim}

  # Body force vector
  f::Array{Float64,2}

end

function Soln(dom)

  t = 0.0
  w = zeros(dom.grid.cell)
  f = zeros(dom.nbodypts,Whirl2d.ndim)

  Soln(t,dom,w,f)
end

function Base.show(io::IO, s::Soln)
    println(io, "Solution: t = $(s.t)")
end




end
