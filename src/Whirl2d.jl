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

include("Systems.jl")
using .Systems

include("ddf.jl")
using .DDF



end
