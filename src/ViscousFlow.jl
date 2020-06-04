"""
The `ViscousFlow` module is here

"""
module ViscousFlow

using Reexport


#== Imports/Exports ==#

include("Utils.jl")
@reexport using .Utils

@reexport using CartesianGrids

include("Bodies.jl")

@reexport using .Bodies

include("RigidBodyMotions.jl")

@reexport using .RigidBodyMotions

include("Layers.jl")

@reexport using .Layers

include("SaddlePointSystems.jl")

@reexport using .SaddlePointSystems

include("TimeMarching.jl")
@reexport using .TimeMarching

include("Systems.jl")
@reexport using .Systems



#== Plot Recipes ==#

include("plot_recipes.jl")



end
