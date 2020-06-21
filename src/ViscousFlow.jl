"""
The `ViscousFlow` module is here

"""
module ViscousFlow

using Reexport

@reexport using CartesianGrids
@reexport using RigidBodyTools
@reexport using ConstrainedSystems

#== Imports/Exports ==#

include("utils.jl")
include("systems.jl")


#== Plot Recipes ==#

include("plot_recipes.jl")



end
