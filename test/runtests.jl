using ViscousFlow
using Test
##using TestSetExtensions
using Literate


include("pointforce.jl")
#include("navierstokes.jl")

outputdir = "../examples"
litdir = "./literate"

for (root, dirs, files) in walkdir(litdir)
    for file in files
        endswith(file,".jl") && Literate.notebook(joinpath(root, file),outputdir)
    end
end
