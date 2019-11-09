#using Compat.Test
#using Compat

using ViscousFlow
using Test
##using TestSetExtensions


#@test isempty(detect_ambiguities(ViscousFlow))

include("fields.jl")
include("points.jl")
include("layers.jl")
include("timemarching.jl")
include("saddle.jl")
include("systems.jl")


#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
