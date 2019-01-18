using Compat.Test
using Compat

#using Test
##using TestSetExtensions

using ViscousFlow

#@test isempty(detect_ambiguities(ViscousFlow))

include("fields.jl")
include("points.jl")
include("timemarching.jl")
include("saddle.jl")
include("systems.jl")


#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
