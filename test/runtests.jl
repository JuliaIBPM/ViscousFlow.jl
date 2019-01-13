using Compat.Test
using Compat

#using Test
##using TestSetExtensions
##using Distributed

using Whirl
#import Whirl

#@test isempty(detect_ambiguities(Whirl))

include("fields.jl")
#include("points.jl")
#include("timemarching.jl")
#include("saddle.jl")
#include("systems.jl")


#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
