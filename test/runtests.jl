using Test
#using TestSetExtensions
#using Distributed

using Whirl
#import Whirl

#@test isempty(detect_ambiguities(Whirl))

include("fields.jl")
#include("Polygons.jl")

#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#include("fields.jl")
#include("gaussian_vortex.jl")

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
