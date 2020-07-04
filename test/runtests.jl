#using Compat.Test
#using Compat

using ViscousFlow
using Test
##using TestSetExtensions


#@test isempty(detect_ambiguities(ViscousFlow))
include("pointforce.jl")
include("navierstokes.jl")


#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
