using Base.Test
using TestSetExtensions

using Whirl

@test isempty(detect_ambiguities(Whirl))

@testset ExtendedTestSet "All tests" begin
    @includetests ARGS
end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
