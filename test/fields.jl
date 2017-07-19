include(joinpath(Pkg.dir("Whirl2d"), "src/Fields.jl"))
include(joinpath(Pkg.dir("Whirl2d"), "src/Grids.jl"))
using Fields

@testset "Fields" begin
    @testset "Shifting Primal Edges to Dual Edges" begin

        edges = Edges(Primal, (3, 5), 1)
        edges.u .= reshape(1:30, 5, 6)
        edges.v .= reshape(1:28, 4, 7)

        dual = Fields.shift(edges)
        @test dual.u == [0.0  4.0   9.0  14.0  19.0  24.0  0.0
                         0.0  5.0  10.0  15.0  20.0  25.0  0.0
                         0.0  6.0  11.0  16.0  21.0  26.0  0.0
                         0.0  7.0  12.0  17.0  22.0  27.0  0.0]

        @test dual.v == [0.0  0.0   0.0   0.0   0.0   0.0
                         3.5  7.5  11.5  15.5  19.5  23.5
                         4.5  8.5  12.5  16.5  20.5  24.5
                         5.5  9.5  13.5  17.5  21.5  25.5
                         0.0  0.0   0.0   0.0   0.0   0.0]

        edges.u .= reshape(rand(30), 5, 6)
        edges.v .= reshape(rand(28), 4, 7)

        vx = zeros(dual.u)
        vy = zeros(dual.v)
        Grids.shift!(vx, vy, 1+(1:3), 1+(1:5), edges.u, edges.v)
        Fields.shift!(dual, edges)
        @test dual.u ≈ vx
        @test dual.v ≈ vy
    end

    @testset "Shifting Dual Nodes to Dual Edges" begin

        ω = Nodes(Dual, (3, 5), 1)
        ω .= reshape(1:35, 5, 7)

        dual = shift(ω)

        @test dual.u == [0.0  6.5  11.5  16.5  21.5  26.5  0.0
                         0.0  7.5  12.5  17.5  22.5  27.5  0.0
                         0.0  8.5  13.5  18.5  23.5  28.5  0.0
                         0.0  9.5  14.5  19.5  24.5  29.5  0.0]
        @test dual.v == [0.0   0.0   0.0   0.0   0.0   0.0
                         4.5   9.5  14.5  19.5  24.5  29.5
                         5.5  10.5  15.5  20.5  25.5  30.5
                         6.5  11.5  16.5  21.5  26.5  31.5
                         0.0   0.0   0.0   0.0   0.0   0.0]

        ω .= rand(size(ω))
        vx = zeros(dual.u)
        vy = zeros(dual.v)

        Grids.shift!(vx, vy, 1+(1:3), 1+(1:5), ω.data)
        Fields.shift!(dual, ω)
        @test dual.u ≈ vx
        @test dual.v ≈ vy
    end
end
