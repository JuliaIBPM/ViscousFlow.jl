include(joinpath(Pkg.dir("Whirl2d"), "src/Fields.jl"))
include(joinpath(Pkg.dir("Whirl2d"), "src/Grids.jl"))
using Fields

@testset "Fields" begin
    @testset "Discrete Curl" begin
        # Dual Node to Edge

        s = Nodes(Dual, (3, 5), 1)
        s.data .= reshape((1:35).^2, 5, 7)

        f = Edges(Primal, (3, 5), 1)
        Fields.curl!(f, s)

        @test f.u == [ 0.0    0.0    0.0    0.0    0.0    0.0
                       0.0   95.0  145.0  195.0  245.0  295.0
                       0.0  105.0  155.0  205.0  255.0  305.0
                       0.0  115.0  165.0  215.0  265.0  315.0
                       0.0  125.0  175.0  225.0  275.0  325.0]

        @test f.v == [ 0.0    0.0    0.0    0.0    0.0    0.0    0.0
                       0.0  -15.0  -25.0  -35.0  -45.0  -55.0  -65.0
                       0.0  -17.0  -27.0  -37.0  -47.0  -57.0  -67.0
                       0.0  -19.0  -29.0  -39.0  -49.0  -59.0  -69.0]

        vx = zeros(f.u)
        vy = zeros(f.v)

        s.data .= reshape(rand(35), 5, 7)
        f.u .= 0
        f.v .= 0
        Fields.curl!(f, s)
        Grids.curl!(vx, vy, 1+(1:3), 1+(1:5), s.data)
        @test vx ≈ f.u
        @test vy ≈ f.v

        # Edge to Dual Node
        f.u .= reshape((1:30).^2, 5, 6)
        f.v .= reshape((1:28).^2, 4, 7)
        s .= 0

        Fields.curl!(s, f)

        @test s == [0.0    0.0    0.0     0.0     0.0     0.0  0.0
                    0.0  -34.0  -76.0  -118.0  -160.0  -202.0  0.0
                    0.0  -42.0  -84.0  -126.0  -168.0  -210.0  0.0
                    0.0  -50.0  -92.0  -134.0  -176.0  -218.0  0.0
                    0.0    0.0    0.0     0.0     0.0     0.0  0.0]

        f.u .= reshape((1:30).^2, 5, 6)
        f.v .= reshape((1:28).^2, 4, 7)
        s .= 0
        cell = zeros(s.data)

        Grids.curl!(cell, 1+(1:3), 1+(1:5), f.u, f.v)
        Fields.curl!(s, f)
        @test s == cell
    end
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

        dual = Fields.shift(ω)

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
