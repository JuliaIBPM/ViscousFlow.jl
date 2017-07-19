include(joinpath(Pkg.dir("Whirl2d"), "src/Fields.jl"))
include(joinpath(Pkg.dir("Whirl2d"), "src/Grids.jl"))
using Fields

@testset "Fields" begin
    @testset "Discrete Divergence" begin
        s = DualNodes(5, 4)
        s .= rand(5, 4)

        @test iszero(divergence(Fields.shift(curl(s))))

        q′ = Edges(Dual, (5, 4))
        q′.u .= reshape(1:8, 4, 2)
        q′.v .= reshape(1:9, 3, 3)

        # Not sure if this is the behavior we want yet
        # Currently, the ghost cells are not affected
        # by the divergence operator
        s .= 1.0
        divergence!(s, q′)
        @test s == [ 1.0  1.0  1.0  1.0
                     1.0  4.0  4.0  1.0
                     1.0  4.0  4.0  1.0
                     1.0  4.0  4.0  1.0
                     1.0  1.0  1.0  1.0 ]
    end

    @testset "Discrete Curl" begin
        s = DualNodes(5, 4)
        s .= reshape(1:20, 4, 5)'

        q = curl(s)

        @test q.u == [ 1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0
                       1.0  1.0  1.0 ]

        @test q.v == [ -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0
                       -4.0  -4.0  -4.0  -4.0 ]
    end

    @testset "Shifting Primal Edges to Dual Edges" begin

        q = Edges(Primal, (5, 4))
        q.u .= reshape(1:15, 5, 3)
        q.v .= reshape(1:16, 4, 4)

        Qq = Fields.shift(q)
        @test Qq.u == [ 4.0   9.0
                        5.0  10.0
                        6.0  11.0
                        7.0  12.0 ]

        @test Qq.v == [ 3.5  7.5  11.5
                        4.5  8.5  12.5
                        5.5  9.5  13.5 ]

    end

    @testset "Shifting Dual Nodes to Dual Edges" begin

        w = DualNodes(5, 4)
        w .= reshape(1:20, 5, 4)

        Ww = Fields.shift(w)

        @test Ww.u == [ 6.5  11.5
                        7.5  12.5
                        8.5  13.5
                        9.5  14.5 ]
        @test Ww.v == [ 4.5   9.5  14.5
                        5.5  10.5  15.5
                        6.5  11.5  16.5 ]
    end
end
