include(joinpath(Pkg.dir("Whirl2d"), "src/Fields.jl"))
using Fields

@testset "Fields" begin
    @testset "Hadamard Product" begin
        edges_p  = Edges{Primal, 30, 40}()
        edges_p.u .= rand(size(edges_p.u))

        # Should be safe for the output to be the same as the input
        edges_p2 = edges_p ∘ edges_p
        product!(edges_p, edges_p, edges_p)
        @test edges_p2.u == edges_p.u
        @test edges_p2.v == edges_p.v

        edges_d  = Edges{Dual, 30, 40}()
        @test_throws MethodError (edges_p ∘ edges_d)
    end

    @testset "Discrete Laplacian" begin
        s = DualNodes{30, 40}()
        s[3:end-2, 3:end-2] .= rand(26, 36)

        L = Laplacian(30, 40)

        @test L*s ≈ -curl(curl(s))

        @test_throws MethodError (L \ s)

        L = Laplacian(30, 40, with_inverse = true, fftw_flags = FFTW.PATIENT)
        @test L \ (L*s) ≈ s
    end

    @testset "Discrete Divergence" begin
        s = DualNodes{5, 4}()
        s .= rand(5, 4)

        @test iszero(divergence(Fields.shift(curl(s))))

        q′ = Edges{Dual, 5, 4}()
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
        s = DualNodes{5, 4}()
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

        q = Edges{Primal, 5, 4}()
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

        w = DualNodes{5, 4}()
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
