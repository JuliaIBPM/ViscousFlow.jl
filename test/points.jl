import Whirl: Fields
using Fields

@testset "Point-Field Routines" begin

  @testset "Point creation" begin
    @test_throws AssertionError VectorData([1,2,3],[1,2])

    f = ScalarData(10)
    f .= rand(10)

  end

  n = 10
  x = 0.5 + 0.2*rand(n)
  y = 0.5 + 0.2*rand(n)
  X = VectorData(x,y)

  nx = 12; ny = 12

  dx = 0.1
  H = Regularize(X,dx)

  @testset "Regularize vector to primal edges" begin

    f = VectorData(X)
    f.u .= rand(n)
    f.v .= rand(n)

    q = Edges(Primal,(nx,ny))

    H(q,f)
    @test sum(f.u) ≈ sum(q.u)*dx*dx
    @test sum(f.v) ≈ sum(q.v)*dx*dx

  end

  @testset "Regularize vector to dual edges" begin

  f = VectorData(X)
  f.u .= rand(n)
  f.v .= rand(n)

  p = Edges(Dual,(nx,ny))

  H(p,f)
  @test sum(f.u) ≈ sum(p.u)*dx*dx
  @test sum(f.v) ≈ sum(p.v)*dx*dx

  end

  @testset "Regularize scalar to primal nodes" begin

  f = ScalarData(X)
  f .= rand(n)

  w = Nodes(Primal,(nx,ny))

  H(w,f)
  @test sum(f) ≈ sum(w)*dx*dx

  end

  @testset "Regularize scalar to dual nodes" begin

  f = ScalarData(X)
  f .= rand(n)

  w = Nodes(Dual,(nx,ny))

  H(w,f)
  @test sum(f) ≈ sum(w)*dx*dx

  end



end
