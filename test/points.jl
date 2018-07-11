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
  H̃ = Regularize(X,dx,filter=true)


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

  @testset "Regularize vector to dual and primal nodes" begin

  f = VectorData(X)
  f.u .= rand(n)
  f.v .= rand(n)

  w = NodePair(Dual,(nx,ny))

  H(w,f)
  @test sum(f.u) ≈ sum(w.u)*dx*dx
  @test sum(f.v) ≈ sum(w.v)*dx*dx

  w = NodePair(Primal,(nx,ny))

  H(w,f)
  @test sum(f.u) ≈ sum(w.u)*dx*dx
  @test sum(f.v) ≈ sum(w.v)*dx*dx

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


  @testset "Matrix representation" begin

  f = ScalarData(X)
  w = Nodes(Dual,(nx,ny))
  Hmat = RegularizationMatrix(H,f,w)
  Emat = InterpolationMatrix(H,w,f)

  f .= rand(n)

  w2 = Nodes(Dual,(nx,ny))
  A_mul_B!(w,Hmat,f)
  H(w2,f)
  @test w ≈ w2

  @test_throws MethodError A_mul_B!(f,Hmat,w)

  f2 = ScalarData(f)
  A_mul_B!(f,Emat,w)
  H(f2,w)
  @test f ≈ f2

  w .= rand(nx,ny)
  Ẽmat = InterpolationMatrix(H̃,w,f)

  A_mul_B!(f,Ẽmat,w)
  H̃(f2,w)
  @test f ≈ f2

  w = Nodes(Primal,(nx,ny))
  Hmat = RegularizationMatrix(H,f,w)
  Emat = InterpolationMatrix(H,w,f)
  Ẽmat = InterpolationMatrix(H̃,w,f)


  w2 = Nodes(Primal,(nx,ny))
  A_mul_B!(w,Hmat,f)
  H(w2,f)
  @test w ≈ w2

  w .= rand(size(w))
  f2 = ScalarData(f)
  A_mul_B!(f,Emat,w)
  H(f2,w)
  @test f ≈ f2

  f2 = ScalarData(f)
  A_mul_B!(f,Ẽmat,w)
  H̃(f2,w)
  @test f ≈ f2

  f = VectorData(X)
  f.u .= rand(n)
  f.v .= rand(n)

  p = Edges(Dual,(nx,ny))
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)

  p2 = Edges(Dual,(nx,ny))
  A_mul_B!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(size(p.u))
  p.v .= rand(size(p.v))
  f2 = VectorData(f)
  A_mul_B!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  A_mul_B!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  p = NodePair(Dual,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)

  A_mul_B!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(size(p.u))
  p.v .= rand(size(p.v))
  A_mul_B!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  A_mul_B!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  p = NodePair(Primal,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)
  A_mul_B!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(size(p.u))
  p.v .= rand(size(p.v))
  A_mul_B!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  A_mul_B!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  end


end
