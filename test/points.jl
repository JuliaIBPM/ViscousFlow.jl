import ViscousFlow: Fields
using Compat
using Compat.LinearAlgebra

#if VERSION < v"0.7-"
#  import Base: A_mul_B!
#  mul!(x,B,y) = A_mul_B!(x,B,y)
#end

@testset "Point-Field Routines" begin

  @testset "Point creation" begin
    @test_throws AssertionError VectorData([1,2,3],[1,2])

    @test_throws AssertionError TensorData([1,2,3],[1,2],[2,3],[4,5])

    f = ScalarData(10)
    f .= rand(10)

    ft = TensorData(f)
    ft[25] = 4
    @test ft.dvdx[5] == 4

  end

  @testset "Point operations" begin

    Y = VectorData(4)
    X = Y + (1,2)
    @test X.v[4] == 2.0

    X2 = (1,2) + Y
    @test X2 == X

    X = Y - (1,2)
    @test X.u[4] == -1.0

    X = (1,2) - Y
    @test X.u[4] == 1.0

    Z = 2.0 × X
    @test typeof(Z) <: VectorData
    @test Z.u[1] == -4.0

    X = TensorData(Y)
    fill!(X,1)
    Z = (1,2) ⋅ X
    @test typeof(Z) <: VectorData
    @test Z.u[3] == 3.0

  end

  n = 10
  x = 0.5 .+ 0.2*rand(n)
  y = 0.5 .+ 0.2*rand(n)
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

  @testset "Regularize tensor to edge gradient" begin

  f = TensorData(X)
  f.dudx .= rand(n)
  f.dudy .= rand(n)
  f.dvdx .= rand(n)
  f.dvdy .= rand(n)

  gradq = EdgeGradient(Dual,(nx,ny))

  H(gradq,f)

  @test sum(f.dudx) ≈ sum(gradq.dudx)*dx*dx
  @test sum(f.dudy) ≈ sum(gradq.dudy)*dx*dx
  @test sum(f.dvdx) ≈ sum(gradq.dvdx)*dx*dx
  @test sum(f.dvdy) ≈ sum(gradq.dvdy)*dx*dx

  gradq = EdgeGradient(Primal,(nx,ny))

  H(gradq,f)

  @test sum(f.dudx) ≈ sum(gradq.dudx)*dx*dx
  @test sum(f.dudy) ≈ sum(gradq.dudy)*dx*dx
  @test sum(f.dvdx) ≈ sum(gradq.dvdx)*dx*dx
  @test sum(f.dvdy) ≈ sum(gradq.dvdy)*dx*dx


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
  mul!(w,Hmat,f)
  H(w2,f)
  @test w ≈ w2

  @test_throws MethodError mul!(f,Hmat,w)

  f2 = ScalarData(f)
  mul!(f,Emat,w)
  H(f2,w)
  @test f ≈ f2

  w .= rand(nx,ny)
  Ẽmat = InterpolationMatrix(H̃,w,f)

  mul!(f,Ẽmat,w)
  H̃(f2,w)
  @test f ≈ f2

  w = Nodes(Primal,(nx,ny))
  Hmat = RegularizationMatrix(H,f,w)
  Emat = InterpolationMatrix(H,w,f)
  Ẽmat = InterpolationMatrix(H̃,w,f)


  w2 = Nodes(Primal,(nx,ny))
  mul!(w,Hmat,f)
  H(w2,f)
  @test w ≈ w2

  w .= rand(Float64,size(w))
  f2 = ScalarData(f)
  mul!(f,Emat,w)
  H(f2,w)
  @test f ≈ f2

  f2 = ScalarData(f)
  mul!(f,Ẽmat,w)
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
  mul!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(Float64,size(p.u))
  p.v .= rand(Float64,size(p.v))
  f2 = VectorData(f)
  mul!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  mul!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  p = NodePair(Dual,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)

  mul!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(Float64,size(p.u))
  p.v .= rand(Float64,size(p.v))
  mul!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  mul!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  p = NodePair(Primal,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)
  mul!(p,Hmat,f)
  H(p2,f)
  @test p.u ≈ p2.u && p.v ≈ p2.v

  p.u .= rand(Float64,size(p.u))
  p.v .= rand(Float64,size(p.v))
  mul!(f,Emat,p)
  H(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  mul!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.u ≈ f2.u && f.v ≈ f2.v

  f = TensorData(X)
  f.dudx .= rand(n)
  f.dudy .= rand(n)
  f.dvdx .= rand(n)
  f.dvdy .= rand(n)
  f2 = TensorData(f)

  p = EdgeGradient(Dual,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)

  mul!(p,Hmat,f)
  H(p2,f)
  @test p.dudx ≈ p2.dudx && p.dudy ≈ p2.dudy && p.dvdx ≈ p2.dvdx && p.dvdy ≈ p2.dvdy

  p.dudx .= rand(Float64,size(p.dudx))
  p.dudy .= rand(Float64,size(p.dudy))
  p.dvdx .= rand(Float64,size(p.dvdx))
  p.dvdy .= rand(Float64,size(p.dvdy))
  mul!(f,Emat,p)
  H(f2,p)
  @test f.dudx ≈ f2.dudx && f.dudy ≈ f2.dudy && f.dvdx ≈ f2.dvdx && f.dvdy ≈ f2.dvdy

  mul!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.dudx ≈ f2.dudx && f.dudy ≈ f2.dudy && f.dvdx ≈ f2.dvdx && f.dvdy ≈ f2.dvdy

  p = EdgeGradient(Primal,(nx,ny))
  p2 = deepcopy(p)
  Hmat = RegularizationMatrix(H,f,p)
  Emat = InterpolationMatrix(H,p,f)
  Ẽmat = InterpolationMatrix(H̃,p,f)

  mul!(p,Hmat,f)
  H(p2,f)
  @test p.dudx ≈ p2.dudx && p.dudy ≈ p2.dudy && p.dvdx ≈ p2.dvdx && p.dvdy ≈ p2.dvdy

  p.dudx .= rand(Float64,size(p.dudx))
  p.dudy .= rand(Float64,size(p.dudy))
  p.dvdx .= rand(Float64,size(p.dvdx))
  p.dvdy .= rand(Float64,size(p.dvdy))
  mul!(f,Emat,p)
  H(f2,p)
  @test f.dudx ≈ f2.dudx && f.dudy ≈ f2.dudy && f.dvdx ≈ f2.dvdx && f.dvdy ≈ f2.dvdy

  mul!(f,Ẽmat,p)
  H̃(f2,p)
  @test f.dudx ≈ f2.dudx && f.dudy ≈ f2.dudy && f.dvdx ≈ f2.dvdx && f.dvdy ≈ f2.dvdy

  end


end
