import ViscousFlow: SaddlePointSystems

using LinearAlgebra

@testset "Saddle-Point Systems" begin

  @testset "Matrix tests" begin

    A1 = Float64[1 2; 2 1]
    B2 = Float64[2 3;-1 -1]
    B1 = B2'
    C = Matrix{Float64}(undef,2,2)
    C .= [5 -2; 3 -4];

    A = SaddleSystem(A1,B2,B1,C)

    # test inputs and outputs as tuples
    rhs = ([1.0,2.0],[3.0,4.0])
    sol = (zeros(2),zeros(2))

    ldiv!(sol,A,rhs)

    @test norm((A*sol)[1]-rhs[1]) < 1e-14
    @test norm((A*sol)[2]-rhs[2]) < 1e-14

    sol2 = A\rhs

    @test norm((A*sol2)[1]-rhs[1]) < 1e-14
    @test norm((A*sol2)[2]-rhs[2]) < 1e-14

    # test inputs and outputs as vectors
    rhs1, rhs2 = rhs
    rhsvec = [rhs1;rhs2]
    solvec = similar(rhsvec)

    solvec = A\rhsvec

    @test norm(A*solvec-rhsvec) < 1e-14

    # test with the other constructor
    Aother = SaddleSystem(A1,B2,B1,C,zeros(Float64,2),zeros(Float64,2))

    sol_other = Aother\rhs

    @test norm((Aother*sol_other)[1]-rhs[1]) < 1e-14
    @test norm((Aother*sol_other)[2]-rhs[2]) < 1e-14

    # Test with SaddleVector
    rhs = SaddleVector(rhs1,rhs2)
    As = SaddleSystem(A1,B2,B1,C,rhs)

    sol = As\rhs

    @test norm(state(As*sol)-state(rhs)) < 1e-14
    @test norm(constraint(As*sol)-constraint(rhs)) < 1e-14

  end

  nx = 130; ny = 130
  Lx = 2.0
  dx = Lx/(nx-2)
  w = Nodes(Dual,(nx,ny))

  L = plan_laplacian(size(w),with_inverse=true)

  n = 128
  θ = range(0,stop=2π,length=n+1)
  R = 0.5
  xb = 1.0 .+ R*cos.(θ[1:n])
  yb = 1.0 .+ R*sin.(θ[1:n])
  ds = (2π/n)*R
  X = VectorData(xb,yb)
  f = ScalarData(X)

  E = Regularize(X,dx;issymmetric=true)
  Hmat,Emat = RegularizationMatrix(E,f,w)

  @testset "Construction of linear maps" begin

    u = similar(w)
    wvec = vec(w)

    w .= rand(size(w)...)

    uvec = zeros(length(w))

    Lop = SaddlePointSystems.linear_map(L,w)

    u = L*w

    uvec = Lop*wvec

    @test SaddlePointSystems._wrap_vec(uvec,u) == u

    Linv = SaddlePointSystems.linear_inverse_map(L,w)

    yvec = Linv*wvec

    y = L\w

    @test SaddlePointSystems._wrap_vec(yvec,y) == y

    # point-wise operators
    fvec = vec(f)

    f[10] = 1.0
    Hop = SaddlePointSystems.linear_map(Hmat,f,w);

    y = Hmat*f

    yvec = Hop*fvec

    @test SaddlePointSystems._wrap_vec(yvec,y) == y

    Eop = SaddlePointSystems.linear_map(Emat,w,f)

    g = Emat*w

    gvec = Eop*vec(w)

    @test SaddlePointSystems._wrap_vec(gvec,g) == g

  end

  ψb = ScalarData(X)
  w = Nodes(Dual,(nx,ny))
  ψb .= -(xb .- 1)
  f .= ones(Float64,n)*ds
  ψ = Nodes(Dual,w)

  nada = empty(f)

  @testset "Field operators" begin

    rhs = SaddleVector(w,ψb)

    A = SaddleSystem(L,Emat,Hmat,rhs)

    A.S*vec(f)

    A.S⁻¹*vec(f)

    sol = deepcopy(rhs)

    ldiv!(sol,A,rhs)

    sol2 = A\rhs

    @test state(sol2) == state(sol)
    @test constraint(sol2) == constraint(sol)

    rhs2 = A*sol2

    @test norm(constraint(rhs2)-ψb) < 1e-14

    fex = -2*cos.(θ[1:n])
    @test norm(constraint(sol)-fex*ds) < 0.02
    @test ψ[nx,65] ≈ -ψ[1,65]
    @test ψ[65,ny] ≈ ψ[65,1]


  end

  fv = VectorData(X)
  q = Edges(Primal,w)
  Hvmat,Evmat = RegularizationMatrix(E,fv,q)

  @testset "Vector force data" begin

    B₁ᵀ(f) = Curl()*(Hvmat*f)
    B₂(w) = -(Evmat*(Curl()*(L\w)))

    rhsv = SaddleVector(w,fv)
    A = SaddleSystem(I,B₂,B₁ᵀ,rhsv)

    sol = A\rhsv
    @test state(sol) == zero(w) && constraint(sol) == zero(fv)

  end

  Ẽ = Regularize(X,dx;weights=ds,filter=true)
  H̃mat = RegularizationMatrix(Ẽ,ψb,w)
  Ẽmat = InterpolationMatrix(Ẽ,w,ψb)
  Pmat = Ẽmat*H̃mat

  @testset "Filtering" begin

    rhs = SaddleVector(w,ψb)
    Afilt = SaddleSystem(L,Emat,Hmat,rhs,filter=Pmat)

    sol = Afilt\rhs

    fex = -2*cos.(θ[1:n])
    @test norm(constraint(sol)-fex*ds) < 0.01

  end

  @testset "Reduction to unconstrained system" begin

    op = SaddlePointSystems.linear_map(nothing,nada)

    @test size(op) == (0,0)
    @test op*nada == ()

    op = SaddlePointSystems.linear_map(nothing,w,nada)

    @test size(op) == (0,length(w))
    @test op*vec(w) == ()

    op = SaddlePointSystems.linear_map(nothing,nada,w)

    @test size(op) == (length(w),0)
    @test op*nada == vec(zero(w))

    rhsnc = SaddleVector(w,nada)

    Anc = SaddleSystem(L,rhsnc)

    sol = Anc\rhsnc
    ψ = state(sol)
    fnull = constraint(sol)

    @test ψ == L\w
    @test fnull == nada

    q.u .= 3;
    rhsqnc = SaddleVector(q,nada)
    Aqnc = SaddleSystem(I,rhsqnc)

    sol = Aqnc\rhsqnc

    q2 = state(sol)
    fnull = constraint(sol)

    @test all(q2.u .== 3.0) && all(q2.v .== 0.0)

  end

  @testset "Tuple of saddle point systems" begin

    rhs = SaddleVector(w,ψb)
    A = SaddleSystem(L,Emat,Hmat,rhs)
    sol = A\rhs
    ψ = state(sol)
    f = constraint(sol)

    rhsnc = SaddleVector(w,nada)
    Anc = SaddleSystem(L,rhsnc)

    sys = (A,Anc)
    sol1, sol2 = sys\(rhs,rhsnc)


    @test state(sol1) == ψ && constraint(sol1) == f && state(sol2) == zero(w) && constraint(sol2) == nada

    newrhs1,newrhs2 = sys*(sol1,sol2)

    @test newrhs1 == A*sol1 && newrhs2 == Anc*sol2

  end

  @testset "Recursive saddle point with vectors" begin
    A1 = Float64[1 2; 2 1]
    B21 = Float64[2 3]
    B11 = B21'
    C1 = Matrix{Float64}(undef,1,1)
    C1.= 5

    B22 = Float64[-1 -1 3]
    B12 = Float64[-1 -1 -2]'
    C2 = Matrix{Float64}(undef,1,1)
    C2.= -4

    rhs11 = [1.0,2.0]
    rhs12 = Vector{Float64}(undef,1)
    rhs12 .= 3.0
    #rhs1 = (rhs11,rhs12)
    rhs1 = [rhs11;rhs12]

    rhs2 = Vector{Float64}(undef,1)
    rhs2 .= 4.0

    rhs = (rhs1,rhs2)

    A = SaddleSystem(A1,B21,B11,C1)

    Abig = SaddleSystem(A,B22,B12,C2,rhs1,rhs2)

    sol = Abig\rhs

    out = Abig*sol

    @test norm(out[1]-rhs1) < 1e-14


  end


end
