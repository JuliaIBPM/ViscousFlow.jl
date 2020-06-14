import ViscousFlow: Layers

Δx = 0.02
xlim = (-5.98,5.98)
ylim = (-5.98,5.98)
g = PhysicalGrid(xlim,ylim,Δx)

w = Nodes(Dual,size(g))
dq = Edges(Dual,w)
n = 150
radius = 1.0
body = Circle(radius,n)
X = VectorData(collect(body))
f = VectorData(X)
ϕ = ScalarData(f)
fill!(ϕ,1.0)

regop = Regularize(X,Δx,I0=origin(g),issymmetric=true,ddftype=CartesianGrids.Yang3)
Hv, Ev = RegularizationMatrix(regop,f,dq)
Hs, Es = RegularizationMatrix(regop,ϕ,w)

@testset "Double Layer" begin

  dlayer = DoubleLayer(body,Hv,weight=1/Δx)
  ϕ .= rand(length(ϕ))

  @test abs(sum(dlayer(ϕ))) < 100*eps(1.0)

  dlayer2 = DoubleLayer(body,regop,w)

  @test dlayer2(ϕ) == dlayer(ϕ)

end

@testset "Single Layer" begin

  slayer = SingleLayer(body,Hs)
  ϕ .= 1.0

  @test abs(sum(slayer(ϕ))-2π*radius) < 2e-3

  slayer2 = SingleLayer(body,regop,w)

  @test slayer2(ϕ) == slayer(ϕ)

end

@testset "Masks" begin

  inner = Mask(body,regop,w)

  fill!(w,1)

  @test abs(sum(inner(w))*cellsize(g)^2 - π*radius^2) < 1e-3

  outer = ComplementaryMask(inner)

end
