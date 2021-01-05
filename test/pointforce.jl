using LinearAlgebra
import ViscousFlow: DGaussian, SpatialDGaussian

@testset "PointForce" begin

  Re = 200
  xlim, ylim = (-2,2), (-2,4)
  Δx,Δt = setstepsizes(Re,gridRe=4)

  sys = NavierStokes(Re,Δx,xlim,ylim,Δt)
  w₀ = Nodes(Dual,size(sys))

  x0 = [1.5,0.0]
  f0 = 10.0
  t0 = 4
  σt = 1
  pulse = PointForce(w₀,x0,f0,t0,σt,sys)

  @test isapprox(abs(sum(pulse(t0))-f0),0,atol=1e-2)

  x0 = (1.5,0.0)
  pulse = PointForce(w₀,x0,f0,t0,σt,sys)

  @test isapprox(abs(sum(pulse(t0))-f0),0,atol=1e-2)


end

@testset "Pulses" begin

  σ = 0.2
  x0 = 0
  A = 1
  dg = DGaussian(σ,x0,A)
  x = 0.1
  @test dg(x) == -2*A*x/sqrt(π)/σ^3*exp(-x^2/σ^2)

  σx = 0.5
  σy = 0.5
  x0 = 0
  y0 = 0
  A = 1
  dgaussx = SpatialDGaussian(σx,σy,x0,y0,A,1)
  dgaussy = SpatialDGaussian(σx,σy,x0,y0,A,2)
  x = 0.1
  y = 0.2
  @test dgaussx(0.1,0.2) == dgaussy(0.2,0.1) ≈ -2*A*x/π/σx^3/σy*exp(-x^2/σx^2)*exp(-y^2/σy^2)

end
