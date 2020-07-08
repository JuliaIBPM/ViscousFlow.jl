using LinearAlgebra

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
