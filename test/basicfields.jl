using ImmersedLayers
using LinearAlgebra

@testset "Basic fields" begin

  my_params = Dict()
  xlim = (-2.0,2.0)
  ylim = (-2.0,2.0)
  my_params["Re"] = 200
  my_params["grid Re"] = 4.0
  g = setup_grid(xlim,ylim,my_params)

  sys = viscousflow_system(g,phys_params=my_params)

  σ = 0.2
  x0 = 0.0
  y0 = 0.0
  A = 1
  gauss = SpatialGaussian(σ,σ,x0,y0,A)

  u = init_sol(gauss,sys)

  t = 0.0
  w = vorticity(u,sys,t)
  vel = velocity(u,sys,t)
  vdv = convective_acceleration(u,sys,t)
  p = pressure(u,sys,t)
  ψ = streamfunction(u,sys,t)
  Q = Qcrit(u,sys,t)
  
  tspan = (0.0,10.0)
  integrator = init(u,tspan,sys)
  step!(integrator,1.0)

  oseen_exact(t) = SpatialGaussian(sqrt(σ^2+2*t/my_params["Re"]),sqrt(σ^2+2*t/my_params["Re"]),x0,y0,A)
  exactsol(t) = init_sol(oseen_exact(t),sys)

  @test norm(vorticity(exactsol(integrator.t),sys,integrator.t)-vorticity(integrator),sys) < 0.005

end

@testset "Inviscid system" begin

  xlim = (-2.0,2.0)
  ylim = (-2.0,2.0)
  g = setup_grid(xlim,ylim,0.02)

  sys = viscousflow_system(g)

  σ = 0.2
  x0 = 0.0
  y0 = 0.0
  A = 1
  gauss = SpatialGaussian(σ,σ,x0,y0,A)
  u = init_sol(gauss,sys)



end
