
using LinearAlgebra
using SpecialFunctions

@testset "Navier-Stokes" begin

  @testset "Nonlinear convective term" begin

  xlim = (-2,2)
  ylim = (-2,2)

  Re = 100.0

  Δx, Δt = setstepsizes(Re,gridRe=1.0)
  sys = NavierStokes(Re,Δx,xlim,ylim,Δt)
  w = Nodes(Dual,size(sys))
  dw = Nodes(Dual,w)

  σx = 0.15; σy = 0.3
  gauss = SpatialGaussian(σx,σy,0,0,1.0)
  gfw = GeneratedField(w,gauss,sys.grid)
  xg, yg = coordinates(w,sys.grid)
  w .= cellsize(sys)*gfw().*2.0.*(1/σx^2 .+ 1/σy^2 .- 2.0.*(xg.^2/σx^4 .+ (yg').^2/σy^4))

  ns_rhs!(dw,w,sys,0.0)

  nlex = Nodes(Dual,w)
  gauss2 = SpatialGaussian(σx/sqrt(2),σy/sqrt(2),0,0,1.0)
  gfnl = GeneratedField(w,gauss2,sys.grid)
  nlex .= -cellsize(sys).*gfnl().*(8/π/σx^3/σy^3).*xg.*yg'.*(1/σx^2 .- 1/σy^2)

  @test norm(dw-nlex,sys.grid)/norm(nlex,sys.grid) < 1e-3

  end

  struct OseenPressure <: CartesianGrids.AbstractSpatialField
    σ :: Real
    x0 :: Real
    y0 :: Real
    A :: Real
  end
  function (p::OseenPressure)(x,y)
    r = sqrt((x-p.x0)^2+(y-p.y0)^2)
    p.A^2/(4π^2)*(-0.5/r^2 + 1/r^2*exp(-r^2/p.σ^2) -
                  expint(r^2/p.σ^2)/p.σ^2 -
                  0.5/r^2*exp(-2*r^2/p.σ^2) +
                  expint(2*r^2/p.σ^2)/p.σ^2)
   end


  @testset "Lamb-Oseen vortex" begin

    Re = 200 + 50rand()
    U = 0.2randn()
    U∞ = (U,0.0)

    xlim = (-2.0,2.0)
    ylim = (-2.0,2.0)

    Δx, Δt = setstepsizes(Re,gridRe=4)

    sys = NavierStokes(Re,Δx,xlim,ylim,Δt,freestream = U∞)

    @test isnothing(sys.motions)
    @test isnothing(sys.bodies)

    σ = 0.2; x0 = 0.0; y0 = 0.0; A = 1
    gauss = SpatialGaussian(σ,x0,y0,A)

    u0 = newstate(gauss,sys)

    tspan = (0.0,0.5)
    integrator = init(u0,tspan,sys)

    step!(integrator,0.5)

    press = pressure(integrator)

    oseen_vorticity(t) = SpatialGaussian(sqrt(σ^2+4*t/Re),x0+U*t,y0,A)
    #exactvorticity(t) = newstate(oseen_exact(t),sys)
    exactvorticity(t) = CartesianGrids.GeneratedField(state(u0),oseen_vorticity(t),sys.grid)

    oseen_pressure(t) = OseenPressure(sqrt(σ^2+4*t/Re),x0+U*t,y0,A)
    exactpressure(t) = CartesianGrids.GeneratedField(press,oseen_pressure(t),sys.grid)

    @test norm(vorticity(integrator)-exactvorticity(integrator.t)()) < 1e-2

    @test norm(press - exactpressure(integrator.t)()) < 5e-2


  end

  @testset "Flow past cylinder" begin

    Re = 200
    U = 1.0
    U∞ = (U,0.0)

    xlim = (-1.0,3.0)
    ylim = (-1.0,1.0)

    Δx, Δt = setstepsizes(Re,gridRe=4)

    @test Δx == 0.02
    @test Δt == min(0.5*Δx,0.5*Δx^2*Re)

    body = Circle(0.5,1.5Δx)


    sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,freestream = U∞)

    motion = RigidBodyTools.RigidBodyMotion(0.0,0.0)
    @test sys.motions[1].kin == motion.kin

    @test size(sys,1) == 208
    @test size(sys,2) == 104
    @test size(sys) == (208,104)

    @test origin(sys) == (54,52)

    u0 = newstate(sys)

    w₀ = Nodes(Dual,size(sys))
    τ = VectorData(sys.points)
    @test state(u0) == w₀
    @test constraint(u0) == τ

    du = deepcopy(u0)
    sys.f(du,u0,sys,0.0)
    @test maximum(state(du)) == minimum(state(du)) == 0.0
    @test all(x -> (x ≈ -U),constraint(du).u)
    @test all(x -> (x ≈ 0.0),constraint(du).v)

    tspan = (0.0,1.0)
    integrator = init(u0,tspan,sys)

    @test integrator.alg isa LiskaIFHERK{Direct}

    step!(integrator,0.2)

  end


end
