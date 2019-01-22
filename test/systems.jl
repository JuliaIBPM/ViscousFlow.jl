import ViscousFlow: Systems

using Compat
using Compat.LinearAlgebra

@testset "Navier-Stokes" begin

  @testset "Lamb-Oseen vortex" begin

    woseen(x::Tuple{Real,Real},t;Re=1.0,x0::Tuple{Real,Real}=(0,0),t0=1) =
                          exp(-((x[1]-x0[1])^2+(x[2]-x0[2])^2)/(4(t+t0)/Re))/(1+t/t0)

    #wexact(x, y, t, Re, σ, U, Δx) = Δx*exp(-((x-U*t)^2+y^2)/(σ^2+4t/Re))/(1+4t/(Re*σ^2))

    #nx = 260; ny = 130
    Re = 200 + 50rand()
    U = 1.0 + 0.2randn()
    U∞ = (U,0.0)
    #Lx = 2.0
    Δx = 0.015  #Lx/(ny-1)
    Δt = min(0.5*Δx,0.5*Δx^2*Re)

    #sys = Systems.NavierStokes((nx,ny),Re,Δx,Δt,U∞ = U∞)
    sys = Systems.NavierStokes(Re,Δx,(0.0,3.0),(0.0,2.0),Δt,U∞ = U∞)

    w₀ = Nodes(Dual,size(sys))
    xg,yg = coordinates(w₀,dx=Δx)
    x0 = (1,1)
    t0 = 1
    wexact(t) = [Δx*woseen((x,y),t;Re=Re,x0=x0.+U∞.*t,t0=t0) for x in xg, y in yg]

    ifrk = IFRK(w₀,sys.Δt,
                (t,w) -> plan_intfact(t,w,sys),
                (w,t) -> r₁(w,t,sys) ,rk=TimeMarching.RK31)

    t = 0.0
    w₀ .= wexact(t)

    w = deepcopy(w₀)
    tf = 1.0
    T = 0:Δt:tf

    for ti in T
      t, w = ifrk(t,w)
    end

    @test norm(w-wexact(t),Inf) < 1e-1
    @test abs(sum(w)-sum(wexact(t))) < 1e-3

  end

  @testset "Flow past cylinder" begin

    Re = 200
    U = 1.0
    U∞ = (U,0.0)

    n = 100
    body = Bodies.Ellipse(0.5,n)

    xlim = (-1.0,3.0)
    ylim = (-1.0,1.0)

    X = VectorData(body.x,body.y)

    Δx = 0.02
    Δt = min(0.5*Δx,0.5*Δx^2*Re)

    sys = Systems.NavierStokes(Re,Δx,xlim,ylim,Δt,U∞ = U∞, X̃ = X, isstore = true)

    @test size(sys,1) == 202
    @test size(sys,2) == 102
    @test size(sys) == (202,102)

    @test Systems.origin(sys) == (51,51)

    wf = Systems.PointForce(Nodes(Dual,size(sys)),(1.5,0.0),10.0,1.5,1.0,sys)
    @test maximum(wf(1.5))==2.5

    qf = Systems.PointForce(Edges(Primal,size(sys)),(1.5,0.0),(10.0,-10.0),1.5,1.0,sys)
    @test maximum(qf(1.5).u)==3.333333333333333
    @test minimum(qf(1.5).v)==-3.333333333333333



  end

end
