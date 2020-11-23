
using LinearAlgebra

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

    #sys = NavierStokes((nx,ny),Re,Δx,Δt,U∞ = U∞)
    sys = NavierStokes(Re,Δx,(0.0,3.0),(0.0,2.0),Δt,freestream = U∞)

    w₀ = Nodes(Dual,size(sys))
    xg,yg = coordinates(w₀,dx=Δx)
    x0 = (1,1)
    t0 = 1
    wexact(t) = [Δx*woseen((x,y),t;Re=Re,x0=x0.+U∞.*t,t0=t0) for x in xg, y in yg]

    ifrk = IFRK(w₀,sys.Δt,
                (t,w) -> plan_intfact(t,w,sys),
                (w,t) -> r₁(w,t,sys) ,rk=ConstrainedSystems.RK31)

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
    body = Circle(0.5,n)

    xlim = (-1.0,3.0)
    ylim = (-1.0,1.0)

    X = VectorData(body.x,body.y)

    Δx = 0.02
    Δt = min(0.5*Δx,0.5*Δx^2*Re)
    Δx2, Δt2 = setstepsizes(Re,gridRe=4)
    @test Δx == Δx2 && Δt == Δt2

    sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,freestream = U∞)

    @test size(sys,1) == 208
    @test size(sys,2) == 104
    @test size(sys) == (208,104)

    @test origin(sys) == (54,52)

    #=
    wf = PointForce(Nodes(Dual,size(sys)),(1.5,0.0),10.0,1.5,1.0,sys)
    @test isapprox(sum(wf(1.5)),10,atol=1e-2)

    qf = PointForce(Edges(Primal,size(sys)),(1.5,0.0),(10.0,-10.0),1.5,1.0,sys)
    @test isapprox(sum(qf(1.5).u),10,atol=1e-2)
    @test isapprox(sum(qf(1.5).v),10,atol=1e-2)
    =#

    w₀ = Nodes(Dual,size(sys))
    f = VectorData(X)

    plan_intfact(t,u) = CartesianGrids.plan_intfact(t,u,sys)
    plan_constraints(u,t) = ConstrainedSystems.plan_constraints(u,t,sys)
    r₁(u,t) = ConstrainedSystems.r₁(u,t,sys)
    r₂(u,t) = ConstrainedSystems.r₂(u,t,sys)


    solver = IFHERK(w₀,f,timestep(sys),plan_intfact,plan_constraints,(r₁,r₂),rk=ConstrainedSystems.RK31)

    t = 0.0
    u = zero(w₀)

    tsim = 0.2
    @test timerange(tsim,sys) == Δt:Δt:tsim

    for ti in timerange(tsim,sys)
      t, u, f = solver(t,u)
    end

  end



end
