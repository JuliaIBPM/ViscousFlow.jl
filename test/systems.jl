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
    body = Bodies.Circle(0.5,n)

    xlim = (-1.0,3.0)
    ylim = (-1.0,1.0)

    X = VectorData(body.x,body.y)

    Δx = 0.02
    Δt = min(0.5*Δx,0.5*Δx^2*Re)

    sys = Systems.NavierStokes(Re,Δx,xlim,ylim,Δt,U∞ = U∞, X̃ = X, isstore = true)

    @test size(sys,1) == 208
    @test size(sys,2) == 104
    @test size(sys) == (208,104)

    @test Systems.origin(sys) == (54,52)

    wf = Systems.PointForce(Nodes(Dual,size(sys)),(1.5,0.0),10.0,1.5,1.0,sys)
    @test sum(wf(1.5)) ≈ 10

    qf = Systems.PointForce(Edges(Primal,size(sys)),(1.5,0.0),(10.0,-10.0),1.5,1.0,sys)
    @test sum(qf(1.5).u) ≈ 10
    @test sum(qf(1.5).v) ≈ -10



  end

  @testset "Histories" begin

    a = 1.0
    b = 0.0
    c = 5.0

    h = History(Nodes(Dual,(5,5)))
    d = Nodes(Dual,(5,5))
    fill!(d,a)
    push!(h,deepcopy(d))
    fill!(d,b)
    push!(h,deepcopy(d))
    fill!(d,c)
    push!(h,deepcopy(d))

    @test length(h.r) == length(h)

    dh = diff(h)

    @test typeof(dh) <: History

    @test dh[2][1,1] == c-b

    # test periodic history
    hp = History(h.vec,htype=PeriodicHistory)
    @test hp[4] == hp[1]

    # test the differencing of periodic history
    dhp = diff(hp)
    @test typeof(dhp) <: History
    @test dhp[1][1,1] == b-a
    @test dhp[5] == dhp[2]

    # test staggered indexing
    hp2 = History(hp[1:2:5],htype=PeriodicHistory)
    @test hp2[2][1,1] == c

    # test circular shifting
    hpshift = circshift(hp,1)
    @test hpshift[1][1,1] == c

    # testing arithmetic
    h3 = 3*h
    @test h3[3][1,1] == 3*c

    h4 = h+h3
    @test h4[1][1,1] == 4*a

    h5 = -h
    @test h5[3][1,1] == -c

    h6 = h/2
    @test h6[1][1,1] == a/2

    # test setting ghosts on regular history
    h_pre = deepcopy(h)
    h_post = deepcopy(h)
    set_first_ghost!(h,h_pre)
    @test h[1][1,1] == c

    set_last_ghost!(h,h_post)
    @test h[end][1,1] == a

    @test length(h) == length(h.r)+2
    @test h.r.start == 2
    @test h.r.stop == 4



  end

end
