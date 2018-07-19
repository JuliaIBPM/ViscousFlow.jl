import Whirl: Systems
using Systems


@testset "Lamb-Oseen vortex" begin

    woseen(x::Tuple{Real,Real},t;Re=1.0,x0::Tuple{Real,Real}=(0,0),t0=1) =
                          exp(-((x[1]-x0[1])^2+(x[2]-x0[2])^2)/(4(t+t0)/Re))/(1+t/t0)

    #wexact(x, y, t, Re, σ, U, Δx) = Δx*exp(-((x-U*t)^2+y^2)/(σ^2+4t/Re))/(1+4t/(Re*σ^2))

    nx = 260; ny = 130
    Re = 200 + 50rand()
    U = 1.0 + 0.2randn()
    U∞ = (U,0.0)
    Lx = 2.0
    Δx = Lx/(ny-1)
    Δt = min(0.5*Δx,0.5*Δx^2*Re)
    w₀ = Nodes(Dual,(nx,ny))

    xg,yg = coordinates(w₀,dx=Δx)
    x0 = (1,1)
    t0 = 1
    wexact(t) = [woseen((x,y),t;Re=Re,x0=x0.+U∞.*t,t0=t0) for x in xg, y in yg]

    sys = Systems.NavierStokes((nx,ny),Re,Δx,Δt,U∞ = U∞)

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
    @test sum(w) ≈ sum(wexact(t))

end
