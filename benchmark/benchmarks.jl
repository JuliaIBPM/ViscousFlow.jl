using PkgBenchmark

include(joinpath(Pkg.dir("Whirl2d"), "src/Whirl2d.jl"))
import Whirl2d
import Whirl2d:@get
@get Whirl2d (Systems, Grids, DualPatch, Bodies, TimeMarching, NavierStokes);

@benchgroup "Navier Stokes" begin
    xmin = [-1.0,-1.0]
    xmax = [3.0,1.0]
    dom = Systems.DualDomain(xmin,xmax)

    Δx = 0.02
    dom = Systems.add_grid(dom,Δx)

    Re = 200
    physparams = NavierStokes.set_freestream([1.0,0.0])
    NavierStokes.set_Re!(physparams,Re)
    Δt = min(0.5*Δx,0.5*Δx^2*Re)

    α = Δt/(Re*Δx^2)
    tparams = TimeMarching.TimeParams(Δt,TimeMarching.RK31())

    params = (physparams,α)
    ops = NavierStokes.set_operators!(dom,params);

    s = NavierStokes.Soln(dom)

    x = Grids.xcell(dom.grid)
    y = Grids.ycell(dom.grid)
    σ = 0.2
    wexact(t) = [exp(-((x[i]-t)^2+y[j]^2)/(σ^2+4t/Re))/(1+4t/(Re*σ^2)) for i = 1:length(x), j = 1:length(y)];
    s.u[dom.grid.cellint[1],dom.grid.cellint[2]] = dom.grid.Δx*wexact(0.0);


    @bench "Single-Step NS Solver" TimeMarching.ifrk!($s,$tparams,$ops)
end
