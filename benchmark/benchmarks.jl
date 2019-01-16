using PkgBenchmark
using BenchmarkHelper

import Whirl: Systems, TimeMarching, Fields
using Fields
using Systems
using TimeMarching
#using SaddlePointSystems

#include(joinpath(Pkg.dir("Whirl"), "src/Whirl.jl"))
#import Whirl
#import Whirl:@get, Soln, IntFactSystems, TimeMarching, Fields

#using IntFactSystems
#include(joinpath(Pkg.dir("Whirl"), "src/systems/navier_stokes.jl"))
#import TimeMarching: ifrk!, RK31

@benchgroup "Navier Stokes" begin
  woseen(x::Tuple{Real,Real},t;Re=1.0,x0::Tuple{Real,Real}=(0,0),t0=1) =
                      exp(-((x[1]-x0[1])^2+(x[2]-x0[2])^2)/(4(t+t0)/Re))/(1+t/t0)


  nx = 260; ny = 130
  Re = 200
  U = 1.0
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
  @bench "Single-step NS solver" $ifrk($t,$w)

#=
    wexact(x, y, t) = exp(-((x-t)^2+y^2)/(0.2^2+4t/200))/(1+4t/(200*0.2^2))

    Re = 200
    Δx = 0.02
    x = linspace(-1-0.5Δx, 3+0.5Δx, 202)
    y = linspace(-1-0.5Δx, 1+0.5Δx, 102)
    Δt = min(0.5*Δx,0.5*Δx^2*Re)
    dims = length.((x,y))

    ns = NavierStokes(dims, Re, Δx, Δt, (1.0, 0.0))
    s = Soln(0.0, Fields.DualNodes(dims));
    s.u .= [wexact(xᵢ, yᵢ, 0)*Δx for xᵢ in x, yᵢ in y]
    s₊ = Soln(0.0, Fields.DualNodes(dims));


    @bench "Single-Step NS Solver" ifrk!($s₊, $s, $Δt, $RK31, $ns)
    =#
end

#=
@benchgroup "Fields"  begin
    celldims = (500, 500)

    nodes = DualNodes(celldims)
    nodes .= rand(Float64,size(nodes))

    edges = Edges(Primal, nodes)
    edges.u .= reshape(1:length(edges.u), size(edges.u))
    edges.v .= reshape(1:length(edges.v), size(edges.v))

    dual = Edges(Dual, nodes)

    L = Laplacian(celldims, with_inverse = true, fftw_flags = FFTW.PATIENT)

    @bench "Shift Edges to Edges" Fields.shift!($dual, $edges)
    @bench "Shift Nodes to Edges" Fields.shift!($dual, $nodes)

    @bench "Node to Edge Curl" curl!($edges, $nodes)
    @bench "Edge to Node Curl" curl!($nodes, $edges)

    out = DualNodes(celldims)
    @bench "Laplacian" A_mul_B!($out, $L, $nodes)
    @bench "Inverse Laplacian" A_ldiv_B!($out, $L, $nodes)
end
=#
