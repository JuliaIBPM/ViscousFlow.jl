using PkgBenchmark

include(joinpath(Pkg.dir("Whirl2d"), "src/Whirl2d.jl"))
import Whirl2d
import Whirl2d:@get, Soln, IntFactSystems, TimeMarching, Fields

using IntFactSystems
include(joinpath(Pkg.dir("Whirl2d"), "src/systems/navier_stokes.jl"))
import TimeMarching: ifrk!, RK31

@benchgroup "Navier Stokes" begin
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
end

@benchgroup "Fields"  begin
    celldims = (500, 500)

    nodes = DualNodes(celldims)
    nodes .= rand(size(nodes))

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
