import Whirl: IntFactSystems, TimeMarching

using IntFactSystems
include(joinpath(Pkg.dir("Whirl"), "src/systems/navier_stokes.jl"))
import TimeMarching: ifrk!, RK31
import Whirl: Soln

using Fields

@testset "Gaussian Vorticity Patch" begin
    wexact(x, y, t, Re, σ, U, Δx) = Δx*exp(-((x-U*t)^2+y^2)/(σ^2+4t/Re))/(1+4t/(Re*σ^2))

    Re = 200 + 50rand()
    Δx = 0.02
    x = linspace(-1-0.5Δx, 3+0.5Δx, 202)
    y = linspace(-1-0.5Δx, 1+0.5Δx, 102)
    Δt = min(0.5*Δx,0.5*Δx^2*Re)
    dims = length.((x,y))

    U = 1.0 + 0.2randn()
    σ = 0.2 + 0.05randn()

    ns = NavierStokes(dims, Re, Δx, Δt, (U, 0.0))
    s = Soln(0.0, Fields.DualNodes(dims));
    s.u .= [wexact(xᵢ, yᵢ, 0, Re, σ, U, Δx) for xᵢ in x, yᵢ in y]
    s₊ = Soln(0.0, Fields.DualNodes(dims));

    for i in 1:10
        ifrk!(s₊, s, Δt, RK31, ns)
        s₊, s = s, s₊
    end

    exact_soln = [wexact(xᵢ, yᵢ, s.t, Re, σ, U, Δx) for xᵢ in x, yᵢ in y]
    @test norm(s.u - exact_soln) < 1e-3
    @test sum(s.u) ≈ sum(exact_soln)
end
