using Fields
import IntFactSystems: A⁻¹, r₁

struct NavierStokes{NX, NY} <: System{Unconstrained}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Tuple{Float64, Float64}

    # Discretization
    "Grid spacing"
    Δx
    "Time step (used to determine integrating factor diffusion rate)"
    Δt

    # Operators
    "Convolution operator for the integrating factor"
    E::CircularConvolution{NX, NY}
    "Convolution operator for LGF-based inverse Laplacian"
    invlap::CircularConvolution{NX, NY}

    # Scratch space

    ## Required structures for time marching
    A⁻¹g::Nodes{Dual,NX, NY}
    Ñ::Nodes{Dual,NX, NY}
    q::Nodes{Dual,NX, NY}
    w::Vector{Nodes{Dual,NX, NY}}

    ## Pre-allocated space for intermediate values
    Cs::Edges{Primal, NX, NY}
    Ww::Edges{Dual, NX, NY}
    QCs::Edges{Dual, NX, NY}

end

function NavierStokes(dims::Tuple{Int, Int}, Re, Δx, Δt, U∞ = (0.0, 0.0), nw = 2)
    NX, NY = dims

    α = Δt/(Re*Δx^2)
    qtab = [intfact(x, y, 0.5α) for x in 0:NX-1, y in 0:NY-1]
    E = CircularConvolution(qtab, FFTW.PATIENT)
    invlap = CircularConvolution(view(Fields.LGF_TABLE, 1:NX, 1:NY), FFTW.PATIENT)

    A⁻¹g = Nodes{Dual,NX, NY}()
    Ñ    = Nodes{Dual,NX, NY}()
    q = Nodes{Dual,NX, NY}()
    w = [Nodes{Dual,NX, NY}() for i in 1:nw]

    Cs   = Edges{Primal, NX, NY}()
    Ww   = Edges{Dual, NX, NY}()
    QCs  = Edges{Dual, NX, NY}()


    NavierStokes{NX, NY}(Re, U∞, Δx, Δt, E, invlap, A⁻¹g, Ñ, q, w, Cs, Ww, QCs)
end

function r₁(Ñ::Nodes{Dual,NX, NY}, w, t, sys::NavierStokes{NX, NY}) where {NX, NY}
    Cs           = sys.Cs
    Ww = Ww_QCs  = sys.Ww
    QCs          = sys.QCs
    s = Ñ

    A_mul_B!(s, sys.invlap, w)
    curl!(Cs, s)
    Fields.shift!(QCs, Cs)
    QCs.u .-= sys.U∞[1]
    QCs.v .-= sys.U∞[2]

    Fields.shift!(Ww, w)
    Ww.u[:,1] = 0; Ww.u[:,end] = 0
    Ww.v[1,:] = 0; Ww.v[end,:] = 0

    #fill!(Ñ, 0)
    Ñ[:,1] = 0; Ñ[:,end] = 0
    Ñ[1,:] = 0; Ñ[end,:] = 0
    product!(Ww_QCs, Ww, QCs)
    divergence!(Ñ, Ww_QCs)
    Δx⁻¹ = 1/sys.Δx
    scale!(Ñ, Δx⁻¹)

    Ñ
end

A⁻¹(out, u, sys::NavierStokes) = A_mul_B!(out, sys.E, u)
