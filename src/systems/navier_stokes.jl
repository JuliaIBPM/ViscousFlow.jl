struct NavierStokes{NX, NY}  #<: System{Unconstrained}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Tuple{Float64, Float64}

    # Discretization
    "Grid spacing"
    Δx::Float64
    "Time step (used to determine integrating factor diffusion rate)"
    Δt::Float64
    "Runge-Kutta method"
    rk::TimeMarching.RKParams

    # Operators
    "Laplacian operator"
    L::Fields.Laplacian{NX,NY}

    # Scratch space

    ## Pre-allocated space for intermediate values
    Ww::Edges{Dual, NX, NY}
    Qq::Edges{Dual, NX, NY}

end

function NavierStokes(dims::Tuple{Int, Int}, Re, Δx, Δt;
                       U∞ = (0.0, 0.0),
                       rk::TimeMarching.RKParams=TimeMarching.RK31)
    NX, NY = dims

    α = Δt/(Re*Δx^2)

    L = plan_laplacian((NX,NY),with_inverse=true)

    Ww   = Edges{Dual, NX, NY}()
    Qq  = Edges{Dual, NX, NY}()

    NavierStokes{NX, NY}(Re, U∞, Δx, Δt, rk, L, Ww, Qq)
end

function Base.show(io::IO, sys::NavierStokes{NX,NY}) where {NX,NY}
    print(io, "Navier-Stokes system on a grid of size $NX x $NY")
end

function TimeMarching.r₁(w,t,sys::NavierStokes{NX,NY}) where {NX,NY}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/sys.Δx

  shift!(Qq,curl(L\w)) # -velocity, on dual edges
  Qq.u .-= sys.U∞[1]
  Qq.v .-= sys.U∞[2]

  return scale!(divergence(Qq∘shift!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

Fields.plan_intfact(t,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
        Fields.plan_intfact(t/(sys.Re*sys.Δx^2),w)


#=

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

=#
