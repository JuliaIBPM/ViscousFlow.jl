mutable struct NavierStokes{NX, NY, N, isstatic}  #<: System{Unconstrained}
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

    # Body coordinate data, if present
    # if a static problem, these coordinates are in inertial coordinates
    # if a non-static problem, in their own coordinate systems
    X̃::VectorData{N}

    # Pre-stored regularization and interpolation matrices (if present)
    Hmat::Nullable{RegularizationMatrix}
    Emat::Nullable{InterpolationMatrix}


    # Scratch space

    ## Pre-allocated space for intermediate values
    Vb::VectorData{N}
    Fq::Edges{Primal, NX, NY}
    Ww::Edges{Dual, NX, NY}
    Qq::Edges{Dual, NX, NY}

    # Flags
    _isstore :: Bool

end

function NavierStokes(dims::Tuple{Int, Int}, Re, Δx, Δt;
                       U∞ = (0.0, 0.0), X̃ = VectorData{0}(),
                       isstore = false,
                       isstatic = true,
                       rk::TimeMarching.RKParams=TimeMarching.RK31)
    NX, NY = dims

    α = Δt/(Re*Δx^2)

    L = plan_laplacian((NX,NY),with_inverse=true)

    Vb = VectorData(X̃)
    Fq = Edges{Primal,NX,NY}()
    Ww   = Edges{Dual, NX, NY}()
    Qq  = Edges{Dual, NX, NY}()
    N = length(X̃)÷2

    if length(N) > 0 && isstore && isstatic
      # in this case, X̃ is assumed to be in inertial coordinates
      regop = Regularize(X̃,Δx;issymmetric=true)
      Hmat, Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())
    else
      Hmat = Nullable{RegularizationMatrix}()
      Emat = Nullable{InterpolationMatrix}()
    end

    # should be able to set up time marching operator here...

    NavierStokes{NX, NY, N, isstatic}(Re, U∞, Δx, Δt, rk, L, X̃, Hmat, Emat, Vb, Fq, Ww, Qq, isstore)
end

function Base.show(io::IO, sys::NavierStokes{NX,NY,N,isstatic}) where {NX,NY,N,isstatic}
    print(io, "Navier-Stokes system on a grid of size $NX x $NY")
end


Fields.plan_intfact(t,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
        Fields.plan_intfact(t/(sys.Re*sys.Δx^2),w)

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

function TimeMarching.r₁(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY},
                            motion::RigidBodyMotions.RigidBodyMotion) where {NX,NY}
    _,ċ,_,_,α̇,_ = motion(t)
    return TimeMarching.r₁(u[1],t,sys), [real(ċ),imag(ċ),α̇]
end

function TimeMarching.r₂(w::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,N}
    ΔV = VectorData(sys.X̃)
    ΔV.u .-= sys.U∞[1]
    ΔV.v .-= sys.U∞[2]
    return ΔV
end

function TimeMarching.r₂(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false},
                            motion::RigidBodyMotions.RigidBodyMotion) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.
  xc, yc, α = u[2]
  T = Bodies.RigidTransform((xc,yc),α)
  x, y = T(sys.X̃.u,sys.X̃.v)

  ΔV = VectorData(sys.X̃)
  _,ċ,_,_,α̇,_ = motion(t)
  for i = 1:N
      Δz = (x[i]-xc)+im*(y[i]-yc)
      ċi = ċ + im*α̇*Δz
      ΔV.u[i] = real(ċi)
      ΔV.v[i] = imag(ċi)
  end
  ΔV.u .-= sys.U∞[1]
  ΔV.v .-= sys.U∞[2]
  return ΔV, Vector{Float64}()
end


TimeMarching.B₁ᵀ(f,sys::NavierStokes{NX,NY,N,C}) where {NX,NY,N,C} = Curl()*(get(sys.Hmat)*f)
TimeMarching.B₂(w,sys::NavierStokes{NX,NY,N,C}) where {NX,NY,N,C} = -(get(sys.Emat)*(Curl()*(sys.L\w)))

TimeMarching.B₁ᵀ(f::VectorData{N},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N} = Curl()*regop(sys.Fq,f)
TimeMarching.B₂(w::Nodes{Dual,NX,NY},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N} = -(regop(sys.Vb,Curl()*(sys.L\w)))



TimeMarching.plan_constraints(w::Nodes{Dual,NX,NY},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,N} =
                    (f -> TimeMarching.B₁ᵀ(f,sys),w -> TimeMarching.B₂(w,sys))

function TimeMarching.plan_constraints(u::Tuple{Nodes{Dual,NX,NY},Vector{Float64}},t,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N}

  # for now, just assume that there is only one body. will fix later.

  xc, yc, α = u[2]
  T = Bodies.RigidTransform((xc,yc),α)
  # should be able to save some time and memory allocation here...
  x, y = T(sys.X̃.u,sys.X̃.v)
  X = VectorData(x,y)
  regop = Regularize(X,sys.Δx;issymmetric=true)
  if sys._isstore
    Hmat, Emat = RegularizationMatrix(regop,VectorData{N}(),Edges{Primal,NX,NY}())
    sys.Hmat = Hmat
    sys.Emat = Emat
    return (f->TimeMarching.B₁ᵀ(f,sys), f->zeros(Float64,size(u[2]))),
           (w->TimeMarching.B₂(w,sys), u->Vector{Float64}())
  else
    return (f->TimeMarching.B₁ᵀ(f,regop,sys), f->zeros(Float64,size(u[2]))),
           (w->TimeMarching.B₂(w,regop,sys), u->Vector{Float64}())
  end



end

# Here, need to write plan_constraints functions that compute B₁ᵀand B₂ for
# cases with and without static bodies

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
