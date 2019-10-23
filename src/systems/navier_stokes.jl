import Base: size


"""
$(TYPEDEF)

A system type that utilizes a grid of `NX` x `NY` dual cells and `N` Lagrange forcing
points to solve the discrete Navier-Stokes equations in vorticity form. The
parameter `isstatic` specifies whether the forcing points remain static in the
grid.

# Fields
- `Re`: Reynolds number
- `U∞`: Tuple of components of free-stream velocity
- `Δx`: Size of each side of a grid cell
- `I0`: Tuple of indices of the primal node corresponding to physical origin
- `Δt`: Time step
- `rk`: Runge-Kutta coefficients
- `L`: Pre-planned discrete Laplacian operator and inverse
- `X̃`: Lagrange point coordinate data (if present), expressed in inertial coordinates
        (if static) or in body-fixed coordinates (if moving)
- `Hmat`: Pre-computed regularization matrix (if present)
- `Emat`: Pre-computed interpolation matrix (if present)
- `Vb`: Buffer space for vector data on Lagrange points
- `Fq`: Buffer space for primal cell edge data
- `Ww`: Buffer space for dual cell edge data
- `Qq`: More buffer space for dual cell edge data
- `_isstore`: flag to specify whether to store regularization/interpolation matrices

# Constructors:

`NavierStokes(Re,Δx,xlimits,ylimits,Δt
              [,U∞ = (0.0, 0.0)][,X̃ = VectorData{0}()]
              [,isstore=false][,isstatic=true][,isfilter=false]
              [,rk=TimeMarching.RK31]
              [,ddftype=Fields.Yang3])` specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. Note that `isstore` set to `true`
              would store matrix versions of the operators. This makes the method
              faster, at the cost of storage. If `isfilter` is set to true, then
              the regularization relies on a filtered version.

"""
mutable struct NavierStokes{NX, NY, N, isstatic}  #<: System{Unconstrained}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Tuple{Float64, Float64}

    # Discretization
    "Grid metadata"
    grid::Fields.PhysicalGrid{2}
    #"Grid spacing"
    #Δx::Float64
    #"Indices of the primal node corresponding to the physical origin"
    #I0::Tuple{Int,Int}
    "Time step"
    Δt::Float64
    "Runge-Kutta method"
    rk::TimeMarching.RKParams

    # Operators
    "Laplacian operator"
    L::Fields.Laplacian{NX,NY}

    # Body coordinate data, if present
    # if a static problem, these coordinates are in inertial coordinates
    # if a non-static problem, in their own coordinate systems
    X̃::VectorData{N,Float64}

    # Pre-stored regularization and interpolation matrices (if present)
    Hmat::Union{RegularizationMatrix,Nothing}
    Emat::Union{InterpolationMatrix,Nothing}
    Hmat_grad::Union{RegularizationMatrix,Nothing}
    Emat_grad::Union{InterpolationMatrix,Nothing}

    # Conditioner matrices
    Cmat::Union{AbstractMatrix,Nothing}
    Cmat_grad::Union{AbstractMatrix,Nothing}

    # Scratch space

    ## Pre-allocated space for intermediate values
    Vb::VectorData{N,Float64}
    Fq::Edges{Primal, NX, NY, Float64}
    Ww::Edges{Dual, NX, NY,Float64}
    Qq::Edges{Dual, NX, NY,Float64}

    # Flags
    _isstore :: Bool

end

function NavierStokes(Re, Δx, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt;
                       U∞ = (0.0, 0.0), X̃ = VectorData{0,Float64}(),
                       isstore = false,
                       isstatic = true,
                       isasymptotic = false,
                       isfilter = false,
                       rk::TimeMarching.RKParams=TimeMarching.RK31,
                       ddftype=Fields.Yang3)

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    α = Δt/(Re*Δx^2)

    L = plan_laplacian((NX,NY),with_inverse=true)

    Vb = VectorData(X̃)
    Fq = Edges{Primal,NX,NY,Float64}()
    Ww = Edges{Dual, NX, NY,Float64}()
    Qq = Edges{Dual, NX, NY,Float64}()
    N = length(X̃)÷2

    Hmat = nothing
    Emat = nothing
    Cmat = nothing

    Hmat_grad = nothing
    Emat_grad = nothing
    Cmat_grad = nothing

    if length(N) > 0 && isstore && isstatic
      # in this case, X̃ is assumed to be in inertial coordinates

      regop = Regularize(X̃,Δx;I0=Fields.origin(g),issymmetric=true,ddftype=ddftype)
      Hmat, Emat = RegularizationMatrix(regop,Vb,Fq)
      if isfilter
        regopfilt = Regularize(X̃,Δx;I0=Fields.origin(g),filter=true,weights=Δx^2,ddftype=ddftype)
        Ẽmat = InterpolationMatrix(regopfilt,Fq,Vb)
        Cmat = sparse(Ẽmat*Hmat)
      end
      if isasymptotic
        Hmat_grad, Emat_grad = RegularizationMatrix(regop,TensorData{N}(),grad(Fq))
        if isfilter
          Ẽmat_grad = InterpolationMatrix(regopfilt,grad(Fq),TensorData{N}())
          Cmat_grad = sparse(Ẽmat_grad*Hmat_grad)
        end
      end
    end

    # should be able to set up time marching operator here...

    #NavierStokes{NX, NY, N, isstatic}(Re, U∞, Δx, I0, Δt, rk, L, X̃, Hmat, Emat, Vb, Fq, Ww, Qq, isstore)
    NavierStokes{NX, NY, N, isstatic}(Re, U∞, g, Δt, rk, L, X̃, Hmat, Emat,
                                        Hmat_grad, Emat_grad,
                                        Cmat,Cmat_grad,
                                        Vb, Fq, Ww, Qq, isstore)
end

function Base.show(io::IO, sys::NavierStokes{NX,NY,N,isstatic}) where {NX,NY,N,isstatic}
    print(io, "Navier-Stokes system on a grid of size $NX x $NY")
end

# some convenience functions
"""
    size(sys::NavierStokes,d::Int) -> Int

Return the number of indices of the grid used by `sys` along dimension `d`.
"""
size(sys::NavierStokes{NX,NY},d::Int) where {NX,NY} = d == 1 ? NX : NY

"""
    size(sys::NavierStokes) -> Tuple{Int,Int}

Return a tuple of the number of indices of the grid used by `sys`
"""
size(sys::NavierStokes{NX,NY}) where {NX,NY} = (size(sys,1),size(sys,2))

"""
    cellsize(sys::NavierStokes) -> Float64

Return the grid cell size of system `sys`
"""
Fields.cellsize(sys::NavierStokes) = cellsize(sys.grid)

"""
    origin(sys::NavierStokes) -> Tuple{Int,Int}

Return a tuple of the indices of the primal node that corresponds to the
physical origin of the coordinate system used by `sys`. Note that these
indices need not lie inside the range of indices occupied by the grid.
For example, if the range of physical coordinates occupied by the grid
is (1.0,3.0) x (2.0,4.0), then the origin is not inside the grid.
"""
Fields.origin(sys::Systems.NavierStokes) = origin(sys.grid)

# Basic operators for any Navier-Stokes system

# Integrating factor -- rescale the time-step size
Fields.plan_intfact(Δt,w,sys::NavierStokes{NX,NY}) where {NX,NY} =
        Fields.plan_intfact(Δt/(sys.Re*cellsize(sys)^2),w)

# RHS of Navier-Stokes (non-linear convective term)
function TimeMarching.r₁(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY}) where {NX,NY,T}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/cellsize(sys)

  interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
  Qq.u .-= sys.U∞[1]
  Qq.v .-= sys.U∞[2]

  return rmul!(divergence(Qq∘interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# RHS of Navier-Stokes (non-linear convective term)
function TimeMarching.r₁(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY,T}

  Ww = sys.Ww
  Qq = sys.Qq
  L = sys.L
  Δx⁻¹ = 1/cellsize(sys)

  interpolate!(Qq,curl(L\w)) # -velocity, on dual edges
  _,ċ,_,_,_,_ = U∞(t)
  Qq.u .-= real(ċ)
  Qq.v .-= imag(ċ)

  return rmul!(divergence(Qq∘interpolate!(Ww,w)),Δx⁻¹) # -∇⋅(wu)

end

# Operators for a system with a body

# RHS of a stationary body with no surface velocity
function TimeMarching.r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,N,T}
    ΔV = VectorData(sys.X̃)
    ΔV.u .-= sys.U∞[1]
    ΔV.v .-= sys.U∞[2]
    return ΔV
end

function TimeMarching.r₂(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true},U∞::RigidBodyMotions.RigidBodyMotion) where {NX,NY,N,T}
    ΔV = VectorData(sys.X̃)
    _,ċ,_,_,_,_ = U∞(t)
    ΔV.u .-= real(ċ)
    ΔV.v .-= imag(ċ)
    return ΔV
end

# Constraint operators, using stored regularization and interpolation operators
# B₁ᵀ = CᵀEᵀ, B₂ = -ECL⁻¹
TimeMarching.B₁ᵀ(f::VectorData{N},sys::NavierStokes{NX,NY,N,C}) where {NX,NY,N,C} = Curl()*(sys.Hmat*f)
TimeMarching.B₂(w::Nodes{Dual,NX,NY,T},sys::NavierStokes{NX,NY,N,C}) where {NX,NY,T,N,C} = -(sys.Emat*(Curl()*(sys.L\w)))

# Constraint operators, using non-stored regularization and interpolation operators
TimeMarching.B₁ᵀ(f::VectorData{N},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,N} = Curl()*regop(sys.Fq,f)
TimeMarching.B₂(w::Nodes{Dual,NX,NY,T},regop::Regularize,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,T,N} = -(regop(sys.Vb,Curl()*(sys.L\w)))

# Constraint operator constructors
# Constructor using stored operators
TimeMarching.plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,true}) where {NX,NY,T,N} =
                    (f -> TimeMarching.B₁ᵀ(f,sys),w -> TimeMarching.B₂(w,sys))

# Constructor using non-stored operators
function TimeMarching.plan_constraints(w::Nodes{Dual,NX,NY,T},t,sys::NavierStokes{NX,NY,N,false}) where {NX,NY,T,N}
  regop = Regularize(sys.X̃,Fields.cellsize(sys);I0=Fields.origin(sys),issymmetric=true)

  return f -> TimeMarching.B₁ᵀ(f,regop,sys),w -> TimeMarching.B₂(w,regop,sys)
end

# compute physical values of fields
vorticity(w::Nodes{Dual},sys) = w/cellsize(sys)
velocity(w::Nodes{Dual},sys) = -curl(sys.L\w)
streamfunction(w::Nodes{Dual},sys) = -cellsize(sys)*(sys.L\w)
nl(w::Nodes{Dual},sys) = TimeMarching.r₁(w,0.0,sys)/cellsize(sys)



include("navierstokes/systemutils.jl")

include("navierstokes/movingbody.jl")
