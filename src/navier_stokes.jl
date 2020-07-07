import Base: size

import ConstrainedSystems: r₁, r₂, B₁ᵀ, B₂, plan_constraints
import CartesianGrids: cellsize, origin

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
              [,rk=ConstrainedSystems.RK31]
              [,ddftype=CartesianGrids.Yang3])` specifies the Reynolds number `Re`, the grid
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
    grid::CartesianGrids.PhysicalGrid{2}
    #"Grid spacing"
    #Δx::Float64
    #"Indices of the primal node corresponding to the physical origin"
    #I0::Tuple{Int,Int}
    "Time step"
    Δt::Float64
    "Runge-Kutta method"
    rk::ConstrainedSystems.RKParams

    # Operators
    "Laplacian operator"
    L::CartesianGrids.Laplacian{NX,NY}

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
                       U∞ = (0.0, 0.0), X̃ = VectorData(0),
                       isstore = false,
                       isstatic = true,
                       isasymptotic = false,
                       isfilter = false,
                       rk::ConstrainedSystems.RKParams=ConstrainedSystems.RK31,
                       ddftype=CartesianGrids.Yang3)

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

      regop = Regularize(X̃,Δx;I0=CartesianGrids.origin(g),issymmetric=true,ddftype=ddftype)
      Hmat, Emat = RegularizationMatrix(regop,Vb,Fq)
      if isfilter
        regopfilt = Regularize(X̃,Δx;I0=CartesianGrids.origin(g),filter=true,weights=Δx^2,ddftype=ddftype)
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

"""
    setstepsizes(Re[,gridRe=2][,cfl=0.5][,fourier=0.5]) -> Float64, Float64

Set the grid cell spacing and time step size based on the Reynolds number `Re`,
the grid Reynolds number `gridRe`, cfl number `cfl`, and grid Fourier number `fourier`.
The last three parameters all have default values.

# Example

Here is an example of setting parameters based on Reynolds number 100 (with
  default choices for grid Reynolds number, CFL number, and Fourier number):
```jldoctest
julia> Δx, Δt = setstepsizes(100)
(0.02, 0.01)
```
"""
function setstepsizes(Re::Real; gridRe = 2.0, cfl = 0.5, fourier = 0.5)
    Δx = gridRe/Re
    Δt = min(fourier*Δx,cfl*Δx^2*Re)
    return Δx, Δt
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
cellsize(sys::NavierStokes) = cellsize(sys.grid)

"""
    timestep(sys::NavierStokes) -> Float64

Return the time step size of system `sys`
"""
timestep(sys::NavierStokes) = sys.Δt

"""
    origin(sys::NavierStokes) -> Tuple{Int,Int}

Return a tuple of the indices of the primal node that corresponds to the
physical origin of the coordinate system used by `sys`. Note that these
indices need not lie inside the range of indices occupied by the grid.
For example, if the range of physical coordinates occupied by the grid
is (1.0,3.0) x (2.0,4.0), then the origin is not inside the grid.
"""
origin(sys::NavierStokes) = origin(sys.grid)


"""
    timerange(tf,sys::NavierStokes)

Create a range of times, starting at the t = Δt (the time step of `sys`),
and ending at t = `tf`.
"""
timerange(tf,sys) = timestep(sys):timestep(sys):tf


# Other functions
_hasfilter(sys::NavierStokes) = !(sys.Cmat == nothing)


include("navierstokes/fields.jl")
include("navierstokes/pulses.jl")
include("navierstokes/basicoperators.jl")
include("navierstokes/rigidbodyoperators.jl")
include("navierstokes/movingbodyoperators.jl")
