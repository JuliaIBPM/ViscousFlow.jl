import Base: size

import ConstrainedSystems: r₁, r₂, B₁ᵀ, B₂, plan_constraints
import CartesianGrids: cellsize, origin

"""
$(TYPEDEF)

A system type that utilizes a grid of `NX` x `NY` dual cells and `N` Lagrange forcing
points to solve the discrete Navier-Stokes equations in vorticity form. The
parameter `static_points` specifies whether the forcing points remain static in the
grid.

# Constructors:

`NavierStokes(Re,Δx,xlimits,ylimits,Δt
              [,freestream = (0.0, 0.0)][,points = VectorData(0)]
              [,store_operators=true][,static_points=true]
              [,rk=ConstrainedSystems.RK31]
              [,ddftype=CartesianGrids.Yang3])`specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. Note that `store_operators` set to `true`
              stores matrix versions of the operators. This makes the method
              faster, at the cost of storage.

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList)` passes the body information
              directly. The other keywords can be supplied, although `points`
              would be ignored.

"""
mutable struct NavierStokes{NX, NY, N, MT<:MotionType}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Tuple{Float64, Float64}
    "Body motions"
    motions::Union{Vector{RigidBodyMotion},Nothing}

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
    points::VectorData{N,Float64}

    # Pre-stored regularization and interpolation matrices (if present)
    Hmat::Union{RegularizationMatrix,Nothing}
    Emat::Union{InterpolationMatrix,Nothing}

    # Conditioner matrices
    Cmat::Union{AbstractMatrix,Nothing}

    # Scratch space

    ## Pre-allocated space for intermediate values
    Vb::VectorData{N,Float64}
    Fq::Edges{Primal, NX, NY, Float64}
    Ww::Edges{Dual, NX, NY,Float64}
    Qq::Edges{Dual, NX, NY,Float64}

    # Flags
    _isstored :: Bool

end

function NavierStokes(Re::Real, Δx::Real, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt::Real;
                       freestream = (0.0, 0.0), points = VectorData(0),
                       motions::Union{Vector{RigidBodyMotion},Nothing} = nothing,
                       store_operators = true,
                       static_points = true,
                       rk::ConstrainedSystems.RKParams=ConstrainedSystems.RK31,
                       ddftype=CartesianGrids.Yang3)

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    α = Δt/(Re*Δx^2)

    L = plan_laplacian((NX,NY),with_inverse=true)

    Vb = VectorData(points)
    Fq = Edges{Primal,NX,NY,Float64}()
    Ww = Edges{Dual, NX, NY,Float64}()
    Qq = Edges{Dual, NX, NY,Float64}()
    N = length(points)÷2

    Hmat = nothing
    Emat = nothing
    Cmat = nothing

    if N > 0 && store_operators && static_points
      # in this case, points are assumed to be in inertial coordinates

      regop = Regularize(points,Δx;I0=CartesianGrids.origin(g),issymmetric=true,ddftype=ddftype)
      Hmat, Emat = RegularizationMatrix(regop,Vb,Fq)
      regopfilt = Regularize(points,Δx;I0=CartesianGrids.origin(g),filter=true,weights=Δx^2,ddftype=ddftype)
      Ẽmat = InterpolationMatrix(regopfilt,Fq,Vb)
      Cmat = sparse(Ẽmat*Hmat)
    end

    NavierStokes{NX, NY, N, _motiontype(static_points)}(Re, freestream, motions,
                                        g, Δt, rk,
                                        L,
                                        points, Hmat, Emat,Cmat,
                                        Vb, Fq, Ww, Qq, store_operators)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::Union{Body,BodyList};kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt;points=VectorData(collect(bodies)),kwargs...)

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,motion::RigidBodyMotion;kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt;motions=[motion],points=VectorData(collect(body)),kwargs...)

function NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList,motions::Vector{RigidBodyMotion};kwargs...)
    length(bodies) == length(motions) || error("Inconsistent lengths of bodies and motions lists")
    NavierStokes(Re,Δx,xlim,ylim,Δt;motions=motions,points=VectorData(collect(bodies)),kwargs...)
end


function Base.show(io::IO, sys::NavierStokes{NX,NY,N,MT}) where {NX,NY,N,MT}
    mtype = (MT == StaticBodies) ? "static" : "moving"
    print(io, "Navier-Stokes system on a grid of size $NX x $NY and $N $mtype immersed points")
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
_motiontype(isstatic::Bool) = isstatic ? StaticBodies : MovingBodies

include("navierstokes/fields.jl")
include("navierstokes/pulses.jl")
include("navierstokes/basicoperators.jl")
include("navierstokes/rigidbodyoperators.jl")
include("navierstokes/movingbodyoperators.jl")
