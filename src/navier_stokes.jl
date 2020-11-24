import Base: size

const NDIM = 2

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
              [,freestream = (0.0, 0.0)]
              [,store_operators=true][,static_points=true]
              [,rk=ConstrainedSystems.RK31]
              [,ddftype=CartesianGrids.Yang3])`specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. Note that `store_operators` set to `true`
              stores matrix versions of the operators. This makes the method
              faster, at the cost of storage.

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList)` passes the body
              information. The other keywords can also be supplied. This
              sets the motions of the body/ies to be stationary.

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList,
              motions::RigidBodyMotion/RigidMotionList)`
              passes the body and associated motion information.
              The other keywords can also be supplied. The list of
              motions must be the same length as the list of bodies.

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList,
              motions::RigidBodyMotion/RigidMotionList)`
              passes the body and associated motion information.
              The other keywords can also be supplied. The list of
              motions must be the same length as the list of bodies.`

"""
mutable struct NavierStokes{NX, NY, N, MT<:PointMotionType, FS<:FreestreamType}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::Union{Tuple{Float64, Float64},RigidBodyMotion}
    "Bodies"
    bodies::Union{BodyList,Nothing}
    "Body motions"
    motions::Union{RigidMotionList,Nothing}

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
    L::CartesianGrids.Laplacian
    Lc::CartesianGrids.Laplacian

    # Layers
    df::Union{DoubleLayer,Nothing}
    sc::Union{SingleLayer,Nothing}
    sn::Union{SingleLayer,Nothing}

    # Body coordinate data, if present
    # if a static problem, these coordinates are in inertial coordinates
    # if a non-static problem, in their own coordinate systems
    points::VectorData{N,Float64}

    # Pre-stored regularization and interpolation matrices (if present)
    Rf::Union{RegularizationMatrix,Nothing} # faces (edges)
    Ef::Union{InterpolationMatrix,Nothing}
    Rc::Union{RegularizationMatrix,Nothing} # cell centers
    Ec::Union{InterpolationMatrix,Nothing}
    Rn::Union{RegularizationMatrix,Nothing} # cell nodes
    En::Union{InterpolationMatrix,Nothing}

    # Conditioner matrices
    Cf::Union{AbstractMatrix,Nothing}

    # Scratch space

    ## Pre-allocated space for intermediate values
    Vb::VectorData{N,Float64}
    Sb::ScalarData{N,Float64}
    Ff::Edges{Primal, NX, NY, Float64}
    Ww::Edges{Dual, NX, NY,Float64}
    Qq::Edges{Dual, NX, NY,Float64}
    Fc::Nodes{Primal, NX, NY,Float64}
    Fn::Nodes{Dual, NX, NY,Float64}

    # Flags
    _isstored :: Bool

end

function NavierStokes(Re::Real, Δx::Real, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt::Real;
                       freestream::F = (0.0, 0.0), bodies::Union{BodyList,Nothing} = nothing,
                       motions::Union{RigidMotionList,Nothing} = nothing,
                       store_operators = true,
                       static_points = true,
                       rk::ConstrainedSystems.RKParams=ConstrainedSystems.RK31,
                       ddftype=CartesianGrids.Yang3) where {F}

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    α = Δt/(Re*Δx^2)

    Ff = Edges{Primal,NX,NY,Float64}()
    Ww = Edges{Dual, NX, NY,Float64}()
    Qq = Edges{Dual, NX, NY,Float64}()
    Fc = Nodes{Primal,NX,NY,Float64}()
    Fn = Nodes{Dual,NX,NY,Float64}()

    L = plan_laplacian(Fn,with_inverse=true)
    Lc = plan_laplacian(Fc,with_inverse=true)

    # Regularization and interpolation operators assumed empty
    Rf = nothing
    Ef = nothing
    Rc = nothing
    Ec = nothing
    Rn = nothing
    En = nothing
    Cf = nothing
    df = nothing
    sn = nothing
    sc = nothing

    points = isnothing(bodies) ? VectorData(0) : VectorData(collect(bodies))
    Vb = VectorData(points)
    Sb = ScalarData(points)
    N = length(points)÷NDIM


    if N > 0 && store_operators && static_points
      # in this case, points are assumed to be in inertial coordinates

      body_areas = areas(bodies)
      body_normals = normals(bodies)

      df = DoubleLayer(bodies,g,Ff)
      sc = SingleLayer(bodies,g,Fc)
      sn = SingleLayer(bodies,g,Fn)

      regop = Regularize(points,Δx;I0=CartesianGrids.origin(g),weights=body_areas.data,ddftype=ddftype)
      Rf = RegularizationMatrix(regop,Vb,Ff)
      Ef = InterpolationMatrix(regop,Ff,Vb)
      Rc = RegularizationMatrix(regop,Sb,Fc)
      Ec = InterpolationMatrix(regop,Fc,Sb)
      Rn = RegularizationMatrix(regop,Sb,Fn)
      En = InterpolationMatrix(regop,Fn,Sb)

      regopfilt = Regularize(points,Δx;I0=CartesianGrids.origin(g),filter=true,weights=Δx^2,ddftype=ddftype)
      Ẽf = InterpolationMatrix(regopfilt,Ff,Vb)
      Cf = sparse(Ẽf*Rf)

    end

    NavierStokes{NX, NY, N, _motiontype(static_points), _fstype(F)}(
                          Re, freestream, bodies, motions,
                          g, Δt, rk,
                          L, Lc,
                          df,sc,sn,
                          points, Rf, Ef, Rc, Ec, Rn, En, Cf,
                          Vb, Sb, Ff, Ww, Qq, Fc, Fn,
                          store_operators)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList;
        motions=RigidMotionList(map(x -> RigidBodyMotion(complex(0.0),0.0),bodies)),kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt;bodies=bodies,motions=motions,kwargs...)

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body;kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]);kwargs...)

function NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList,motions::RigidMotionList;kwargs...)
    length(bodies) == length(motions) || error("Inconsistent lengths of bodies and motions lists")
    NavierStokes(Re,Δx,xlim,ylim,Δt,bodies;motions=motions,kwargs...)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,motion::RigidBodyMotion;kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),RigidMotionList([motion]);kwargs...)


function Base.show(io::IO, sys::NavierStokes{NX,NY,N,MT,FS}) where {NX,NY,N,MT,FS}
    mtype = (MT == StaticPoints) ? "static" : "moving"
    fsmsg = (FS == StaticFreestream) ? "Static freestream = $(sys.U∞)" : "Variable freestream"
    bdmsg = (length(sys.bodies) == 1) ? "1 body" : "$(length(sys.bodies)) bodies"
    println(io, "Navier-Stokes system on a grid of size $NX x $NY and $N $mtype immersed points")
    println(io, "   $fsmsg")
    println(io, "   $bdmsg")
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
_motiontype(isstatic::Bool) = isstatic ? StaticPoints : MovingPoints
_fstype(F) = F == RigidBodyMotion ? VariableFreestream : StaticFreestream



include("navierstokes/fields.jl")
include("navierstokes/pulses.jl")
include("navierstokes/basicoperators.jl")
include("navierstokes/rigidbodyoperators.jl")
include("navierstokes/movingbodyoperators.jl")
