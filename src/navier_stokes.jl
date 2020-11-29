import Base: size

const NDIM = 2

import ConstrainedSystems: r₁, r₂, B₁ᵀ, B₂, plan_constraints
import CartesianGrids: cellsize, origin
import RigidBodyTools: assign_velocity!
import ImmersedLayers: normals, areas

"""
$(TYPEDEF)

A system type that utilizes a grid of `NX` x `NY` dual cells and `N` Lagrange forcing
points to solve the discrete Navier-Stokes equations in vorticity form. The
parameter `static_points` specifies whether the forcing points remain static in the
grid. It should be set to `false` if a supplied motion requires that the points move.

# Constructors:

`NavierStokes(Re,Δx,xlimits,ylimits,Δt
              [,freestream = (0.0, 0.0)]
              [,store_operators=true][,static_points=true]
              [,rk=ConstrainedSystems.RK31]
              [,flow_side=ExternalInternalFlow]
              [,ddftype=CartesianGrids.Yang3])`specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. Note that `store_operators` set to `true`
              stores matrix versions of the operators. This makes the method
              faster, at the cost of storage. The `freestream` argument can be
              passed as either a tuple (a static freestream) or a `RigidBodyMotion`
              for a time-varying freestream. The `flow_side` can be set to
              `ExternalFlow`, `InternalFlow`, or `ExternalInternalFlow` (default).

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList)` passes the body
              information. The other keywords can also be supplied. This
              sets the motions of the body/ies to be stationary.

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList,
              motions::RigidBodyMotion/RigidMotionList)`
              passes the body and associated motion information.
              The other keywords can also be supplied. The list of
              motions must be the same length as the list of bodies.

"""
mutable struct NavierStokes{NX, NY, N, MT<:PointMotionType, FS<:FreestreamType, SD<:FlowSide}
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

    # Layers
    dlf::Union{DoubleLayer,Nothing}
    slc::Union{SingleLayer,Nothing}
    sln::Union{SingleLayer,Nothing}

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
    Δus::VectorData{N,Float64}
    τ::VectorData{N,Float64}
    Vf::Edges{Primal, NX, NY, Float64}
    Vv::Edges{Primal, NX, NY, Float64}
    Sc::Nodes{Primal, NX, NY,Float64}
    Sn::Nodes{Dual, NX, NY,Float64}
    Wn::Nodes{Dual, NX, NY,Float64}
    Vtf::EdgeGradient{Primal,Dual,NX,NY,Float64}
    DVf::EdgeGradient{Primal,Dual,NX,NY,Float64}
    VDVf::EdgeGradient{Primal,Dual,NX,NY,Float64}

    # Flags
    _isstored :: Bool

end

function NavierStokes(Re::Real, Δx::Real, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt::Real;
                       freestream::F = (0.0, 0.0), bodies::Union{BodyList,Nothing} = nothing,
                       motions::Union{RigidMotionList,Nothing} = nothing,
                       store_operators = true,
                       static_points = true,
                       flow_side::Type{SD} = ExternalInternalFlow,
                       rk::ConstrainedSystems.RKParams=ConstrainedSystems.RK31,
                       ddftype=CartesianGrids.Yang3) where {F,SD<:FlowSide}

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    α = Δt/(Re*Δx^2)

    Vf = Edges{Primal,NX,NY,Float64}()
    Vv = Edges{Primal,NX,NY,Float64}()
    Sc = Nodes{Primal,NX,NY,Float64}()
    Sn = Nodes{Dual,NX,NY,Float64}()
    Wn = Nodes{Dual,NX,NY,Float64}()
    Vtf = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    DVf = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    VDVf = EdgeGradient{Primal,Dual,NX,NY,Float64}()

    L = plan_laplacian(Sn,with_inverse=true)

    # Regularization and interpolation operators assumed empty
    Rf = nothing
    Ef = nothing
    Rc = nothing
    Ec = nothing
    Rn = nothing
    En = nothing
    Cf = nothing
    dlf = nothing
    sln = nothing
    slc = nothing

    points = isnothing(bodies) ? VectorData(0) : VectorData(collect(bodies))
    Vb = VectorData(points)
    Sb = ScalarData(points)
    Δus = VectorData(points)
    τ = VectorData(points)
    N = length(points)÷NDIM


    if N > 0 && store_operators && static_points
      # in this case, points are assumed to be in inertial coordinates

      body_areas = areas(bodies)
      body_normals = normals(bodies)

      if !(flow_side==ExternalInternalFlow)
        dlf = DoubleLayer(bodies,g,Vf)
        slc = SingleLayer(bodies,g,Sc)
        sln = SingleLayer(bodies,g,Sn)
      end

      regop = Regularize(points,Δx;I0=CartesianGrids.origin(g),weights=body_areas.data,ddftype=ddftype)
      # May not need all of these...
      Rf = RegularizationMatrix(regop,Vb,Vf) # Used by B₁ᵀ
      Ef = InterpolationMatrix(regop,Vf,Vb) # Used by constraint_rhs! and B₂
      Rc = RegularizationMatrix(regop,Sb,Sc)
      Ec = InterpolationMatrix(regop,Sc,Sb)
      Rn = RegularizationMatrix(regop,Sb,Sn)
      En = InterpolationMatrix(regop,Sn,Sb)

      regopfilt = Regularize(points,Δx;I0=CartesianGrids.origin(g),filter=true,weights=Δx^2,ddftype=ddftype)
      Ẽf = InterpolationMatrix(regopfilt,Vf,Vb)
      Cf = sparse(Ẽf*Rf)

    end

    NavierStokes{NX, NY, N, _motiontype(static_points), _fstype(F), flow_side}(
                          Re, freestream, bodies, motions,
                          g, Δt, rk,
                          L,
                          dlf,slc,sln,
                          points, Rf, Ef, Rc, Ec, Rn, En, Cf,
                          Vb, Sb, Δus, τ,
                          Vf, Vv, Sc, Sn, Wn, Vtf, DVf, VDVf,
                          store_operators)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList;
        motions=RigidMotionList(map(x -> RigidBodyMotion(0.0,0.0),bodies)),kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt;bodies=bodies,motions=motions,kwargs...)

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body;kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]);kwargs...)

function NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList,motions::RigidMotionList;kwargs...)
    length(bodies) == length(motions) || error("Inconsistent lengths of bodies and motions lists")
    NavierStokes(Re,Δx,xlim,ylim,Δt,bodies;motions=motions,kwargs...)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,motion::RigidBodyMotion;kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),RigidMotionList([motion]);kwargs...)


function Base.show(io::IO, sys::NavierStokes{NX,NY,N,MT,FS,SD}) where {NX,NY,N,MT,FS,SD}
    mtype = (MT == StaticPoints) ? "static" : "moving"
    fsmsg = (FS == StaticFreestream) ? "Static freestream = $(sys.U∞)" : "Variable freestream"
    bdmsg = (length(sys.bodies) == 1) ? "1 body" : "$(length(sys.bodies)) bodies"
    sdmsg = (SD == ExternalFlow) ? "External flow" : ((SD == InternalFlow) ? "Internal flow" : "External/internal")
    println(io, "$sdmsg Navier-Stokes system on a grid of size $NX x $NY and $N $mtype immersed points")
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

# Wrap the output of the motion evaluation in VectorData
@inline assign_velocity!(u::VectorData,a...) = (assign_velocity!(u.u,u.v,a...); u)


@inline normals(sys::NavierStokes) = (!isnothing(sys.bodies)) ? normals(sys.bodies) : nothing

"""
    freestream(t,sys) -> Tuple

Return the value of the freestream in `sys` at time `t` as a tuple.
"""
freestream(t::Real,sys::NavierStokes{NX,NY,N,MT,StaticFreestream}) where {NX,NY,N,MT} = sys.U∞

function freestream(t::Real,sys::NavierStokes{NX,NY,N,MT,VariableFreestream}) where {NX,NY,N,MT}
    _,ċ,_,_,_,_ = sys.U∞(t)
    return reim(ċ)
end

# Other functions
_hasfilter(sys::NavierStokes) = !(sys.Cmat == nothing)
_motiontype(isstatic::Bool) = isstatic ? StaticPoints : MovingPoints
_fstype(F) = F == RigidBodyMotion ? VariableFreestream : StaticFreestream


include("navierstokes/surfacevelocities.jl")
include("navierstokes/fields.jl")
include("navierstokes/pulses.jl")
include("navierstokes/basicoperators.jl")
include("navierstokes/rigidbodyoperators.jl")
include("navierstokes/movingbodyoperators.jl")
