import Base: size

const NDIM = 2


import CartesianGrids: cellsize, origin
import RigidBodyTools: assign_velocity!
import ImmersedLayers: normals, areas

include("navierstokes/linesource.jl")

#=
The goal of this should be to change NavierStokes into
a more general type for any PDEs, similar to the structures
used for OrdinaryDiffEq.

Rule of thumb: keep the typing generic!

Other things to keep in mind:
- For MovingPoints problems, the immersion operators will get changed. These
   need to be easily accessible from the overall system. We can endow
   the term type with a parameter for types that must be updated if the points
   update. Each such type will need to be accompanied by an operator for
   setting it up (and possibly a separate one for updating it).
- Some operators, such as Rf and Ef, may get used by more than one term.
=#

"""
$(TYPEDEF)

A system type that utilizes a grid of `NX` x `NY` dual cells and `N` Lagrange forcing
points to solve the discrete Navier-Stokes equations in vorticity form. The
parameter `static_points` specifies whether the forcing points remain static in the
grid. It should be set to `false` if a supplied motion requires that the points move.

# Constructors:

`NavierStokes(Re,Δx,xlimits,ylimits,Δt
              [,freestream = (0.0, 0.0)]
              [,pulses=nothing])` specifies the Reynolds number `Re`, the grid
              spacing `Δx`, the dimensions of the domain in the tuples `xlimits`
              and `ylimits` (excluding the ghost cells), and the time step size `Δt`.
              The other arguments are optional. The `freestream` argument can be
              passed as either a tuple (a static freestream) or a `RigidBodyMotion`
              for a time-varying freestream. The `pulses` argument can be
              used to pass in one or more spatiotemporal pulses.


`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList
              [,flow_side=ExternalInternalFlow]
              [,ddftype=CartesianGrids.Yang3])` passes the body
              information. This constructor
              sets the motions of the body/ies to be stationary.
              The same optional arguments used for the basic constructor
              also apply for this one. In addition, the `flow_side` can be set to
              `ExternalFlow` (default), `InternalFlow`, or `ExternalInternalFlow`.
              However, it is forced to `ExternalInternalFlow` for open Bodies
              (like `Plate` type).

`NavierStokes(Re,Δx,xlimits,ylimits,Δt,bodies::Body/BodyList,
              motions::Motion/MotionList
              [,static_points=false])`
              passes the body and associated motion information.
              The list of motions must be the same length as the list of bodies.
              The same optional arguments used for the other constructors
              also apply for this one. In addition, `static_points` can
              be set to `true` if the supplied motion should not cause the
              points to move.

"""
mutable struct NavierStokes{NX, NY, N, MT<:PointMotionType, FS<:FreestreamType, SD<:FlowSide, DDF<:CartesianGrids.DDFType, FT,
                            SP, WP, VP, SPP, GVP, VBP, SBP, FST,
                            RHST, CDT, COT, DLT, VCT, PCT, SCT, SVCT, STCT}
    # Physical Parameters
    "Reynolds number"
    Re::Float64
    "Free stream velocities"
    U∞::FST
    "Bodies"
    bodies::Union{BodyList,Nothing}
    "Body motions"
    motions::Union{RigidMotionList,DirectlySpecifiedMotionList,Nothing}
    "Pulses"
    pulses::Union{Vector{ModulatedField},Nothing}
    
    ## ADD τ_bc HERE (PRESCRIBED LINE SOURCE)
    τ_bc::Union{Vector{PrescribedLineSource},Nothing}
    
    # Discretization
    "Grid metadata"
    grid::CartesianGrids.PhysicalGrid{2}
    "Time step"
    Δt::Float64

    # Operators
    "Laplacian operator"
    L::CartesianGrids.Laplacian

    # Layers
    dlf::Union{DoubleLayer,Nothing} # used for viscous surface terms
    slc::Union{SingleLayer,Nothing} # used for scalar potential field in velocity
    sln::Union{SingleLayer,Nothing} # might not be used

    # Body coordinate data, if present
    # if a static problem, these coordinates are in inertial coordinates
    # if a non-static problem, in their own coordinate systems
    points::VectorData{N,Float64}

    # Pre-stored regularization and interpolation matrices (if present)
    Rf::Union{RegularizationMatrix,Nothing} # faces (edges)
    Ef::Union{InterpolationMatrix,Nothing}
    Cf::Union{AbstractMatrix,Nothing}

    # Operators
    f :: FT

    # prototypes
    state_prototype :: SP
    vorticity_prototype :: WP
    velocity_prototype :: VP
    scapot_prototype :: SPP
    gradv_prototype :: GVP
    surfacevel_prototype :: VBP
    surfacescalar_prototype :: SBP

    # Caches
    rhs_cache :: RHST
    convderiv_cache :: CDT
    constraintop_cache :: COT
    doublelayer_cache :: DLT
    velocity_cache :: VCT
    pressure_cache :: PCT
    streamfunction_cache :: SCT
    surfacevel_cache :: SVCT
    surfacetraction_cache :: STCT

end

function NavierStokes(Re::Real, Δx::Real, xlimits::Tuple{Real,Real},ylimits::Tuple{Real,Real}, Δt::Real;
                       freestream::FST = (0.0, 0.0),
                       bodies::Union{BodyList,Nothing} = nothing,
                       motions::Union{RigidMotionList,DirectlySpecifiedMotionList,Nothing} = nothing,
                       pulses::PT = nothing,
                       τ_bc = nothing,
                       static_points = true,
                       flow_side::Type{SD} = ExternalFlow,
                       ddftype=CartesianGrids.Yang3) where {FST,PT,SD<:FlowSide}

    g = PhysicalGrid(xlimits,ylimits,Δx)
    NX, NY = size(g)

    N = numpts(bodies)

    α = Δt/(Re*Δx^2)

    vorticity_prototype = Nodes{Dual,NX,NY,Float64}()
    velocity_prototype = Edges{Primal,NX,NY,Float64}()
    scapot_prototype = Nodes{Primal,NX,NY,Float64}()
    gradv_prototype = EdgeGradient{Primal,Dual,NX,NY,Float64}()
    surfacevel_prototype = VectorData(N)
    surfacescalar_prototype = ScalarData(N)

    L = plan_laplacian(vorticity_prototype,with_inverse=true)

    pulsefields = _process_pulses(pulses,vorticity_prototype,g)


    # for now, if there are any bodies that are Open,
    # then force flow_side to ExternalInternalFlow.
    # but should be more flexible here
    flow_side_internal = _any_open_bodies(bodies) ? ExternalInternalFlow : flow_side


    points, dlf, slc, sln, Rf, Ef, Cf =
              _immersion_operators(bodies,g,flow_side_internal,ddftype,
                        velocity_prototype,scapot_prototype,surfacevel_prototype)


    viscous_L = plan_laplacian(similar(vorticity_prototype),
                               factor=1/(Re*Δx^2))

    rhs_cache = RHSCache(similar(velocity_prototype))

    convderiv_cache = ConvectiveDerivativeCache(similar(velocity_prototype),
                                                similar(gradv_prototype),
                                                similar(gradv_prototype),
                                                similar(gradv_prototype))
    constraintop_cache = ConstraintOperatorCache(similar(velocity_prototype),
                                                 similar(surfacevel_prototype))

    doublelayer_cache = DoubleLayerCache(similar(velocity_prototype),
                                         similar(surfacevel_prototype))

    velocity_cache = VelocityCache(similar(vorticity_prototype),
                                   similar(velocity_prototype),
                                   similar(scapot_prototype),
                                   VectorData(N),
                                   ScalarData(N))

    pressure_cache = PressureCache(similar(velocity_prototype),
                                   similar(velocity_prototype),
                                   similar(scapot_prototype))

    streamfunction_cache = StreamfunctionCache(similar(vorticity_prototype))

    surfacevel_cache = SurfaceVelocityCache(similar(surfacevel_prototype))

    surfacetraction_cache = SurfaceTractionCache(similar(surfacevel_prototype),
                                                similar(surfacevel_prototype),
                                                similar(surfacescalar_prototype))

    if isnothing(bodies)
      state_prototype = solvector(state=vorticity_prototype)
      f = ConstrainedODEFunction(ns_rhs!,viscous_L,_func_cache=state_prototype)
    else
      if static_points
        state_prototype = solvector(state=vorticity_prototype,constraint=surfacevel_prototype)
        f = ConstrainedODEFunction(ns_rhs!,bc_constraint_rhs!,
                                      ns_op_constraint_force!,bc_constraint_op!,
                                      viscous_L,_func_cache=state_prototype)
      else
        state_prototype = solvector(state=vorticity_prototype,constraint=surfacevel_prototype,aux_state=zero_body_state(bodies))
        rhs! = ConstrainedSystems.r1vector(state_r1 = ns_rhs!,aux_r1 = rigid_body_rhs!)
        f = ConstrainedODEFunction(rhs!,bc_constraint_rhs!,
                                      ns_op_constraint_force!,bc_constraint_op!,
                                      viscous_L,_func_cache=state_prototype,
                                      param_update_func=update_immersion_operators!)
      end

    end





    NavierStokes{NX, NY, N, _motiontype(static_points), _fstype(FST), flow_side_internal,
                  ddftype, typeof(f), typeof(state_prototype),typeof(vorticity_prototype),
                  typeof(velocity_prototype),typeof(scapot_prototype),typeof(gradv_prototype),
                  typeof(surfacevel_prototype),typeof(surfacescalar_prototype),
                  FST,typeof(rhs_cache),
                  typeof(convderiv_cache), typeof(constraintop_cache), typeof(doublelayer_cache),
                  typeof(velocity_cache), typeof(pressure_cache),
                  typeof(streamfunction_cache),typeof(surfacevel_cache),typeof(surfacetraction_cache)}(
                          Re, freestream, bodies, motions, pulsefields,
                          _line_source(τ_bc,Edges(Primal,size(g)),g),
                          g, Δt, # rk,
                          L,
                          dlf,slc,sln,
                          points, Rf, Ef, Cf,
                          f,state_prototype,
                          vorticity_prototype, velocity_prototype, scapot_prototype,
                          gradv_prototype, surfacevel_prototype, surfacescalar_prototype,
                          rhs_cache, convderiv_cache, constraintop_cache,
                          doublelayer_cache,
                          velocity_cache, pressure_cache, streamfunction_cache,
                          surfacevel_cache,surfacetraction_cache)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList;
        motions=RigidMotionList(map(x -> RigidBodyMotion(0.0,0.0),bodies)),kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt;bodies=bodies,motions=motions,kwargs...)

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),kwargs...)

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,τ_bc,kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),τ_bc=τ_bc,kwargs...)


function NavierStokes(Re,Δx,xlim,ylim,Δt,bodies::BodyList,motions::Union{RigidMotionList,DirectlySpecifiedMotionList};static_points=false,kwargs...)
    length(bodies) == length(motions) || error("Inconsistent lengths of bodies and motions lists")
    NavierStokes(Re,Δx,xlim,ylim,Δt,bodies;motions=motions,static_points=static_points,kwargs...)
end

NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,motion::RigidBodyMotion;static_points=false,kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),RigidMotionList([motion]);static_points=static_points,kwargs...)


NavierStokes(Re,Δx,xlim,ylim,Δt,body::Body,motion::DirectlySpecifiedMotion;static_points=false,kwargs...) =
        NavierStokes(Re,Δx,xlim,ylim,Δt,BodyList([body]),DirectlySpecifiedMotionList([motion]);static_points=static_points,kwargs...)


function Base.show(io::IO, sys::NavierStokes{NX,NY,N,MT,FS,SD}) where {NX,NY,N,MT,FS,SD}
    mtype = (MT == StaticPoints) ? "static" : "moving"
    fsmsg = (FS == StaticFreestream) ? "Static freestream = $(sys.U∞)" : "Variable freestream"
    sdmsg = (N == 0) ? "Unbounded" : ((SD == ExternalFlow) ? "External flow" : ((SD == InternalFlow) ? "Internal flow" : "External/internal"))
    println(io, "$sdmsg Navier-Stokes system on a grid of size $NX x $NY and $N $mtype immersed points")
    println(io, "   $fsmsg")
    if N > 0
      bdmsg = (length(sys.bodies) == 1) ? "1 body" : "$(length(sys.bodies)) bodies"
      println(io, "   $bdmsg")
    end
end

# Routines to set up the immersion operators

function _immersion_operators(bodies::BodyList,g::PhysicalGrid,flow_side::Type{SD},ddftype,
                              gridvector_prototype::GridData{NX,NY},gridscalar_prototype::GridData{NX,NY},
                              surfvector_prototype::PointData{N}) where {NX,NY,N,SD<:ViscousFlow.FlowSide}

  points = VectorData(collect(bodies))
  numpts(bodies) == N || error("Inconsistent size of bodies")

  body_areas = areas(bodies)
  body_normals = normals(bodies)

  if !(flow_side==ExternalInternalFlow)
    dlf = DoubleLayer(bodies,g,gridvector_prototype)
    slc = SingleLayer(bodies,g,gridscalar_prototype)
    #sln = SingleLayer(bodies,g,Sn)
    sln = nothing
  else
    dlf = nothing
    slc = nothing
    sln = nothing
  end

  #regop = Regularize(points,cellsize(g);I0=CartesianGrids.origin(g),weights=body_areas.data,ddftype=ddftype)
  regop = _regularization(points,g,bodies,ddftype)

  Rf = RegularizationMatrix(regop,surfvector_prototype,gridvector_prototype) # Used by B₁ᵀ
  Ef = InterpolationMatrix(regop,gridvector_prototype,surfvector_prototype) # Used by constraint_rhs! and B₂

  regopfilt = Regularize(points,cellsize(g);I0=CartesianGrids.origin(g),filter=true,weights=cellsize(g)^2,ddftype=ddftype)

  Ẽf = InterpolationMatrix(regopfilt,gridvector_prototype,surfvector_prototype)
  Cf = sparse(Ẽf*Rf)

  return points, dlf, slc, sln, Rf, Ef, Cf
end

_immersion_operators(::Nothing,a...) =
    VectorData(0), nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing


# For updating the system with body data

function update_immersion_operators!(sys::NavierStokes{NX,NY,N,MT,FS,SD,DDF},bodies::BodyList) where {NX,NY,N,MT,FS,SD,DDF<:CartesianGrids.DDFType}
    @unpack velocity_prototype, scapot_prototype, surfacevel_prototype = sys
    sys.bodies = deepcopy(bodies)
    sys.points, sys.dlf, sys.slc, sys.sln, sys.Rf, sys.Ef, sys.Cf =
      _immersion_operators(sys.bodies,sys.grid,SD,DDF,velocity_prototype,scapot_prototype,surfacevel_prototype)
    return sys
end

update_immersion_operators!(sys::NavierStokes,body::Body) = update_immersion_operators!(sys,BodyList([body]))

function update_immersion_operators!(sys::NavierStokes,x::AbstractVector)
    tl! = RigidTransformList(x)
    tl!(sys.bodies)
    update_immersion_operators!(sys,sys.bodies)
end

# The form passed to ConstrainedODEFunction
update_immersion_operators!(sys::NavierStokes,u,sys_old::NavierStokes,t) =
    update_immersion_operators!(sys,aux_state(u))




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
    Δt = min(cfl*Δx,fourier*Δx^2*Re)
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
    freestream(t,sys::NavierStokes) -> Tuple

Return the value of the freestream in `sys` at time `t` as a tuple.
"""
freestream(t::Real,sys::NavierStokes{NX,NY,N,MT,StaticFreestream}) where {NX,NY,N,MT} = sys.U∞

function freestream(t::Real,sys::NavierStokes{NX,NY,N,MT,VariableFreestream}) where {NX,NY,N,MT}
    _,ċ,_,_,_,_ = sys.U∞(t)
    return reim(ċ)
end

"""
    newstate(sys::NavierStokes)

Return a new (zero) instance of the state vector for `sys`.
"""
newstate(sys::NavierStokes) = zero(sys.state_prototype)

"""
    newstate(s::AbstractSpatialField,sys::NavierStokes)

Return an instance of the state vector for `sys`, assigned the
data in the spatial field `s`.
"""
function newstate(s::AbstractSpatialField,sys::NavierStokes)
  u = newstate(sys)
  gf = GeneratedField(state(u),s,sys.grid)
  state(u) .= cellsize(sys)*gf()
  return u
end

"""
    flowside(sys::NavierStokes) -> FlowSide

Returns the side of the body on which the flow is to be computed.
"""
flowside(::NavierStokes{NX,NY,N,MT,FS,SD}) where {NX,NY,N,MT,FS,SD} = SD

# Other functions
_hasfilter(sys::NavierStokes) = !isnothing(sys.Cf)
_motiontype(isstatic::Bool) = isstatic ? StaticPoints : MovingPoints
  _fstype(F) = F <: Union{RigidBodyMotion,Kinematics} ? VariableFreestream : StaticFreestream

_body_closure_type(b::T) where {T<:Body{N,C}} where {N,C} = C

_any_open_bodies(nothing) = false
_any_open_bodies(bodies::BodyList) =  any(b -> _body_closure_type(b) == RigidBodyTools.OpenBody,bodies)

_regularization(sys::NavierStokes{NX, NY, N, MT, FS, SD, DDF}) where {NX,NY,N,MT,FS,SD,DDF} =
        _regularization(sys.points,sys.grid,sys.bodies,DDF)

_regularization(points,g,bodies,ddftype) = Regularize(points,cellsize(g),
                                I0=CartesianGrids.origin(g),weights=areas(bodies).data,ddftype=ddftype)


function _line_source(lineparams::Vector{<:LineSourceParams},u::VectorGridData,g::PhysicalGrid)
  τ_bc = PrescribedLineSource[]
  for p in lineparams
    push!(τ_bc,PrescribedLineSource(p,u,g))
  end
  τ_bc
end

_line_source(lineparams::LineSourceParams,u::VectorGridData,g::PhysicalGrid) =
    _line_source([lineparams],u,g)

_line_source(::Nothing,u,g) = nothing



function set_linesource_strength!(sys::NavierStokes,q::Vector{T}) where {T<:Real}
    @unpack τ_bc = sys
    if isnothing(τ_bc)
      return sys
    else
      total_len = mapreduce(τl -> length(τl.τ),+,τ_bc)
      @assert total_len == length(q)
      qcpy = copy(q)
      first_index = 0
      for τl in τ_bc
        τl.τ .= view(q,first_index+1:first_index+length(τl.τ))
        first_index += length(τl.τ)
      end
    end
    return sys
end


include("navierstokes/surfacevelocities.jl")
include("navierstokes/fields.jl")
include("navierstokes/pointforce.jl")
include("navierstokes/basicoperators.jl")
include("navierstokes/rigidbodyoperators.jl")
include("navierstokes/movingbodyoperators.jl")
include("navierstokes/timemarching.jl")
