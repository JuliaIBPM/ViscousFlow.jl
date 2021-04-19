### Computing fields and other physical quantities ###
using LinearAlgebra
import LinearAlgebra: transpose!
import CartesianGrids: convective_derivative!

"""
    vorticity(sol::ODESolution,sys::NavierStokes,t)

Return the vorticity field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
@inline vorticity(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} =
            vorticity!(sys.Wn,w,sys,t)

function vorticity!(w::Nodes{Dual,NX,NY},win::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  w .= 0.0
  _vorticity!(w,win,sys,t)
  _vorticity_sheet!(w,sys,t)
  return w
end

@inline _vorticity!(w::Nodes{Dual,NX,NY},win::Nodes{Dual,NX,NY},
                    sys::NavierStokes{NX,NY},t::Real) where {NX,NY} = w .+= win/cellsize(sys)

_vorticity_sheet!(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,0},t::Real) where {NX,NY} = w

_vorticity_sheet!(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,
                                              ExternalInternalFlow},t::Real) where {NX,NY,N,MT,FS} = w

function _vorticity_sheet!(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,SD},t::Real) where {NX,NY,N,MT,FS,SD}
  if isnothing(sys.Rn)
    regop = _regularization(sys)
    Rn = RegularizationMatrix(regop,sys.Sb,sys.Sn)
  else
    Rn = sys.Rn
  end
  surface_velocity_jump!(sys.Δus,sys,t)
  cross!(sys.Sb,sys.Δus,normals(sys))
  w .+= Rn*sys.Sb
  return w
end


function velocity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    u .= 0.0
    _velocity_vorticity!(u,w,sys)
    _velocity_single_layer!(u,sys,t)
    _velocity_freestream!(u,sys,t)
    return u
end

function _velocity_vorticity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    sys.Wn .= 0.0
    _unscaled_streamfunction_vorticity!(sys.Wn,w,sys)
    sys.Vf .= 0.0
    curl!(sys.Vf,sys.Wn)
    u .+= sys.Vf
    return u
end

@inline _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,
                                      ExternalInternalFlow},t::Real) where {NX,NY,N,MT,FS} = u

@inline _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,0},t::Real) where {NX,NY} = u


function _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,SD},t::Real) where {NX,NY,N,MT,FS,SD}
    surface_velocity_jump!(sys.Δus,sys,t)
    pointwise_dot!(sys.Sb,sys.Δus,normals(sys))
    sys.Sb .*= cellsize(sys)
    sys.slc(sys.Sc,sys.Sb)
    _unscaled_scalarpotential!(sys.Sc,sys.Sc,sys)
    sys.Vf .= 0.0
    grad!(sys.Vf,sys.Sc)
    u .+= sys.Vf
    return u
end

function _velocity_freestream!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,SD},t::Real) where {NX,NY,N,MT,FS,SD}
  U∞, V∞ = freestream(t,sys)
  u.u .+= U∞
  u.v .+= V∞
  return u
end

"""
    velocity(sol::ODESolution,sys::NavierStokes,t)

Return the velocity field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
velocity(w::Nodes{Dual,NX,NY},a...) where {NX,NY} = velocity!(Edges(Primal,(NX,NY)),w,a...)


function streamfunction!(ψ::Nodes{Dual,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  ψ .= 0.0
  _unscaled_streamfunction_vorticity!(ψ,w,sys)
  ψ .*= cellsize(sys)
  _streamfunction_freestream!(ψ,sys,t)
  return ψ
end

function _unscaled_streamfunction_vorticity!(ψ::Nodes{Dual,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY}
  ldiv!(sys.Sn,sys.L,w)
  ψ .-= sys.Sn
  return ψ
end

function _streamfunction_freestream!(ψ::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  xg, yg = coordinates(ψ,sys.grid)
  U∞, V∞ = freestream(t,sys)
  ψ .+= U∞*yg' .- V∞*xg
  return ψ
end

"""
    streamfunction(sol::ODESolution,sys::NavierStokes,t)

Return the streamfunction field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
streamfunction(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} = streamfunction!(Nodes(Dual,(NX,NY)),w,sys,t)


function scalarpotential!(ϕ::Nodes{Primal,NX,NY},dil::Nodes{Primal,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  _unscaled_scalarpotential!(ϕ,dil,sys)
  # This assumes that dil is the physical dilation times the grid spacing
  ϕ .*= cellsize(sys)
  return ϕ
end

function _unscaled_scalarpotential!(ϕ::Nodes{Primal,NX,NY},dil::Nodes{Primal,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY}
  ldiv!(ϕ,sys.L,dil)
  return ϕ
end

"""
    scalarpotential(sol::ODESolution,sys::NavierStokes,t)

Return the scalar potential field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
scalarpotential(dil::Nodes{Primal,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} = scalarpotential!(Nodes(Primal,(NX,NY)),dil,sys,t)


function convective_derivative!(out::Edges{Primal,NX,NY},u::Edges{Primal,NX,NY},
                        sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  out .= u
  _unscaled_convective_derivative!(out,sys)
  Δx⁻¹ = 1/cellsize(sys)
  out .*= Δx⁻¹
  return out
end

"""
    convective_derivative(sol::ODESolution,sys::NavierStokes,t)

Return the convective derivative associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
convective_derivative(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} =
      convective_derivative!(Edges(Primal,(NX,NY)),u,sys,t)

"""
    pressure(sol::ODESolution,sys::NavierStokes,t)

Return the pressure field associated with solution vector `sol` on the grid in `sys`,
at time(s) `t`.
"""
function pressure(u::ConstrainedSystems.ArrayPartition,sys::NavierStokes{NX,NY},t) where {NX,NY}
    w, τ = state(u), constraint(u)
    fill!(sys.Vn,0.0)
    _vel_ns_rhs_convectivederivative!(sys.Vn,w,sys,t)
    _vel_ns_rhs_double_layer!(sys.Vn,sys,t)
    fill!(sys.Vf,0.0)
    _vel_ns_op_constraint_force!(sys.Vf,τ,sys)
    sys.Vn .-= sys.Vf
    fill!(sys.Sc,0.0)
    divergence!(sys.Sc,sys.Vn)
    press = zero(sys.Sc)
    ldiv!(press,sys.L,sys.Sc)
    press .*= cellsize(sys)
    return press
end


"""
    vorticity(integrator)

Return the vorticity field associated with `integrator` at its current state.
""" function vorticity end

"""
    velocity(integrator)

Return the velocity field associated with `integrator` at its current state.
""" function velocity end

"""
    streamfunction(integrator)

Return the streamfunction field associated with `integrator` at its current state.
""" function streamfunction end

"""
    scalarpotential(integrator)

Return the scalar potential field associated with `integrator` at its current state.
""" function scalarpotential end

"""
    convective_derivative(integrator)

Return the convective derivative associated with `integrator` at its current state.
""" function convective_derivative end

"""
    pressure(integrator)

Return the pressure field associated with `integrator` at its current state.
""" function pressure end


for fcn in (:vorticity,:velocity,:streamfunction,:scalarpotential,:convective_derivative)

  @eval $fcn(s::ConstrainedSystems.ArrayPartition,sys::NavierStokes,t) = $fcn(state(s),sys,t)

end

for fcn in (:vorticity,:velocity,:streamfunction,:scalarpotential,:convective_derivative,:pressure)
  @eval $fcn(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator) = $fcn(integ.u,integ.p,integ.t)

  @eval $fcn(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::NavierStokes,t) = $fcn(sol(t),sys,t)

  @eval $fcn(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::NavierStokes,t::AbstractArray) =
      map(ti -> $fcn(sol(ti),sys,ti),t)

end


# Surface field quantities

"""
    traction(integrator)

Return the traction at all of the surface points associated with `integrator` at
its current state.
""" function traction end

"""
    pressurejump(integrator)

Return the pressure jump at all of the surface points associated with `integrator` at
its current state.
""" function pressurejump end


for fcn in (:traction,:pressurejump)
  @eval $fcn(s::ConstrainedSystems.ArrayPartition,a...;kwargs...) = $fcn(constraint(s),a...;kwargs...)
  @eval $fcn(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator) = $fcn(integ.u,integ.p,integ.t)
  @eval $fcn(sol::ODESolution,sys::NavierStokes,t) = $fcn(sol(t),sys,t)

  @eval $fcn(sol::ODESolution,sys::NavierStokes,t::AbstractArray) =
        map(ti -> $fcn(sol(ti),sys,ti),t)

end

"""
    traction(sol::ODESolution,sys::NavierStokes,t)

Return the traction at all of the surface points associated with solution `sol` on
grid in `sys` at time(s) t.
"""
traction(τ,sys::NavierStokes{NX,NY,0,MT,FS,SD},t) where {NX,NY,MT,FS,SD} = τ

traction(τ::VectorData{N},sys::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow},t) where {NX,NY,N,MT,FS} = τ

function traction(τ::VectorData{N},sys::NavierStokes{NX,NY,N,MT,FS,SD},t) where {NX,NY,N,MT,FS,SD}
    nrm = normals(sys.bodies)
    relative_surface_velocity!(sys.Δus,sys,t)
    pointwise_dot!(sys.Sb,nrm,sys.Δus)
    surface_velocity_jump!(sys.Δus,sys,t)
    product!(sys.Vb,sys.Δus,sys.Sb)
    return τ + sys.Vb
end

"""
    pressurejump(sol::ODESolution,sys::NavierStokes,t)

Return the pressure jump at all of the surface points associated with solution `sol` on
grid in `sys` at time(s) t.
"""
function pressurejump(τ::VectorData{N},sys::NavierStokes{NX,NY,N},t) where {NX,NY,N}

    #_hasfilter(sys) ? (τf = similar(τ); τf .= sys.Cf*force(τ,sys)) : τf = deepcopy(force(τ,sys))
    #τf = deepcopy(τ)
    #return nrm.u∘τf.u + nrm.v∘τf.v # This might need some modification

    press = zero(sys.Sb)
    nrm = normals(sys.bodies)
    pointwise_dot!(press,nrm,traction(τ,sys,t))
    return -press
end

## Total quantities

force(τ,sys::NavierStokes{NX,NY,0},t,bodyi::Int) where {NX,NY} = Vector{Float64}(), Vector{Float64}()

function force(τ::VectorData{N},sys::NavierStokes,t,bodyi::Int) where {N}
    product!(sys.τ,traction(τ,sys,t),areas(sys.bodies))
    fxvec = sum(sys.τ.u,sys.bodies,bodyi)
    fyvec = sum(sys.τ.v,sys.bodies,bodyi)
    return fxvec, fyvec
end


force(s::ConstrainedSystems.ArrayPartition,sys,t,bodyi) = force(constraint(s),sys,t,bodyi)

"""
    force(integ,bodyindex)

Given the state of the system in integrator `integ`, return the current force
on the body with index `bodyindex`.
"""
force(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator,bodyi) = force(integ.u,integ.p,integ.t,bodyi)

"""
    force(sol,sys::NavierStokes,bodyindex) -> Tuple{Vector,Vector}

Given the solution history vector `sol` and the system `sys`, return the force
history on the body with index `bodyindex` as a tuple of vectors, one for
each component.
"""
function force(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::NavierStokes,bodyi::Int)
    fx = map((u,t) -> force(u,sys,t,bodyi)[1],sol.u,sol.t)
    fy = map((u,t) -> force(u,sys,t,bodyi)[2],sol.u,sol.t)
    fx, fy
end


function moment(τ::VectorData{N},bodies::BodyList,sys::NavierStokes,t,bodyi::Int;center=(0.0,0.0)) where {N}
    xc, yc = center
    x, y = collect(bodies)
    dX = VectorData(x .- xc,  y .- yc)
    product!(sys.τ,traction(τ,sys,t),areas(bodies))
    return sum(cross(dX,sys.τ),bodies,bodyi)
end


function moment(s::ConstrainedSystems.ArrayPartition,sys::NavierStokes{NX,NY,N,MovingPoints},t,bodyi;center=(0.0,0.0)) where {NX,NY,N}
    bodies = deepcopy(sys.bodies)
    tl! = RigidTransformList(aux_state(s))
    tl!(bodies)
    moment(constraint(s),bodies,sys,t,bodyi;center=center)
end

moment(s::ConstrainedSystems.ArrayPartition,sys::NavierStokes{NX,NY,N,StaticPoints},t,bodyi;center=(0.0,0.0)) where {NX,NY,N} =
        moment(constraint(s),sys.bodies,sys,t,bodyi;center=center)


"""
    moment(integ,bodyindex[,center=(0,0))

Given the state of the system in integrator `integ`, return the current moment
on the body with index `bodyindex`.  The center of the
moment can be passed as an optional tuple argument.
"""
moment(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator,bodyi;center=(0.0,0.0)) =
      moment(integ.u,integ.p,integ.t,bodyi,center=center)


"""
    moment(sol,sys::NavierStokes,bodyindex[,center=(0,0)]) -> Vector

Given the solution history vector `sol` and the system `sys`, return the moment
history on the body with index `bodyindex` as a vector. The center of the
moment can be passed as an optional tuple argument.
"""
function moment(sol::ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::NavierStokes,bodyi::Int;center=(0.0,0.0))
    mom = map((u,t) -> moment(u,sys,t,bodyi,center=center),sol.u,sol.t)
end
