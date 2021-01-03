### Computing fields and other physical quantities ###
using LinearAlgebra
import LinearAlgebra: transpose!
import CartesianGrids: convective_derivative!

"""
    vorticity(u,sys::NavierStokes,t)

Return the vorticity field associated with state vector `u` on the grid in `sys`,
at time `t`.
"""
@inline vorticity(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} = w/cellsize(sys)

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

function _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow},t::Real) where {NX,NY,N,MT,FS}
    return u
end

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
    velocity(u,sys::NavierStokes,t)

Return the velocity field associated with state vector `u` on the grid in `sys`,
at time `t`.
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
    streamfunction(u,sys::NavierStokes,t)

Return the streamfunction field associated with state vector `u` on the grid in `sys`,
at time `t`.
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
    scalarpotential(u,sys::NavierStokes,t)

Return the scalar potential field associated with state vector `u` on the grid in `sys`,
at time `t`.
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
    convective_derivative(u,sys::NavierStokes,t)

Return the convective derivative associated with state vector `u` on the grid in `sys`,
at time `t`.
"""
convective_derivative(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} =
      convective_derivative!(Edges(Primal,(NX,NY)),u,sys,t)


"""
    vorticity(integrator)

Return the vorticity field associated with `integrator` at its current state.
"""
function vorticity end

"""
    velocity(integrator)

Return the velocity field associated with `integrator` at its current state.
"""
function velocity end

"""
    streamfunction(integrator)

Return the streamfunction field associated with `integrator` at its current state.
"""
function streamfunction end

"""
    scalarpotential(integrator)

Return the scalar potential field associated with `integrator` at its current state.
"""
function scalarpotential end

"""
    convective_derivative(integrator)

Return the convective derivative associated with `integrator` at its current state.
"""
function convective_derivative end


for fcn in (:vorticity,:velocity,:streamfunction,:scalarpotential,:convective_derivative)
  @eval $fcn(s::ConstrainedSystems.ArrayPartition,sys::NavierStokes,t) = $fcn(state(s),sys,t)

  @eval $fcn(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator) = $fcn(integ.u,integ.p,integ.t)
end


function pressure(u::ConstrainedSystems.ArrayPartition,sys::NavierStokes{NX,NY},t) where {NX,NY}
    w, τ = state(u), constraint(u)
    sys.Vv .= 0
    _vel_ns_rhs_convectivederivative!(sys.Vv,w,sys)
    _vel_ns_rhs_double_layer!(sys.Vv,sys,t)
    sys.Vf .= 0
    _vel_ns_op_constraint_force!(sys.Vf,τ,sys)
    sys.Vv .-= sys.Vf
    sys.Sc .= 0
    divergence!(sys.Sc,sys.Vv)
    press = zero(sys.Sc)
    ldiv!(press,sys.L,sys.Sc)
    press .*= cellsize(sys)
    return press
end
pressure(integ::ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator) = pressure(integ.u,integ.p,integ.t)

# Surface field quantities

for fcn in (:force,:pressurejump)
  @eval $fcn(s::ConstrainedSystems.ArrayPartition,a...;kwargs...) = $fcn(constraint(s),a...;kwargs...)
end

force(τ::VectorData{N},sys) where {N} = τ #*cellsize(sys)^2

function pressurejump(τ::VectorData{N},b::Union{Body,BodyList},sys::NavierStokes) where {N}
    @assert N == numpts(b)

    #_hasfilter(sys) ? (τf = similar(τ); τf .= sys.Cf*force(τ,sys)) : τf = deepcopy(force(τ,sys))
    τf = deepcopy(τ)

    nrm = normals(b)
    return nrm.u∘τf.u + nrm.v∘τf.v # This might need some modification
end
