### Basic operators for any Navier-Stokes system ###
# Note that the input vector `w` is vorticity x cell spacing (i.e., it
# has units of velocity), and
# we expect the rhs of the equations to have units of dw/dt

export ns_rhs!


# RHS of Navier-Stokes equations
function ns_rhs!(dw::Nodes{Dual,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  @unpack rhs_cache = sys
  @unpack dv = rhs_cache
  fill!(dw,0.0)
  ns_rhs_velocity!(dv,w,sys,t)
  curl!(dw,dv)
  _ns_rhs_pulses!(dw,sys,t)
  return dw
end

function ns_rhs_velocity!(dv::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  fill!(dv,0.0)
  _vel_ns_rhs_convectivederivative!(dv,w,sys,t)
  _vel_ns_rhs_double_layer!(dv,sys,t)
  return dv
end

function _vel_ns_rhs_convectivederivative!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t) where {NX,NY}
    @unpack convderiv_cache = sys
    @unpack V = convderiv_cache
    Δx⁻¹ = 1/cellsize(sys)
    fill!(V,0.0)
    velocity!(V,w,sys,t)
    _unscaled_convective_derivative!(V,sys)
    V .*= Δx⁻¹
    u .-= V
end

# Operates in-place on `u`, which comes in with the velocity field and
# returns the unscaled convective derivative
function _unscaled_convective_derivative!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY}
    @unpack convderiv_cache = sys
    @unpack Vtf, DVf, VDVf = convderiv_cache
    transpose!(Vtf,grid_interpolate!(DVf,u))
    DVf .= 0.0
    grad!(DVf,u)
    product!(VDVf,Vtf,DVf)
    u .= 0.0
    grid_interpolate!(u,VDVf)
    u
end


_ns_rhs_pulses!(dw::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t) where {NX,NY} = _ns_rhs_pulses!(dw,sys.pulses,cellsize(sys),t)

_ns_rhs_pulses!(dw,::Nothing,Δx,t) = dw

function _ns_rhs_pulses!(dw,pulses::Vector{<:ModulatedField},Δx,t)
  for p in pulses
    dw .+= Δx*p(t)
  end
  dw
end


@inline _vel_ns_rhs_double_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,
                                  ExternalInternalFlow},t::Real) where {NX,NY,N,MT,FS} = u

@inline _vel_ns_rhs_double_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,0},t::Real) where {NX,NY} = u

function _vel_ns_rhs_double_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,SD},t::Real) where {NX,NY,N,MT,FS,SD}
    Δx⁻¹ = 1/cellsize(sys)
    fact = Δx⁻¹/sys.Re
    surface_velocity_jump!(sys.Δus,sys,t)
    fill!(sys.Vf,0.0)
    sys.dlf(sys.Vf,sys.Δus)
    sys.Vf .*= fact
    u .-= sys.Vf
end
