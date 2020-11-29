### Computing fields and other physical quantities ###
using LinearAlgebra
import LinearAlgebra: transpose!
import CartesianGrids: convective_derivative!

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

function _streamfunction_freestream(ψ::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  xg, yg = coordinates(ψ,sys.grid)
  U∞, V∞ = freestream(t,sys)
  ψ .+= U∞*yg' .- V∞*xg
  return ψ
end

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

scalarpotential(dil::Nodes{Primal,NX,NY},sys::NavierStokes{NX,NY},t::Real) where {NX,NY} = scalarpotential!(Nodes(Primal,(NX,NY)),dil,sys,t)


function convective_derivative!(out::Edges{Primal,NX,NY},u::Edges{Primal,NX,NY},
                        sys::NavierStokes{NX,NY},t::Real) where {NX,NY}
  out .= u
  _unscaled_convective_derivative!(out,sys)
  Δx⁻¹ = 1/cellsize(sys)
  out .*= Δx⁻¹
  return out
end


_unscaled_convective_derivative!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY} =
      _unscaled_convective_derivative!(u,sys.Vtf,sys.DVf,sys.VDVf)

# Operates in-place on `u`, which comes in with the velocity field and
# returns the unscaled convective derivative
function _unscaled_convective_derivative!(u::Edges{Primal,NX,NY},
                                          Vtf::EdgeGradient{Primal,Dual,NX,NY},
                                          DVf::EdgeGradient{Primal,Dual,NX,NY},
                                          VDVf::EdgeGradient{Primal,Dual,NX,NY}) where {NX,NY}
    transpose!(Vtf,grid_interpolate!(DVf,u))
    DVf .= 0.0
    grad!(DVf,u)
    product!(VDVf,Vtf,DVf)
    u .= 0.0
    grid_interpolate!(u,VDVf)
    u
end



force(f::VectorData{N},sys) where {N} = f*cellsize(sys)^2

function pressurejump(fds::VectorData{N},b::Union{Body,BodyList},sys::NavierStokes) where {N}
    @assert N == numpts(b)

    _hasfilter(sys) ? (fdsf = similar(fds); fdsf .= sys.Cmat*fds) : fdsf = deepcopy(fds)

    nrm = VectorData(normalmid(b))
    return (nrm.u∘fdsf.u + nrm.v∘fdsf.v)./dlengthmid(b) # This might need some modification
end
