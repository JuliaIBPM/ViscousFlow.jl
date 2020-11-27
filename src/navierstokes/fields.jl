### Computing fields and other physical quantities ###
using LinearAlgebra

@inline vorticity(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY} = w/cellsize(sys)

function velocity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N},t::Real) where {NX,NY,N}
    u .= 0.0
    _velocity_vorticity!(u,w,sys)
    _velocity_single_layer!(u,sys,t)
    _velocity_freestream!(u,sys,t)
    return u
end

function _velocity_vorticity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    _unscaled_streamfunction_vorticity!(sys.Sn,w,sys)
    curl!(sys.Vf,sys.Sn)
    u .+= sys.Vf
    return u
end

function _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,ExternalInternalFlow},t::Real) where {NX,NY,N,MT,FS}
    u .+= 0.0
    return u
end

function _velocity_single_layer!(u::Edges{Primal,NX,NY},sys::NavierStokes{NX,NY,N,MT,FS,SD},t::Real) where {NX,NY,N,MT,FS,SD}
    surface_velocity_jump!(sys.Δus,sys,t)
    pointwise_dot!(sys.Sb,sys.Δus,normals(sys))
    sys.Sb .*= cellsize(sys)
    sys.slc(sys.Sc,sys.Sb)
    _unscaled_scalarpotential!(sys.Sc,sys.Sc,sys)
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


function _unscaled_scalarpotential!(ϕ::Nodes{Primal,NX,NY},dil::Nodes{Primal,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY}
  ldiv!(ϕ,sys.L,dil)
  return ϕ
end



nl(w::Nodes{Dual},sys) = ConstrainedSystems.r₁(w,0.0,sys)/cellsize(sys)

force(f::VectorData{N},sys) where {N} = f*cellsize(sys)^2

function pressurejump(fds::VectorData{N},b::Union{Body,BodyList},sys::NavierStokes) where {N}
    @assert N == numpts(b)

    _hasfilter(sys) ? (fdsf = similar(fds); fdsf .= sys.Cmat*fds) : fdsf = deepcopy(fds)

    nrm = VectorData(normalmid(b))
    return (nrm.u∘fdsf.u + nrm.v∘fdsf.v)./dlengthmid(b) # This might need some modification
end
