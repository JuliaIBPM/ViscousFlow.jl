### Computing fields and other physical quantities ###

@inline vorticity(w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY}) where {NX,NY} =
        w/cellsize(sys)
        
#velocity(w::Nodes{Dual},sys) = -curl(sys.L\w)
function velocity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},Δus::VectorData{N},nrm::VectorData{N},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    u .= -curl(sys.L\w) + cellsize(sys)*grad(sys.Lc\sys.sc(pointwise_dot(Δus,nrm)))
    return u
end
function velocity!(u::Edges{Primal,NX,NY},w::Nodes{Dual,NX,NY},sys::NavierStokes{NX,NY,N}) where {NX,NY,N}
    u .= -curl(sys.L\w)
    return u
end
velocity(w::Nodes{Dual,NX,NY},a...) where {NX,NY} = velocity!(Edges(Primal,(NX,NY)),w,a...)


function streamfunction(w::Nodes{Dual},sys)
  ψ = -cellsize(sys)*(sys.L\w)

  xg, yg = coordinates(ψ,sys.grid)
  ψ .+= sys.U∞[1]*yg' .- sys.U∞[2]*xg
  return ψ
end

nl(w::Nodes{Dual},sys) = ConstrainedSystems.r₁(w,0.0,sys)/cellsize(sys)

force(f::VectorData{N},sys) where {N} = f*cellsize(sys)^2

function pressurejump(fds::VectorData{N},b::Union{Body,BodyList},sys::NavierStokes) where {N}
    @assert N == numpts(b)

    _hasfilter(sys) ? (fdsf = similar(fds); fdsf .= sys.Cmat*fds) : fdsf = deepcopy(fds)

    nrm = VectorData(normalmid(b))
    return (nrm.u∘fdsf.u + nrm.v∘fdsf.v)./dlengthmid(b) # This might need some modification
end
