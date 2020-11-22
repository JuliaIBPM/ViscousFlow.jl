### Computing fields and other physical quantities ###

vorticity(w::Nodes{Dual},sys) = w/cellsize(sys)
velocity(w::Nodes{Dual},sys) = -curl(sys.L\w)

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
