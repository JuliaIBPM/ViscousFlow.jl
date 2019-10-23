# Collections of data

struct EdgeGradient{C <: CellType,D <: CellType, NX,NY, T<: Number} <: GridData{NX,NY,T}
  dudx :: Nodes{C,NX,NY,T}
  dvdy :: Nodes{C,NX,NY,T}
  dudy :: Nodes{D,NX,NY,T}
  dvdx :: Nodes{D,NX,NY,T}
end

function EdgeGradient(T::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    dudxdims = node_inds(T, dualnodedims)
    dudydims = node_inds(othertype(C), dualnodedims)

    EdgeGradient{T, othertype(T), dualnodedims...,dtype}(
            Nodes(T,dualnodedims,dtype=dtype),Nodes(T,dualnodedims,dtype=dtype),
            Nodes(othertype(T),dualnodedims,dtype=dtype),Nodes(othertype(T),dualnodedims,dtype=dtype)
            )
end

(::Type{EdgeGradient{R,S,NX,NY,T}})() where {R<:CellType,S<:CellType,NX,NY,T} = EdgeGradient(R, (NX, NY),dtype=T)

EdgeGradient(C, ::GridData{NX,NY,T}) where {NX, NY,T} = EdgeGradient(C, (NX,NY),dtype=T)

function Base.fill!(g::EdgeGradient, s::Number)
    fill!(g.dudx, s)
    fill!(g.dudy, s)
    fill!(g.dvdx, s)
    fill!(g.dvdy, s)
    g
end

Base.size(A::EdgeGradient{C,D,NX,NY,T}) where {C,D,NX,NY,T} = (2length(A.dudx)+2length(A.dudy),1)
@propagate_inbounds Base.getindex(A::EdgeGradient{C,D,NX,NY,T},i::Int) where {C,D,NX,NY,T} =
   i > length(A.dudx) ? (i > length(A.dudx)+length(A.dudy) ? (i > length(A.dudx)+2length(A.dudy) ? A.dvdy[i-length(A.dudx)-2length(A.dudy)] : A.dvdx[i-length(A.dudx)-length(A.dudy)]) : A.dudy[i-length(A.dudx)]) : A.dudx[i]
@propagate_inbounds Base.setindex!(A::EdgeGradient{C,D,NX,NY,T}, v, i::Int) where {C,D,NX,NY,T} =
   i > length(A.dudx) ? (i > length(A.dudx)+length(A.dudy) ? (i > length(A.dudx)+2length(A.dudy) ? A.dvdy[i-length(A.dudx)-2length(A.dudy)] = convert(T, v) : A.dvdx[i-length(A.dudx)-length(A.dudy)] = convert(T, v)) : A.dudy[i-length(A.dudx)]= convert(T, v)) : A.dudx[i] = convert(T, v)
Base.IndexStyle(::Type{<:EdgeGradient}) = IndexLinear()

function Base.show(io::IO, nodes::EdgeGradient{R, S, NX, NY, T}) where {R, S, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.dudx,1)), ny = $(size(nodes.dudx,2)))"
    vdims = "(nx = $(size(nodes.dudy,1)), ny = $(size(nodes.dudy,2)))"
    println(io, "Edge gradient data in a $nodedims cell grid of type $T data")
    println(io, "  du/dx and dv/dy as $R nodes : $udims")
    println(io, "  du/dy and dv/dx as $S nodes : $vdims")
end

struct NodePair{C <: CellType,D <: CellType, NX,NY,T} <: GridData{NX,NY,T}
  u :: Nodes{C,NX,NY,T}
  v :: Nodes{D,NX,NY,T}
end

function NodePair(R::Type{C},S::Type{D},dualnodedims::Tuple{Int,Int};dtype=Float64) where {C <: CellType, D <: CellType}
    udims = node_inds(R, dualnodedims)
    vdims = node_inds(S, dualnodedims)

    NodePair{R, S, dualnodedims...,dtype}(
            Nodes(R,dualnodedims,dtype=dtype),
            Nodes(S,dualnodedims,dtype=dtype)
            )
end
NodePair(R::Type{C}, U::Type{D},nodes::Nodes{S,NX,NY,T}) where
        {C <: CellType, D <: CellType, S <: CellType, NX,NY,T} = NodePair(R, U, (NX, NY),dtype=T)

NodePair(R::Type{C},dualnodedims::Tuple{Int,Int};dtype=Float64) where {C <: CellType}= NodePair(R,R,dualnodedims,dtype=dtype)

NodePair(R::Type{C}, nodes::Nodes{S,NX,NY,T}) where {C <: CellType, S <: CellType, NX,NY,T} = NodePair(R, (NX, NY),dtype=T)


(::Type{NodePair{R,S,NX,NY,T}})() where {R<:CellType,S<:CellType,NX,NY,T} = NodePair(R, S, (NX, NY),dtype=T)

NodePair(R, S, ::GridData{NX,NY,T}) where {NX, NY,T} = NodePair(R, S, (NX,NY),dtype=T)

function Base.show(io::IO, nodes::NodePair{R, S, NX, NY, T}) where {R, S, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.u,1)), ny = $(size(nodes.u,2)))"
    vdims = "(nx = $(size(nodes.v,1)), ny = $(size(nodes.v,2)))"
    println(io, "($R,$S) node data pair in a $nodedims cell grid of type $T data")
    println(io, "  Number of $R nodes: $udims")
    println(io, "  Number of $S nodes: $vdims")
end

function fill!(nodes::NodePair, s::Number)
    fill!(nodes.u, s)
    fill!(nodes.v, s)
    nodes
end

Base.size(A::NodePair{C,D,NX,NY,T}) where {C,D,NX,NY,T} = (length(A.u)+length(A.v),1)
@propagate_inbounds Base.getindex(A::NodePair{C,D,NX,NY,T},i::Int) where {C,D,NX,NY,T} =
   i > length(A.u) ? A.v[i-length(A.u)] : A.u[i]
@propagate_inbounds Base.setindex!(A::NodePair{C,D,NX,NY,T}, v, i::Int) where {C,D,NX,NY,T} =
   i > length(A.u) ? A.v[i-length(A.u)] = convert(T, v) : A.u[i] = convert(T, v)
Base.IndexStyle(::Type{<:NodePair}) = IndexLinear()
