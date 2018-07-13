# Collections of data

struct EdgeGradient{C <: CellType,D <: CellType, NX,NY}
  dudx :: Nodes{C,NX,NY}
  dvdy :: Nodes{C,NX,NY}
  dudy :: Nodes{D,NX,NY}
  dvdx :: Nodes{D,NX,NY}
end

function EdgeGradient(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    dudxdims = node_inds(T, dualnodedims)
    dudydims = node_inds(othertype(C), dualnodedims)

    EdgeGradient{T, othertype(T), dualnodedims...}(
            Nodes(T,dualnodedims),Nodes(T,dualnodedims),
            Nodes(othertype(T),dualnodedims),Nodes(othertype(T),dualnodedims)
            )
end

function Base.show(io::IO, nodes::EdgeGradient{T, S, NX, NY}) where {T, S, NX, NY}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.dudx,1)), ny = $(size(nodes.dudx,2)))"
    vdims = "(nx = $(size(nodes.dudy,1)), ny = $(size(nodes.dudy,2)))"
    println(io, "Edge gradient data in a $nodedims cell grid")
    println(io, "  du/dx and dv/dy as $T nodes : $udims")
    println(io, "  du/dy and dv/dx as $S nodes : $vdims")
end

struct NodePair{C <: CellType,D <: CellType, NX,NY}
  u :: Nodes{C,NX,NY}
  v :: Nodes{D,NX,NY}
end

function NodePair(T::Type{C},S::Type{D},dualnodedims::Tuple{Int,Int}) where {C <: CellType, D <: CellType}
    udims = node_inds(T, dualnodedims)
    vdims = node_inds(S, dualnodedims)

    NodePair{T, S, dualnodedims...}(
            Nodes(T,dualnodedims),
            Nodes(S,dualnodedims)
            )
end
NodePair(T::Type{C}, U::Type{D},nodes::Nodes{S,NX,NY}) where
        {C <: CellType, D <: CellType, S <: CellType, NX,NY} = NodePair(T, U, (NX, NY))

NodePair(T::Type{C},dualnodedims::Tuple{Int,Int}) where {C <: CellType}= NodePair(T,T,dualnodedims)

NodePair(T::Type{C}, nodes::Nodes{S,NX,NY}) where {C <: CellType, S <: CellType, NX,NY} = NodePair(T, (NX, NY))



(::Type{NodePair{T,S,NX,NY}})() where {T<:CellType,S<:CellType,NX,NY} = NodePair(T, S, (NX, NY))


function Base.show(io::IO, nodes::NodePair{T, S, NX, NY}) where {T, S, NX, NY}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(nodes.u,1)), ny = $(size(nodes.u,2)))"
    vdims = "(nx = $(size(nodes.v,1)), ny = $(size(nodes.v,2)))"
    println(io, "($T,$S) node data pair in a $nodedims cell grid")
    println(io, "  Number of $T nodes: $udims")
    println(io, "  Number of $S nodes: $vdims")
end

function fill!(nodes::NodePair, s::Number)
    fill!(nodes.u, s)
    fill!(nodes.v, s)
    nodes
end

Base.size(A::NodePair{C,D,NX,NY}) where {C,D,NX,NY} = (length(A.u)+length(A.v),1)
@propagate_inbounds Base.getindex(A::NodePair{C,D,NX,NY},i::Int) where {C,D,NX,NY} =
   i > length(A.u) ? A.v[i-length(A.u)] : A.u[i]
@propagate_inbounds Base.setindex!(A::NodePair{C,D,NX,NY}, v, i::Int) where {C,D,NX,NY} =
   i > length(A.u) ? A.v[i-length(A.u)] = convert(Float64, v) : A.u[i] = convert(Float64, v)
