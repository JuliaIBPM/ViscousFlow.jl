struct Nodes{C <: CellType, NX, NY}
    data::Matrix{Float64}
end

# Based on number of dual nodes, return the number of edges
node_inds(::Type{Dual},   dualnodedims) = (dualnodedims[1], dualnodedims[2])
node_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2]-1)

function Nodes(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    dims = node_inds(T, dualnodedims)
    Nodes{T, dualnodedims...}(zeros(dims))
end

Nodes(T, nx::Int, ny::Int) = Nodes{T,nx,ny}(zeros(nx,ny))
(::Type{Nodes{T,NX,NY}})() where {T,NX,NY} = Nodes(T, (NX, NY))

function Base.show(io::IO, nodes::Nodes{T, NX, NY}) where {T, NX, NY}
    nodedims = "(nx = $NX, ny = $NY)"
    dims = "(nx = $(size(nodes.data,1)), ny = $(size(nodes.data,2)))"
    println(io, "$T nodes for a $nodedims cell grid")
    println(io, "  Number of $T nodes: $dims")
end

#### Stuff that will be removed later

struct DualNodes{NX, NY} <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end

DualNodes(dims::Tuple{Int,Int}) = DualNodes{dims...}(zeros(dims))
DualNodes(nx::Int, ny::Int) = DualNodes{nx, ny}(zeros(nx, ny))
(::Type{DualNodes{NX, NY}})() where {NX, NY} = DualNodes(NX, NY)

@wraparray DualNodes data
