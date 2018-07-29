import Base: size, ∘

"""
    Nodes{Dual/Primal}

`Nodes` is a wrapper for scalar-valued data that lie at the centers of either dual cells or
primary cells. A `Nodes` type can be accessed by indexing like any other array,
and allows the use of [size].

# Constructors
- `Nodes(C,dims)` creates a field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid, even if `C` is `Primal`.
- `Nodes(C,w)` performs the same construction, but uses existing field data `w`
  of `Nodes` type to determine the size of the grid.
"""
struct Nodes{C <: CellType, NX, NY} <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end

# This macro allows us to access Nodes.data via just the wrapper itself
@wraparray Nodes data

# Based on number of dual nodes, return the number of nodes
node_inds(::Type{Dual},   dualnodedims) = (dualnodedims[1], dualnodedims[2])
node_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2]-1)

function Nodes(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    dims = node_inds(T, dualnodedims)
    Nodes{T, dualnodedims...}(zeros(dims))
end

# This allows easy construction of nodes of either type from existing nodes of either
# type on the same grid.
Nodes(T, nodes::Nodes{S,NX,NY}) where {S <: CellType, NX, NY} = Nodes(T, (NX, NY) )

Nodes(T, nx::Int, ny::Int) = Nodes(T,(nx,ny))
(::Type{Nodes{T,NX,NY}})() where {T,NX,NY} = Nodes(T, (NX, NY))

Base.similar(::Nodes{T,NX,NY}) where {T,NX,NY} = Nodes(T, (NX, NY))

function Base.show(io::IO, nodes::Nodes{T, NX, NY}) where {T, NX, NY}
    nodedims = "(nx = $NX, ny = $NY)"
    dims = "(nx = $(size(nodes,1)), ny = $(size(nodes,2)))"
    println(io, "$T nodes in a $nodedims cell grid")
    print(io, "  Number of $T nodes: $dims")
end

function product!(out::Nodes{T, NX, NY},
                  p::Nodes{T, NX, NY},
                  q::Nodes{T, NX, NY}) where {T, NX, NY}

    inds = node_inds(T, (NX, NY))
    @inbounds for y in 1:inds[2], x in 1:inds[1]
        out[x,y] = p[x,y] * q[x,y]
    end
    out
end

function product(p::Nodes{T, NX, NY}, q::Nodes{T, NX, NY}) where {T, NX, NY}
    product!(Nodes(T, (NX, NY)), p, q)
end

function (∘)(p::Nodes{T, NX, NY}, q::Nodes) where {T, NX, NY}
    product!(Nodes(T, (NX, NY)), p, q)
end
