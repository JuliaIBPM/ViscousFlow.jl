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
struct Nodes{C <: CellType, NX, NY, T <: Number} <: ScalarGridData{NX,NY,T}
    data::Matrix{T}
end

# Number of indices
# Based on number of dual nodes, return the number of nodes
node_inds(::Type{Dual},   dualnodedims) = (dualnodedims[1], dualnodedims[2])
node_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2]-1)

# Constructors

function Nodes(T::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    dims = node_inds(T, dualnodedims)
    Nodes{T, dualnodedims...,dtype}(zeros(dtype,dims))
end

# This allows easy construction of nodes of either type from existing nodes of either
# type on the same grid.
Nodes(C, ::GridData{NX,NY,T};dtype=T) where {NX, NY,T <: Number} = Nodes(C, (NX, NY),dtype=dtype )


Nodes(C, nx::Int, ny::Int;dtype=Float64) = Nodes(C,(nx,ny),dtype=dtype)
(::Type{Nodes{C,NX,NY,T}})() where {C,NX,NY,T} = Nodes(C, (NX, NY),dtype=T)

Base.similar(::Nodes{C,NX,NY,T};element_type=T) where {C,NX,NY,T} = Nodes(C, (NX, NY),dtype=element_type)

function Base.show(io::IO, nodes::Nodes{C, NX, NY,T}) where {C, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    dims = "(nx = $(size(nodes,1)), ny = $(size(nodes,2)))"
    println(io, "$C nodes in a $nodedims cell grid of type $T data")
    print(io, "  Number of $C nodes: $dims")
end
