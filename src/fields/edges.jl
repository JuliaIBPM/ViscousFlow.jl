import Base: fill!

# EDGE COMPONENT DATA

struct XEdges{C <: CellType, NX, NY, T <: Number} <: ScalarGridData{NX,NY,T}
    data::Matrix{T}
end

struct YEdges{C <: CellType, NX, NY, T <: Number} <: ScalarGridData{NX,NY,T}
    data::Matrix{T}
end

# Number of indices
# Based on number of dual nodes, return the number of edges
xedge_inds(::Type{Dual}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2])
yedge_inds(::Type{Dual}, dualnodedims) = (dualnodedims[1], dualnodedims[2]-1)

xedge_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1], dualnodedims[2]-1)
yedge_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2])

# Constructors

function XEdges(T::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    dims = xedge_inds(T, dualnodedims)
    XEdges{T, dualnodedims...,dtype}(zeros(dtype,dims))
end

function YEdges(T::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    dims = yedge_inds(T, dualnodedims)
    YEdges{T, dualnodedims...,dtype}(zeros(dtype,dims))
end

XEdges(C, ::GridData{NX,NY,T}) where {NX, NY, T <: Number} = XEdges(C, (NX, NY), dtype = T )

YEdges(C, ::GridData{NX,NY,T}) where {NX, NY, T <: Number} = YEdges(C, (NX, NY), dtype = T )


XEdges(C, nx::Int, ny::Int;dtype=Float64) = XEdges(C,(nx,ny),dtype=dtype)
YEdges(C, nx::Int, ny::Int;dtype=Float64) = YEdges(C,(nx,ny),dtype=dtype)

(::Type{XEdges{C,NX,NY,T}})() where {C,NX,NY,T} = XEdges(C, (NX, NY),dtype=T)
(::Type{YEdges{C,NX,NY,T}})() where {C,NX,NY,T} = YEdges(C, (NX, NY),dtype=T)


Base.similar(::XEdges{C,NX,NY,T}) where {C,NX,NY,T} = XEdges(C, (NX, NY),dtype=T)
Base.similar(::YEdges{C,NX,NY,T}) where {C,NX,NY,T} = YEdges(C, (NX, NY),dtype=T)

function Base.show(io::IO, xedges::XEdges{C, NX, NY, T}) where {C, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    dims = "(nx = $(size(xedges,1)), ny = $(size(xedges,2)))"
    println(io, "$C x-edges in a $nodedims cell grid of type $T data")
    print(io, "  Number of $C nodes: $dims")
end

function Base.show(io::IO, yedges::YEdges{C, NX, NY, T}) where {C, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    dims = "(nx = $(size(yedges,1)), ny = $(size(yedges,2)))"
    println(io, "$C y-edges in a $nodedims cell grid of type $T data")
    print(io, "  Number of $C nodes: $dims")
end

# VECTOR EDGE DATA


"""
    Edges{Dual/Primal}

`Edges` is a wrapper for vector-valued data that lie at the faces of either dual cells or
primary cells. `Edges` type data have fields `u` and `v` for the components of the
vector field. These are the normal components of the vector field on the vertical
and horizontal faces of the corresponding cell.

# Constructors
- `Edges(C,dims)` creates a vector field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid.
- `Edges(C,w)` performs the same construction, but uses existing field data `w`
  of `Nodes` type to determine the size of the grid.
"""
struct Edges{C <: CellType, NX, NY, T <: Number} <: VectorGridData{NX,NY,T}
    u::XEdges{C,NX,NY,T}
    v::YEdges{C,NX,NY,T}
end

# Based on number of dual nodes, return the number of edges
edge_inds(T::Type{C},   dualnodedims) where {C <: CellType} =
            xedge_inds(T,dualnodedims), yedge_inds(T,dualnodedims)

function Edges(T::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    udims, vdims = edge_inds(T, dualnodedims)
    u = XEdges(T,dualnodedims...,dtype=dtype)
    v = YEdges(T,dualnodedims...,dtype=dtype)
    Edges{T, dualnodedims...,dtype}(u, v)
end

(::Type{Edges{C,NX,NY,T}})() where {C,NX,NY,T} = Edges(C, (NX, NY),dtype=T)

Edges(C, ::GridData{NX,NY,T}) where {NX, NY,T} = Edges(C, (NX,NY),dtype=T)

Base.similar(::Edges{C,NX,NY,T}) where {C,NX,NY,T} = Edges(C, (NX, NY),dtype=T)

function fill!(edges::Edges, s::Number)
    fill!(edges.u, s)
    fill!(edges.v, s)
    edges
end

Base.size(A::Edges{C,NX,NY}) where {C,NX,NY} = (length(A.u)+length(A.v),1)
@propagate_inbounds Base.getindex(A::Edges{C,NX,NY,T},i::Int) where {C,NX,NY,T <: Number} =
   i > length(A.u) ? A.v[i-length(A.u)] : A.u[i]
@propagate_inbounds Base.setindex!(A::Edges{C,NX,NY,T}, v, i::Int) where {C,NX,NY,T <: Number} =
   i > length(A.u) ? A.v[i-length(A.u)] = convert(T, v) : A.u[i] = convert(T, v)
Base.IndexStyle(::Type{<:Edges}) = IndexLinear()





function Base.show(io::IO, edges::Edges{C, NX, NY, T}) where {C, NX, NY, T}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(edges.u,1)), ny = $(size(edges.u,2)))"
    vdims = "(nx = $(size(edges.v,1)), ny = $(size(edges.v,2)))"
    println(io, "$C edges for a $nodedims cell grid of type $T data")
    println(io, "  Internal u-faces: $udims")
    print(io, "  Internal v-faces: $vdims")
end

function Base.show(io::IO, m::MIME"text/plain", edges::Edges)
    println(io,"$(typeof(edges)) data")
    println(io,"u (in grid orientation)")
    show(io,m,reverse(transpose(edges.u),dims=1))
    println(io)
    println(io,"v (in grid orientation)")
    show(io,m,reverse(transpose(edges.v),dims=1))
end

# function Base.show(io::IO, ::MIME"text/plain", edges::Edges{T, NX, NY}) where {T, NX, NY}
#     println(io,"$(typeof(edges)) data")
#     println(io,"u (in grid orientation):")
#     #Base.showarray(io,flipdim(transpose(edges.u),1),false;header=false)
#     println(io)
#     println(io,"v (in grid orientation):")
#     #Base.showarray(io,flipdim(transpose(edges.v),1),false;header=false)
# end
