import Base: fill!

# EDGE COMPONENT DATA

struct XEdges{C <: CellType, NX, NY} <: ScalarGridData{NX,NY}
    data::Matrix{Float64}
end

struct YEdges{C <: CellType, NX, NY} <: ScalarGridData{NX,NY}
    data::Matrix{Float64}
end

# Number of indices
# Based on number of dual nodes, return the number of edges
xedge_inds(::Type{Dual}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2])
yedge_inds(::Type{Dual}, dualnodedims) = (dualnodedims[1], dualnodedims[2]-1)

xedge_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1], dualnodedims[2]-1)
yedge_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1]-1, dualnodedims[2])

# Constructors

function XEdges(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    dims = xedge_inds(T, dualnodedims)
    XEdges{T, dualnodedims...}(zeros(dims))
end

function YEdges(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    dims = yedge_inds(T, dualnodedims)
    YEdges{T, dualnodedims...}(zeros(dims))
end

XEdges(T, ::ScalarGridData{NX,NY}) where {NX, NY} = XEdges(T, (NX, NY) )
XEdges(T, ::VectorGridData{NX,NY}) where {NX, NY} = XEdges(T, (NX, NY) )

YEdges(T, ::ScalarGridData{NX,NY}) where {NX, NY} = YEdges(T, (NX, NY) )
YEdges(T, ::VectorGridData{NX,NY}) where {NX, NY} = YEdges(T, (NX, NY) )


XEdges(T, nx::Int, ny::Int) = XEdges(T,(nx,ny))
YEdges(T, nx::Int, ny::Int) = YEdges(T,(nx,ny))

(::Type{XEdges{T,NX,NY}})() where {T,NX,NY} = XEdges(T, (NX, NY))
(::Type{YEdges{T,NX,NY}})() where {T,NX,NY} = YEdges(T, (NX, NY))


Base.similar(::XEdges{T,NX,NY}) where {T,NX,NY} = XEdges(T, (NX, NY))
Base.similar(::YEdges{T,NX,NY}) where {T,NX,NY} = YEdges(T, (NX, NY))


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
struct Edges{C <: CellType, NX, NY} <: VectorGridData{NX,NY}
    u::XEdges{C,NX,NY}
    v::YEdges{C,NX,NY}
end

# Based on number of dual nodes, return the number of edges
edge_inds(T::Type{C},   dualnodedims) where {C <: CellType} =
            xedge_inds(T,dualnodedims), yedge_inds(T,dualnodedims)

function Edges(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    udims, vdims = edge_inds(T, dualnodedims)
    u = XEdges(T,dualnodedims...)
    v = YEdges(T,dualnodedims...)
    Edges{T, dualnodedims...}(u, v)
end

Edges(T, nodes::Nodes{S,NX,NY}) where {S <: CellType, NX,NY} = Edges(T, (NX, NY))
(::Type{Edges{T,NX,NY}})() where {T,NX,NY} = Edges(T, (NX, NY))

Base.similar(::Edges{T,NX,NY}) where {T,NX,NY} = Edges(T, (NX, NY))

function fill!(edges::Edges, s::Number)
    fill!(edges.u, s)
    fill!(edges.v, s)
    edges
end

Base.size(A::Edges{C,NX,NY}) where {C,NX,NY} = (length(A.u)+length(A.v),1)
@propagate_inbounds Base.getindex(A::Edges{C,NX,NY},i::Int) where {C,NX,NY} =
   i > length(A.u) ? A.v[i-length(A.u)] : A.u[i]
@propagate_inbounds Base.setindex!(A::Edges{C,NX,NY}, v, i::Int) where {C,NX,NY} =
   i > length(A.u) ? A.v[i-length(A.u)] = convert(Float64, v) : A.u[i] = convert(Float64, v)
Base.IndexStyle(::Type{<:Edges}) = IndexLinear()





function Base.show(io::IO, edges::Edges{T, NX, NY}) where {T, NX, NY}
    nodedims = "(nx = $NX, ny = $NY)"
    udims = "(nx = $(size(edges.u,1)), ny = $(size(edges.u,2)))"
    vdims = "(nx = $(size(edges.v,1)), ny = $(size(edges.v,2)))"
    println(io, "$T edges for a $nodedims cell grid")
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
