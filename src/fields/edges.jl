# VECTOR EDGE DATA


"""
    Edges

`Edges` is a wrapper for vector-valued data that lie at the faces of either dual cells or
primary cells. `Edges` type data have fields `u` and `v` for the components of the
vector field. These are the normal components of the vector field on the vertical
and horizontal faces of the corresponding cell.

# Constructors
- `Edges(C,dims)` creates a vector field of zeros in cells of type `C` (where `C` is
  either `Dual` or `Primal`), on a grid of dimensions `dims`. Note that `dims`
  represent the number of dual cells on the grid.
- `Edges(C,w)` performs the same construction, but uses existing field data `w`
  of `GridData` type to determine the size of the grid.
-  Adding the `dtype=` keyword allows the data type of the field data to be
  changed. The default is `Float64`, but can be changed to, e.g., `ComplexF64`
"""
struct Edges{C <: CellType, NX, NY, T <: Number, DT} <: VectorGridData{NX,NY,T}
    data::DT
    u::XEdges{C,NX,NY,T}
    v::YEdges{C,NX,NY,T}
end

# Based on number of dual nodes, return the number of edges
edge_inds(::Type{C},   dualnodedims) where {C <: CellType} =
            xedge_inds(C,dualnodedims), yedge_inds(C,dualnodedims)


function (::Type{Edges{C,NX,NY,T,DT}})(data::AbstractVector{R}) where {C<: CellType,NX,NY,T<:Number,DT,R}
    udims, vdims = edge_inds(C, (NX, NY))
    n0 = 0
    dims = udims
    dn = prod(dims)
    u = reshape(view(data,n0+1:n0+dn),dims)
    n0 += dn
    dims = vdims
    dn = prod(dims)
    v = reshape(view(data,n0+1:n0+dn),dims)
    Edges{C, NX, NY,R,typeof(data)}(data, XEdges{C,NX,NY,R,typeof(u)}(u),
                                          YEdges{C,NX,NY,R,typeof(v)}(v))
end

function Edges(::Type{C}, dualnodedims::Tuple{Int, Int};dtype=Float64) where {C <: CellType}
    udims, vdims = edge_inds(C, dualnodedims)
    data = zeros(dtype,prod(udims)+prod(vdims))
    Edges{C,dualnodedims...,dtype,typeof(data)}(data)
end

@griddata(Edges,1)

#Base.IndexStyle(::Type{<:Edges}) = IndexLinear() # necessary?

function Base.show(io::IO, edges::Edges{C, NX, NY, T, DT}) where {C, NX, NY, T, DT}
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
