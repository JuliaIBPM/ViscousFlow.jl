import Base: fill!, ∘

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
struct Edges{C <: CellType, NX, NY} <: AbstractMatrix{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
end

# Based on number of dual nodes, return the number of edges
edge_inds(::Type{Dual},   dualnodedims) = (dualnodedims[1]-1, dualnodedims[2]), (dualnodedims[1], dualnodedims[2]-1)
edge_inds(::Type{Primal}, dualnodedims) = (dualnodedims[1], dualnodedims[2]-1), (dualnodedims[1]-1, dualnodedims[2])

function Edges(T::Type{C}, dualnodedims::Tuple{Int, Int}) where {C <: CellType}
    udims, vdims = edge_inds(T, dualnodedims)
    u = zeros(udims)
    v = zeros(vdims)
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

"""
    cellshift!(v::Edges{Dual/Primal},q::Edges{Primal/Dual})

Shift (by linear interpolation) the primal (resp. dual) edge data `q` to the
edges of the dual (resp. primal) cells, and return the result in `v`.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> v = Edges(Dual,(8,6));

julia> Fields.cellshift!(v,q)
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.25  0.25  0.0  0.0  0.0  0.0
 0.0  0.25  0.25  0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function cellshift!(dual::Edges{Dual, NX, NY},
                primal::Edges{Primal, NX, NY}) where {NX, NY}
    uₚ = primal.u
    @inbounds for y in 2:NY-1, x in 1:NX-1
        dual.u[x,y] = (uₚ[x,y] + uₚ[x+1,y] + uₚ[x,y-1] + uₚ[x+1,y-1])/4
    end

    vₚ = primal.v
    @inbounds for y in 1:NY-1, x in 2:NX-1
        dual.v[x,y] = (vₚ[x,y] + vₚ[x-1,y] + vₚ[x,y+1] + vₚ[x-1,y+1])/4
    end
    dual
end

function cellshift(primal::Edges{Primal, NX, NY}) where {NX, NY}
    cellshift!(Edges(Dual, (NX, NY)), primal)
end

function cellshift!(primal::Edges{Primal, NX, NY},
                dual::Edges{Dual, NX, NY}) where {NX, NY}
    uₚ = dual.u
    @inbounds for y in 1:NY-1, x in 2:NX-1
        primal.u[x,y] = (uₚ[x,y] + uₚ[x-1,y] + uₚ[x,y+1] + uₚ[x-1,y+1])/4
    end

    vₚ = dual.v
    @inbounds for y in 2:NY-1, x in 1:NX-1
        primal.v[x,y] = (vₚ[x,y] + vₚ[x+1,y] + vₚ[x,y-1] + vₚ[x+1,y-1])/4
    end
    primal
end

function cellshift(dual::Edges{Dual, NX, NY}) where {NX, NY}
    cellshift!(Edges(Primal, (NX, NY)), dual)
end

"""
    product!(out::Edges/Nodes,p::Edges/Nodes,q::Edges/Nodes)

Compute the Hadamard (i.e. element by element) product of edge or nodal
(primal or dual) data `p` and `q` and return the result in `out`.

# Example

```jldoctest
julia> q = Edges(Dual,(8,6));

julia> out = p = deepcopy(q);

julia> q.u[3,2] = 0.3;

julia> p.u[3,2] = 0.2;

julia> product!(out,p,q)
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.06  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function product!(out::Edges{T, NX, NY},
                  p::Edges{T, NX, NY},
                  q::Edges{T, NX, NY}) where {T, NX, NY}

    uinds, vinds = edge_inds(T, (NX, NY))
    @inbounds for y in 1:uinds[2], x in 1:uinds[1]
        out.u[x,y] = p.u[x,y] * q.u[x,y]
    end

    @inbounds for y in 1:vinds[2], x in 1:vinds[1]
        out.v[x,y] = p.v[x,y] * q.v[x,y]
    end
    out
end

"""
    product(p::Edges/Nodes,q::Edges/Nodes) --> Edges/Nodes

Compute the Hadamard product of edge or nodal (primal or dual) data `p` and `q` and return
the result. This operation can also be carried out with the `∘` operator:

# Example

```jldoctest
julia> q = Edges(Dual,(8,6));

julia> p = deepcopy(q);

julia> q.u[3,2] = 0.3;

julia> p.u[3,2] = 0.2;

julia> p∘q
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
 0.0  0.0  0.06  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function product(p::Edges{T, NX, NY}, q::Edges{T, NX, NY}) where {T, NX, NY}
    product!(Edges(T, (NX, NY)), p, q)
end

function (∘)(p::Edges{T, NX, NY}, q::Edges) where {T, NX, NY}
    product!(Edges(T, (NX, NY)), p, q)
end

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
