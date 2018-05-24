#=
Philosophy and convention:
A node is defined as the center of a grid cell. Grid cells may be those
in a primary grid or those in a dual grid.

The definitions "dual" and "primal" are made based on their conventional use
in holding fluid dynamic data. For example, pressure is held in primal nodes,
vorticity and streamfunction at dual nodes.

The definitions here are made with a "dual grid" in mind. That is, the grid is
defined by an integer number of dual cells in each direction. If a "primal grid"
is needed, then all of the defintions can be swapped (primal -> dual, dual -> primal).

Also, note that there might be dual cells that are "ghosts" (i.e. lie outside
the grid), but these are not distinguished in these basic definitions and operators.
=#


module Fields

import Base: @propagate_inbounds
export Primal, Dual, Edges, Nodes, DualNodes, othertype,
       curl, curl!, divergence, divergence!, gradient, gradient!,
       laplacian, laplacian!, Laplacian,
       product, product!, ∘,
       CircularConvolution

abstract type CellType end
abstract type Primal <: CellType end
abstract type Dual <: CellType end

macro wraparray(wrapper, field)
    T = supertype(eval(wrapper))
    @assert T <: AbstractArray "Wrapped type must be a subtype of AbstractArray"
    el_type, N = T.parameters

    quote
        Base.parent(A::$wrapper) = A.$field
        Base.size(A::$wrapper) = size(A.$field)
        Base.indices(A::$wrapper) = indices(A.$field)

        function Base.show(io::IO, ::MIME"text/plain", A::$wrapper)
          println(io,"$(typeof(A)) data")
          println(io,"Printing in grid orientation (lower left is (1,1)):")
          show(io,"text/plain",flipdim(transpose(A.$field),1))
        end

        @propagate_inbounds Base.getindex(A::$wrapper, i::Int) = A.$field[i]
        @propagate_inbounds Base.getindex(A::$wrapper, I::Vararg{Int, $N}) = A.$field[I...]
        @propagate_inbounds Base.setindex!(A::$wrapper, v, i::Int) = A.$field[i] = convert($el_type, v)
        @propagate_inbounds Base.setindex!(A::$wrapper, v, I::Vararg{Int, $N}) = A.$field[I...] = convert($el_type, v)
    end
end

function othertype end

macro othertype(celltype, k)
    esc(quote
        Fields.othertype(::$celltype) = $k
        Fields.othertype(::Type{$celltype}) = $k
    end)
end

@othertype Primal Dual
@othertype Dual Primal
@othertype CellType CellType

include("fields/nodes.jl")
include("fields/edges.jl")


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

include("fields/operators.jl")


function shift!(dual::Edges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-2, x in 1:NX-1
        dual.u[x,y] = (w[x,y+1] + w[x+1,y+1])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-2
        dual.v[x,y] = (w[x+1,y] + w[x+1,y+1])/2
    end
    dual
end

shift(nodes::Nodes{Dual,NX,NY}) where {NX,NY} = shift!(Edges(Dual, nodes), nodes)

function shift!(dual::Tuple{Nodes{Dual, NX, NY},Nodes{Dual, NX, NY}}, w::Edges{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 1:NX
        dual[1][x,y] = (w.u[x,y-1] + w.u[x,y])/2
    end

    @inbounds for y in 1:NY, x in 2:NX-1
        dual[2][x,y] = (w.v[x-1,y] + w.v[x,y])/2
    end
    dual
end

nodeshift(edges::Edges{Primal,NX,NY}) where {NX,NY} = shift!((Nodes(Dual, (NX,NY)),Nodes(Dual, (NX,NY))),edges)



#### to be removed

function shift!(dual::Edges{Dual, NX, NY}, w::DualNodes{NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-2, x in 1:NX-1
        dual.u[x,y] = (w[x,y+1] + w[x+1,y+1])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-2
        dual.v[x,y] = (w[x+1,y] + w[x+1,y+1])/2
    end
    dual
end

shift(nodes::DualNodes) = shift!(Edges(Dual, nodes), nodes)

end
