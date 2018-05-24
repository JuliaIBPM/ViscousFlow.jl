import Base: fill!, ∘

struct Edges{C <: CellType, NX, NY}
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

Edges(T, nodes::Nodes{Dual,NX,NY}) where {NX,NY} = Edges(T, size(nodes))
(::Type{Edges{T,NX,NY}})() where {T,NX,NY} = Edges(T, (NX, NY))

function fill!(edges::Edges, s::Number)
    fill!(edges.u, s)
    fill!(edges.v, s)
    edges
end

Edges(T, nodes::DualNodes) = Edges(T, size(nodes))


function shift!(dual::Edges{Dual, NX, NY},
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

function shift(primal::Edges{Primal, NX, NY}) where {NX, NY}
    shift!(Edges(Dual, (NX, NY)), primal)
end

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

function Base.show(io::IO, ::MIME"text/plain", edges::Edges{T, NX, NY}) where {T, NX, NY}
    println(io,"$T edge data")
    println(io,"u (in grid orientation):")
    show(io,"text/plain",flipdim(transpose(edges.u),1))
    println(io)
    println(io,"v (in grid orientation):")
    show(io,"text/plain",flipdim(transpose(edges.v),1))
end
