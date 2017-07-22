import Base: ∘

struct Edges{C <: CellType, NX, NY}
    u::Matrix{Float64}
    v::Matrix{Float64}
end

edge_inds(::Type{Dual},   nodedims) = (nodedims[1]-1, nodedims[2]-2), (nodedims[1]-2, nodedims[2]-1)
edge_inds(::Type{Primal}, nodedims) = (nodedims[1], nodedims[2]-1), (nodedims[1]-1, nodedims[2])

function Edges(T::Type{C}, nodedims::Tuple{Int, Int}) where {C <: CellType}
    udims, vdims = edge_inds(T, nodedims)
    u = zeros(udims)
    v = zeros(vdims)
    Edges{T, nodedims...}(u, v)
end

Edges(T, nodes::DualNodes) = Edges(T, size(nodes))
(::Type{Edges{T,NX,NY}})() where {T,NX,NY} = Edges(T, (NX, NY))

function shift!(dual::Edges{Dual, NX, NY},
                primal::Edges{Primal, NX, NY}) where {NX, NY}
    uₚ = primal.u
    @inbounds for y in 1:NY-2, x in 1:NX-1
        dual.u[x,y] = (uₚ[x,y] + uₚ[x+1,y] + uₚ[x,y+1] + uₚ[x+1,y+1])/4
    end

    vₚ = primal.v
    @inbounds for y in 1:NY-1, x in 1:NX-2
        dual.v[x,y] = (vₚ[x,y] + vₚ[x+1,y] + vₚ[x,y+1] + vₚ[x+1,y+1])/4
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
