struct Edges{C <: CellType}
    nodedims::Tuple{Int, Int}
    u::Matrix{Float64}
    v::Matrix{Float64}
end

edge_inds(::Type{Dual},   nodedims) = (nodedims[1]-1, nodedims[2]-2), (nodedims[1]-2, nodedims[2]-1)
edge_inds(::Type{Primal}, nodedims) = (nodedims[1], nodedims[2]-1), (nodedims[1]-1, nodedims[2])

function Edges(T::Type{C}, nodedims::Tuple{Int, Int}) where {C <: CellType}
    udims, vdims = edge_inds(T, nodedims)
    u = zeros(udims)
    v = zeros(vdims)
    Edges{T}(nodedims, u, v)
end

Edges(T, nodes::DualNodes) = Edges(T, size(nodes))

function shift!(dual::Edges{Dual}, primal::Edges{Primal})
    @assert dual.nodedims == primal.nodedims

    uₚ = primal.u
    for y in 1:size(dual.u,2), x in 1:size(dual.u,1)
        dual.u[x,y] = (uₚ[x,y] + uₚ[x+1,y] + uₚ[x,y+1] + uₚ[x+1,y+1])/4
    end

    vₚ = primal.v
    for y in 1:size(dual.v,2), x in 1:size(dual.v,1)
        dual.v[x,y] = (vₚ[x,y] + vₚ[x+1,y] + vₚ[x,y+1] + vₚ[x+1,y+1])/4
    end
    dual
end

shift(primal::Edges{Primal}) = shift!(Edges(Dual, primal.nodedims), primal)

function Base.show(io::IO, edges::Edges{T}) where {T}
    nodedims = "(nx = $(edges.nodedims[1]), ny = $(edges.nodedims[2]))"
    udims = "(nx = $(size(edges.u,1)), ny = $(size(edges.u,2)))"
    vdims = "(nx = $(size(edges.v,1)), ny = $(size(edges.v,2)))"
    println(io, "$T edges for a $nodedims cell grid")
    println(io, "  Internal u-faces: $udims")
    print(io, "  Internal v-faces: $vdims")
end
