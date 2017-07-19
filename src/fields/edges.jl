struct Edges{C <: CellType}
    celldims::Tuple{Int, Int}
    u::Matrix{Float64}
    v::Matrix{Float64}
    ghostlayers::Int
end

function Edges(::Type{Primal}, celldims::Tuple{Int, Int}, ghostlayers=0)
    u = zeros((celldims[1], celldims[2]-1) .+ 2ghostlayers)
    v = zeros((celldims[1]-1, celldims[2]) .+ 2ghostlayers)
    Edges{Primal}(celldims, u, v, ghostlayers)
end

function Edges(::Type{Dual}, celldims::Tuple{Int, Int}, ghostlayers=0)
    u = zeros((celldims[1]-1, celldims[2]) .+ 2ghostlayers)
    v = zeros((celldims[1], celldims[2]-1) .+ 2ghostlayers)
    Edges{Dual}(celldims, u, v, ghostlayers)
end

function shift!(dual::Edges{Dual}, primal::Edges{Primal})
    @assert dual.celldims == primal.celldims

    uₚ = primal.u
    for y in 2:size(dual.u,2)-1, x in 1:size(dual.u,1)
        dual.u[x,y] = (uₚ[x,y-1] + uₚ[x+1,y-1] + uₚ[x,y] + uₚ[x+1,y])/4
    end

    vₚ = primal.v
    for y in 1:size(dual.v,2), x in 2:size(dual.v,1)-1
        dual.v[x,y] = (vₚ[x-1,y] + vₚ[x,y] + vₚ[x-1,y+1] + vₚ[x,y+1])/4
    end
    dual
end

function shift(primal::Edges{Primal})
    dual = Edges(Dual, primal.celldims, primal.ghostlayers)
    shift!(dual, primal)
end

function Base.show(io::IO, edges::Edges{T}) where {T}
    celldims = "(nx = $(edges.celldims[1]), ny = $(edges.celldims[2]))"
    udims = "(nx = $(edges.celldims[1]), ny = $(edges.celldims[2]-1))"
    vdims = "(nx = $(edges.celldims[1]-1), ny = $(edges.celldims[2]))"
    println(io, "$T edges for a $celldims cell grid with $(edges.ghostlayers) ghost layers")
    println(io, "  Internal u-faces: $udims")
    print(io, "  Internal v-faces: $vdims")
end
