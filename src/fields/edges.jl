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
    for j in 2:size(dual.u,2)-1, i in 1:size(dual.u,1)
        dual.u[i,j] = (uₚ[i,j-1] + uₚ[i+1,j-1] + uₚ[i,j] + uₚ[i+1,j])/4
    end

    vₚ = primal.v
    for j in 1:size(dual.v,2), i in 2:size(dual.v,1)-1
        dual.v[i,j] = (vₚ[i-1,j] + vₚ[i,j] + vₚ[i-1,j+1] + vₚ[i,j+1])/4
    end
    dual
end

function shift(primal::Edges{Primal})
    dual = Edges(Dual, primal.celldims, primal.ghostlayers)
    shift!(dual, primal)
end

function Base.show(io::IO, edges::Edges{T}) where {T}
    celldims = "$(edges.celldims[1])×$(edges.celldims[2])"
    udims = "$(edges.celldims[1])×$(edges.celldims[2]-1)"
    vdims = "$(edges.celldims[1]-1)×$(edges.celldims[2])"
    println(io, "$T edges for a $celldims cell grid with $(edges.ghostlayers) ghost layers")
    println(io, "  Internal u-faces: $udims")
    print(io, "  Internal v-faces: $vdims")
end
