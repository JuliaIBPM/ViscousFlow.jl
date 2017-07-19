struct Nodes{C <: CellType} <: AbstractMatrix{Float64}
    celldims::Tuple{Int, Int}
    data::Matrix{Float64}
    ghostlayers::Int
end

@wraparray Nodes data

function Nodes(::Type{C}, celldims::Tuple{Int, Int}, ghostlayers=0) where {C <: CellType}
    data = zeros(celldims .+ 2ghostlayers)
    Nodes{C}(celldims, data, ghostlayers)
end

function Base.show(io::IO, nodes::Nodes{T}) where {T}
    celldims = "(nx = $(nodes.celldims[1]), ny = $(nodes.celldims[2]))"
    println(io, "$T nodes for a $celldims cell grid with $(nodes.ghostlayers) ghost layers")
    print(io, "  Internal nodes: $celldims")
end
