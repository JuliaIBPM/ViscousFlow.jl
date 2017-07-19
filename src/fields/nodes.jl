struct Nodes{C <: CellType} <: AbstractMatrix{Float64}
    celldims::Tuple{Int, Int}
    data::Matrix{Float64}
    ghostlayers::Int
end

@wraparray Nodes data

function Nodes(::Type{Dual}, celldims::Tuple{Int, Int}, ghostlayers=0)
    data = zeros(celldims .+ 2ghostlayers)
    Nodes{Dual}(celldims, data, ghostlayers)
end

function Base.show(io::IO, nodes::Nodes{T}) where {T}
    celldims = "$(nodes.celldims[1])×$(nodes.celldims[2])"
    println(io, "$T nodes for a $celldims cell grid with $(edges.ghostlayers) ghost layers")
    print(io, "  Internal nodes: $celldims")
end
