struct DualNodes <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end

DualNodes(dims::Tuple{Int,Int}) = DualNodes(zeros(dims))
DualNodes(nx::Int, ny::Int) = DualNodes(zeros(nx, ny))

@wraparray DualNodes data
