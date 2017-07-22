struct DualNodes{NX, NY} <: AbstractMatrix{Float64}
    data::Matrix{Float64}
end

DualNodes(dims::Tuple{Int,Int}) = DualNodes{dims...}(zeros(dims))
DualNodes(nx::Int, ny::Int) = DualNodes{nx, ny}(zeros(nx, ny))
(::Type{DualNodes{NX, NY}})() where {NX, NY} = DualNodes(NX, NY)

@wraparray DualNodes data
