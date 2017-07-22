struct FieldPool{NX, NY}
    nodes::Vector{DualNodes{NX, NY}}
    pedges::Vector{Edges{Primal, NX, NY}}
    dedges::Vector{Edges{Dual, NX, NY}}
end

struct PoolIterator{NX, NY, T}
    pool::FieldPool
end

@generated function (itr::PoolIterator{NX, NY, T})(n::Int) where {NX, NY, T}
    if T <: DualNodes
        field = :nodes
    elseif T <: Edges{Dual}
        field = :dedges
    else
        field = :pedges
    end

    quote
        if length(itr.pool.$field) < n
            for i in 1:(n - length(itr.pool.$field))
                push!(itr.pool.$field, T{NX, NY}())
            end
        end
        fill!(itr.pool.$field[n], 0)
    end
end

function FieldPool(dims::Tuple{Int, Int})
    NX, NY = dims
    nodes = Vector{DualNodes{NX, NY}}()
    pedges = Vector{Edges{Primal, NX, NY}}()
    dedges = Vector{Edges{Dual, NX, NY}}()

    FieldPool{NX, NY}(nodes, pedges, dedges)
end
FieldPool(nx::Int, ny::Int) = FieldPool((nx, ny))

function (f::FieldPool{NX, NY})(::Type{T}) where {NX, NY, T}
    return PoolIterator{NX, NY, T}(f)
end

Base.start(pitr::PoolIterator) = 1
Base.next(pitr::PoolIterator, n::Int) = (pitr(n), n+1)
Base.done(pitr::PoolIterator, n::Int) = false
Base.iteratorsize(::Type{PoolIterator}) = IsInfinite()
Base.eltype(::Type{PoolIterator{NX, NY, T}}) where {NX, NY, T} = T{NX, NY}
