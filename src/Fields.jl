module Fields

import Base: @propagate_inbounds
export Primal, Dual, Edges, DualNodes,
       curl, curl!

abstract type CellType end
abstract type Primal <: CellType end
abstract type Dual <: CellType end

macro wraparray(wrapper, field)
    T = supertype(eval(wrapper))
    @assert T <: AbstractArray "Wrapped type must be a subtype of AbstractArray"
    el_type, N = T.parameters

    quote
        Base.parent(A::$wrapper) = A.$field
        Base.size(A::$wrapper) = size(A.$field)
        Base.indices(A::$wrapper) = indices(A.$field)

        @propagate_inbounds Base.getindex(A::$wrapper, i::Int) = A.$field[i]
        @propagate_inbounds Base.getindex(A::$wrapper, I::Vararg{Int, $N}) = A.$field[I...]
        @propagate_inbounds Base.setindex!(A::$wrapper, v, i::Int) = A.$field[i] = convert($el_type, v)
        @propagate_inbounds Base.setindex!(A::$wrapper, v, I::Vararg{Int, $N}) = A.$field[I...] = convert($el_type, v)
    end
end

include("fields/nodes.jl")
include("fields/edges.jl")
include("fields/operators.jl")

function shift!(dual::Edges{Dual}, nodes::DualNodes)
    @assert dual.nodedims == size(nodes)

    w = nodes.data
    for y in 1:size(dual.u,2), x in 1:size(dual.u,1)
        dual.u[x,y] = (w[x,y+1] + w[x+1,y+1])/2
    end

    for y in 1:size(dual.v,2), x in 1:size(dual.v,1)
        dual.v[x,y] = (w[x+1,y] + w[x+1,y+1])/2
    end
    dual
end

shift(nodes::DualNodes) = shift!(Edges(Dual, nodes), nodes)

end
