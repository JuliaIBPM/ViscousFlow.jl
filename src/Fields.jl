module Fields

import Base: @propagate_inbounds
export Primal, Dual, Edges, Nodes,
       shift, shift!

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

include("fields/edges.jl")
include("fields/nodes.jl")

function shift!(dual::Edges{Dual}, nodes::Nodes{Dual})
    ω = nodes.data
    for j in 2:size(dual.u,2)-1, i in 1:size(dual.u,1)
        dual.u[i,j] = (ω[i,j] + ω[i+1,j])/2
    end

    for j in 1:size(dual.v,2), i in 2:size(dual.v,1)-1
        dual.v[i,j] = (ω[i,j] + ω[i,j+1])/2
    end
    dual
end

function shift(nodes::Nodes{Dual})
    dual = Edges(Dual, nodes.celldims, nodes.ghostlayers)
    shift!(dual, nodes)
end

end
