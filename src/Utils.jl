module Utils

export @submodule, @get, MappedVector

macro submodule(mod)
    path = mod * ".jl"
    name = split(mod, '/')[end]
    Expr(:block,
         Expr(:call, :include, path),
         Expr(:using, :., Symbol(name)))
end

"""
A macro for extracting fields from an object.  For example, instead of a statement
like

    (obj.a + obj.b)^obj.c

it is more readable if we first locally bind the variables

    a = obj.a
    b = obj.b
    c = obj.c

    (a + b)^c

Using the `@get` marco, this becomes

    @get obj (a, b, c)

    (a + b)^c

We can also locally assign different names to the binding

    @get obj (a, b, c) (α, β, γ)
    (α + β)^γ
"""
macro get(object, fields...)
    if length(fields) == 1
        try
            @assert typeof(fields[1]) == Expr
            @assert fields[1].head == :tuple
            @assert all([typeof(arg) == Symbol for arg in fields[1].args])
        catch
            throw(ArgumentError("second argument must be a tuple of field names"))
        end
        esc(Expr(:block, [:($sym = $object.$sym) for sym in fields[1].args]...))
    elseif length(fields) == 2
        try
            @assert typeof(fields[1]) == Expr
            @assert typeof(fields[2]) == Expr
            @assert (fields[1].head == :tuple && fields[2].head == :tuple)
            @assert all([typeof(arg) == Symbol for arg in fields[1].args])
            @assert all([typeof(arg) == Symbol for arg in fields[2].args])
        catch
            throw(ArgumentError("second and third argument must be tuples of field names"))
        end

        nargs = length(fields[1].args)
        @assert nargs == length(fields[2].args) "field name tuples must have the same length"
        esc(Expr(:block, [:($(fields[2].args[i]) = $object.$(fields[1].args[i])) for i in 1:nargs]...))
    else
        throw(ArgumentError("Usage: @get <object> (field names...) [(reference names...)]"))
    end
end

# Based on Tim Holy's JuliaCon 2016 Keynote 
"""
A wrapper around an array that applies a function to any index element

# Fields
- `f`: the function to apply
- `data`: the underlying array
- `offset`: an optional offset (0 by default)

# Constructor

    MappedVector(f, data::AbstractVector{T}, T₀ = typeof(f(one(T))), offset = 0) where {T}

where `T₀` should be the type of the mapped values.

# Examples:

```julia
julia> x = [-π, 0.0, π];

julia> y = MappedVector(cos, x, Float64, 1) # apply `cos` on demand and start indexing from 0
Array{Float64,1} → Base.#cos

julia> y[0]
-1.0
```
"""
struct MappedVector{T,A <: AbstractVector,F} <: AbstractVector{T}
    f::F
    data::A
    offset::Int
end

Base.size(A::MappedVector) = size(A.data)
Base.@propagate_inbounds Base.getindex(A::MappedVector, i::Int) = A.f(A.data[i + A.offset])

function MappedVector(f, data::AbstractVector{T}, offset = 0) where {T}
    T₀ = typeof(f(one(T)))
    A  = typeof(data)
    F  = typeof(f)

    MappedVector{T₀, A, F}(f, data, offset)
end

function Base.show(io::IO, M::MIME"text/plain", m::MappedVector{T, A, F}) where {T, A, F}
    print(io, "$A → $F ($(1-m.offset):$(length(m.data)-m.offset))")
end
Base.show(io::IO, m::MappedVector) = Base.show(io, MIME("text/plain"), m)

end

