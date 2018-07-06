import Base: size

#abstract type Points <: AbstractArray{Float64} end

"""
    ScalarData

A wrapper for a one-dimensional array of scalar-valued data. The resulting wrapper
can be indexed in the same way as the array itself.

# Constructors
- `ScalarData(d)` constructs a wrapper for the one-dimensional array of data `d`
- `ScalarData(n::Int)` constructs a wrapper for an array of zeros of length `n`.
- `ScalarData(x::ScalarData)` constructs a wrapper for an array of zeros of the
   same length as that wrapped by `x`.
- `ScalarData(x::VectorData)` constructs a wrapper for an array of zeros of the
    same length as that wrapped by `x`.

# Example

```jldoctest
julia> f = ScalarData(10);

julia> f[5] = 1.0;

julia> f
10 points of scalar-valued data
 0.0
 0.0
 0.0
 0.0
 1.0
 0.0
 0.0
 0.0
 0.0
 0.0
```
"""
struct ScalarData{N} <: AbstractVector{Float64}
    data::Vector{Float64}
end

@wraparray ScalarData data

"""
    VectorData

A wrapper for a one-dimensional array of two-component vector-valued data. The
resulting wrapper can be indexed as though the first component and second
component are stacked on top of each other.

# Constructors
- `VectorData(u,v)` constructs a wrapper for the vector components data `u` and `v`.
- `VectorData(n::Int)` constructs a wrapper with zeros of length `n` for both components.
- `VectorData(x::ScalarData)` constructs a wrapper for zero components of the
   same length as that wrapped by `x`.
- `VectorData(x::VectorData)` constructs a wrapper for zero components of the
    same length as that wrapped by `x`.

# Example

```jldoctest
julia> f = VectorData(10);

julia> f.v[1:5] = 1:5;

julia> f
10 points of vector-valued data
 0.0  1.0
 0.0  2.0
 0.0  3.0
 0.0  4.0
 0.0  5.0
 0.0  0.0
 0.0  0.0
 0.0  0.0
 0.0  0.0
 0.0  0.0

julia> f[7] = 1; f[18] = 0.2;

julia> f
10 points of vector-valued data
 0.0  1.0
 0.0  2.0
 0.0  3.0
 0.0  4.0
 0.0  5.0
 0.0  0.0
 1.0  0.0
 0.0  0.2
 0.0  0.0
 0.0  0.0
```
"""
struct VectorData{N} <: AbstractVector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
end

Points = Union{ScalarData,VectorData}

function ScalarData(data::Vector{T}) where {T <: Real}
  ScalarData{length(data)}(convert.(Float64,data))
end

function VectorData(u::Vector{T},v::Vector{T}) where {T <: Real}
  @assert length(u) == length(v)
  VectorData{length(u)}(convert.(Float64,u),convert.(Float64,v))
end

ScalarData(x::ScalarData) = ScalarData(zeros(x.data))
ScalarData(n::Int) = ScalarData(zeros(Float64,n))
ScalarData(x::VectorData) = ScalarData(zeros(x.u))
VectorData(x::VectorData) = VectorData(zeros(x.u),zeros(x.v))
VectorData(n::Int) = VectorData(zeros(Float64,n),zeros(Float64,n))
VectorData(x::ScalarData) = VectorData(zeros(x.data),zeros(x.data))
(::Type{ScalarData{N}})() where {N} = ScalarData(N)
(::Type{VectorData{N}})() where {N} = VectorData(N)


Base.size(A::VectorData) = size(A.u).+size(A.v)
@propagate_inbounds Base.getindex(A::VectorData{N},i::Int) where {N} =
   i > N ? A.v[i-N] : A.u[i]
@propagate_inbounds Base.setindex!(A::VectorData{N}, v, i::Int) where {N} =
   i > N ? A.v[i-N] = convert(Float64, v) : A.u[i] = convert(Float64, v)


function Base.show(io::IO, pts::ScalarData{N}) where {N}
    println(io, "$N points of scalar-valued data")
end

function Base.show(io::IO, ::MIME"text/plain", pts::ScalarData{N}) where {N}
    println(io, "$N points of scalar-valued data")
    Base.showarray(io,pts.data,false;header=false)
end

function Base.show(io::IO, pts::VectorData{N}) where {N}
    println(io, "$N points of vector-valued data")
end

function Base.show(io::IO, ::MIME"text/plain", pts::VectorData{N}) where {N}
    println(io, "$N points of vector-valued data")
    Base.showarray(io,hcat(pts.u,pts.v),false;header=false)
end
