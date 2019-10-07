import Base: size, show, summary, +, -

abstract type PointData{N} <: AbstractVector{Float64} end

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
10-element Array{Float64,1}:
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
struct ScalarData{N} <: PointData{N}
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
10×2 Array{Float64,2}:
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
10×2 Array{Float64,2}:
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
struct VectorData{N} <: PointData{N}
    u::ScalarData{N}
    v::ScalarData{N}
end

"""
    TensorData

A wrapper for a one-dimensional array of 2x2 tensor-valued data, with fields
`dudx`, `dudy`, `dvdx`, `dvdy`. The resulting wrapper can be indexed as though these four components are stacked
on top of each other.

# Constructors
- `TensorData(dudx,dudy,dvdx,dvdy)` constructs a wrapper for the tensor components data.
- `TensorData(n::Int)` constructs a wrapper with zeros of length `n` for all components.
- `TensorData(x::ScalarData/VectorData/TensorData)` constructs a wrapper for zero components of the
   same length as that wrapped by `x`.

# Example

"""
struct TensorData{N} <: PointData{N}
    dudx::ScalarData{N}
    dudy::ScalarData{N}
    dvdx::ScalarData{N}
    dvdy::ScalarData{N}
end

function ScalarData(data::Vector{T}) where {T <: Real}
  ScalarData{length(data)}(convert.(Float64,data))
end

function VectorData(u::Vector{T},v::Vector{T}) where {T <: Real}
  @assert length(u) == length(v)
  VectorData{length(u)}(ScalarData(u),ScalarData(v))
end

function TensorData(dudx::Vector{T},dudy::Vector{T},
                    dvdx::Vector{T},dvdy::Vector{T}) where {T <: Real}
  @assert length(dudx) == length(dudy) == length(dvdx) == length(dvdy)
  TensorData{length(dudx)}(ScalarData(dudx),ScalarData(dudy),
                           ScalarData(dvdx),ScalarData(dvdy))
end

ScalarData(x::ScalarData) = ScalarData(zero(x.data))
ScalarData(n::Int) = ScalarData(zeros(Float64,n))
ScalarData(x::VectorData) = ScalarData(zero(x.u))
ScalarData(x::TensorData) = ScalarData(zero(x.dudx))
VectorData(x::Tuple{Vector{T},Vector{T}}) where {T <: Real} = VectorData(x[1],x[2])
VectorData(x::VectorData) = VectorData(zero(x.u),zero(x.v))
VectorData(n::Int) = VectorData(zeros(Float64,n),zeros(Float64,n))
VectorData(x::ScalarData) = VectorData(zero(x.data),zero(x.data))
VectorData(x::TensorData) = VectorData(zero(x.dudx),zero(x.dudy))
TensorData(x::TensorData) = TensorData(zero(x.dudx),zero(x.dudy),zero(x.dvdx),zero(x.dvdy))
TensorData(n::Int) = TensorData(zeros(Float64,n),zeros(Float64,n),zeros(Float64,n),zeros(Float64,n))
TensorData(x::ScalarData) = TensorData(zero(x.data),zero(x.data),zero(x.data),zero(x.data))
TensorData(x::VectorData) = TensorData(zero(x.u),zero(x.u),zero(x.v),zero(x.v))

(::Type{ScalarData{N}})() where {N} = ScalarData(N)
(::Type{VectorData{N}})() where {N} = VectorData(N)
(::Type{TensorData{N}})() where {N} = TensorData(N)


Base.similar(::ScalarData{N}) where {N} = ScalarData(N)

Base.similar(::VectorData{N}) where {N} = VectorData(N)

Base.similar(::TensorData{N}) where {N} = TensorData(N)



#Base.size(A::VectorData) = size(A.u).+size(A.v)
"""
    size(A::VectorData,d::Int) -> Int

Return twice the number of vector data points if `d` is 1 (the sum of the length of the `u` and
`v` vectors) and 1 if `d` is 2. This is consistent with the interpretation of VectorData as
a stacked pair of columns, corresponding to the `u` and `v` components, respectively.
"""
Base.size(::VectorData{N},d::Int) where {N} = d == 1 ? 2*N : 1

"""
    size(A::VectorData) -> Tuple

Return a tuple of the number of vector data points by the number of dimensions.
"""
Base.size(A::VectorData) = (size(A,1),)
@propagate_inbounds Base.getindex(A::VectorData{N},i::Int) where {N} =
   i > N ? A.v[i-N] : A.u[i]
@propagate_inbounds Base.setindex!(A::VectorData{N}, v, i::Int) where {N} =
   i > N ? A.v[i-N] = convert(Float64, v) : A.u[i] = convert(Float64, v)

"""
    size(A::TensorData,d::Int) -> Int

Return four times the number of tensor data points if `d` is 1 (the sum of the length of the four components)
and 1 if `d` is 2. This is consistent with the interpretation of TensorData as
a stacked set of columns.
"""
Base.size(::TensorData{N},d::Int) where {N} = d == 1 ? 4*N : 1

"""
    size(A::TensorData) -> Tuple

Return a tuple of the number of tensor data points by the number of dimensions.
"""
Base.size(A::TensorData) = (size(A,1),)

@propagate_inbounds Base.getindex(A::TensorData{N},i::Int) where {N} =
   i > N ? (i > 2*N ? (i > 3*N ? A.dvdy[i-3*N] : A.dvdx[i-2*N]) : A.dudy[i-N] ) : A.dudx[i]
@propagate_inbounds Base.setindex!(A::TensorData{N}, v, i::Int) where {N} =
   i > N ? (i > 2*N ? (i > 3*N ? A.dvdy[i-3*N] = convert(Float64, v) : A.dvdx[i-2*N] = convert(Float64, v) ) : A.dudy[i-N] = convert(Float64, v) ) : A.dudx[i] = convert(Float64, v)



function show(io::IO, m::MIME"text/plain", pts::ScalarData{N}) where {N}
  println(io,"$N points of scalar-valued data")
  show(io,m,pts.data)
end


function show(io::IO, m::MIME"text/plain", pts::VectorData{N}) where {N}
  println(io,"$N points of vector-valued data")
  show(io,m,hcat(pts.u,pts.v))
end

function show(io::IO, m::MIME"text/plain", pts::TensorData{N}) where {N}
  println(io,"$N points of tensor-valued data dudx, dudy, dvdx, dvdy")
  show(io,m,hcat(pts.dudx,pts.dudy,pts.dvdx,pts.dvdy))
end

include("basicpointoperations.jl")
