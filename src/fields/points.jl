import Base: size, show, summary, +, -

abstract type PointData{N,T} <: AbstractVector{T} end

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
struct ScalarData{N,T} <: PointData{N,T}
    data::Vector{T}
end

@wraparray ScalarData data 1

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
struct VectorData{N,T} <: PointData{N,T}
    u::ScalarData{N,T}
    v::ScalarData{N,T}
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
struct TensorData{N,T} <: PointData{N,T}
    dudx::ScalarData{N,T}
    dudy::ScalarData{N,T}
    dvdx::ScalarData{N,T}
    dvdy::ScalarData{N,T}
end

function ScalarData(data::Vector{T}) where {T <: Number}
  ScalarData{length(data),T}(data)
end

function VectorData(u::Vector{T},v::Vector{T}) where {T <: Number}
  @assert length(u) == length(v)
  VectorData{length(u),T}(ScalarData(u),ScalarData(v))
end

function TensorData(dudx::Vector{T},dudy::Vector{T},
                    dvdx::Vector{T},dvdy::Vector{T}) where {T <: Number}
  @assert length(dudx) == length(dudy) == length(dvdx) == length(dvdy)
  TensorData{length(dudx),T}(ScalarData(dudx),ScalarData(dudy),
                           ScalarData(dvdx),ScalarData(dvdy))
end

ScalarData(x::PointData{N,T}) where {N,T} = ScalarData(zeros(T,N))
ScalarData(n::Int;dtype=Float64) = ScalarData(zeros(dtype,n))
#ScalarData(x::VectorData) = ScalarData(zero(x.u))
#ScalarData(x::TensorData) = ScalarData(zero(x.dudx))
VectorData(x::Tuple{Vector{T},Vector{T}}) where {T <: Number} = VectorData(x[1],x[2])
VectorData(x::PointData{N,T}) where {N,T} = VectorData(zeros(T,N),zeros(T,N))
VectorData(n::Int;dtype=Float64) = VectorData(zeros(dtype,n),zeros(dtype,n))
#VectorData(x::ScalarData) = VectorData(zero(x.data),zero(x.data))
#VectorData(x::TensorData) = VectorData(zero(x.dudx),zero(x.dudy))
TensorData(x::PointData{N,T}) where {N,T} = TensorData(zeros(T,N),zeros(T,N),zeros(T,N),zeros(T,N))
TensorData(n::Int;dtype=Float64) = TensorData(zeros(dtype,n),zeros(dtype,n),zeros(dtype,n),zeros(dtype,n))
#TensorData(x::ScalarData) = TensorData(zero(x.data),zero(x.data),zero(x.data),zero(x.data))
#TensorData(x::VectorData) = TensorData(zero(x.u),zero(x.u),zero(x.v),zero(x.v))

(::Type{ScalarData{N,T}})() where {N,T} = ScalarData(N,dtype=T)
(::Type{VectorData{N,T}})() where {N,T} = VectorData(N,dtype=T)
(::Type{TensorData{N,T}})() where {N,T} = TensorData(N,dtype=T)


Base.similar(::ScalarData{N,T}) where {N,T} = ScalarData(N,dtype=T)

Base.similar(::VectorData{N,T}) where {N,T} = VectorData(N,dtype=T)

Base.similar(::TensorData{N,T}) where {N,T} = TensorData(N,dtype=T)



#Base.size(A::VectorData) = size(A.u).+size(A.v)
"""
    size(A::VectorData,d::Int) -> Int

Return twice the number of vector data points if `d` is 1 (the sum of the length of the `u` and
`v` vectors) and 1 if `d` is 2. This is consistent with the interpretation of VectorData as
a stacked pair of columns, corresponding to the `u` and `v` components, respectively.
"""
Base.size(::VectorData{N,T},d::Int) where {N,T} = d == 1 ? 2*N : 1

"""
    size(A::VectorData) -> Tuple

Return a tuple of the number of vector data points by the number of dimensions.
"""
Base.size(A::VectorData) = (size(A,1),)
@propagate_inbounds Base.getindex(A::VectorData{N,T},i::Int) where {N,T} =
   i > N ? A.v[i-N] : A.u[i]
@propagate_inbounds Base.setindex!(A::VectorData{N,T}, v, i::Int) where {N,T} =
   i > N ? A.v[i-N] = convert(T, v) : A.u[i] = convert(T, v)

"""
    size(A::TensorData,d::Int) -> Int

Return four times the number of tensor data points if `d` is 1 (the sum of the length of the four components)
and 1 if `d` is 2. This is consistent with the interpretation of TensorData as
a stacked set of columns.
"""
Base.size(::TensorData{N,T},d::Int) where {N,T} = d == 1 ? 4*N : 1

"""
    size(A::TensorData) -> Tuple

Return a tuple of the number of tensor data points by the number of dimensions.
"""
Base.size(A::TensorData) = (size(A,1),)

@propagate_inbounds Base.getindex(A::TensorData{N,T},i::Int) where {N,T} =
   i > N ? (i > 2*N ? (i > 3*N ? A.dvdy[i-3*N] : A.dvdx[i-2*N]) : A.dudy[i-N] ) : A.dudx[i]
@propagate_inbounds Base.setindex!(A::TensorData{N,T}, v, i::Int) where {N,T} =
   i > N ? (i > 2*N ? (i > 3*N ? A.dvdy[i-3*N] = convert(T, v) : A.dvdx[i-2*N] = convert(T, v) ) : A.dudy[i-N] = convert(T, v) ) : A.dudx[i] = convert(T, v)



function show(io::IO, m::MIME"text/plain", pts::ScalarData{N,T}) where {N,T}
  println(io,"$N points of scalar-valued $T data")
  show(io,m,pts.data)
end


function show(io::IO, m::MIME"text/plain", pts::VectorData{N,T}) where {N,T}
  println(io,"$N points of vector-valued $T data")
  show(io,m,hcat(pts.u,pts.v))
end

function show(io::IO, m::MIME"text/plain", pts::TensorData{N,T}) where {N,T}
  println(io,"$N points of tensor-valued $T data dudx, dudy, dvdx, dvdy")
  show(io,m,hcat(pts.dudx,pts.dudy,pts.dvdx,pts.dvdy))
end

include("basicpointoperations.jl")
