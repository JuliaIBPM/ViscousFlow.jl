import Base: size, show, summary, similar, parent, parentindices

abstract type PointData{N,T} <: AbstractVector{T} end

const NDIM = 2

"""
    ScalarData <: PointData

A wrapper for a one-dimensional array of scalar-valued data. The resulting wrapper
can be indexed in the same way as the array itself.

# Constructors
- `ScalarData(d::AbstractVector[,dtype=Float64])` constructs a wrapper for the one-dimensional array of data `d`
- `ScalarData(n::Int)` constructs a wrapper for an array of zeros of length `n`.
- `ScalarData(x::PointData)` constructs a wrapper for an array of zeros of the
   same length as that wrapped by `x`.
- `ScalarData(n::Int,dtype=ComplexF64)` constructs a wrapper for complex-valued data.
- `ScalarData(x::PointData,dtype=ComplexF64)` constructs a wrapper for an array of
   complex zeros of the same length as that wrapped by `x`.

# Example

```jldoctest
julia> f = ScalarData(10);

julia> f[5] = 1.0;

julia> f
10 points of scalar-valued Float64 data
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
struct ScalarData{N,T,DT<:AbstractVector} <: PointData{N,T}
    data::DT
    ScalarData{N,T,DT}(a::S) where {N,T,DT<:AbstractVector,S<:PointData} = new{length(a),T,typeof(a.data)}(a.data)
    ScalarData{N,T,DT}(a::S) where {N,T,DT<:AbstractVector,S<:AbstractVector} = new{length(a),T,typeof(a)}(a)
end

#@wraparray ScalarData data 1

function ScalarData(data::AbstractVector{T}) where {T <: Number}
  ScalarData{length(data),T,typeof(data)}(data)
end

ScalarData(n::Int;dtype=Float64) = ScalarData(zeros(dtype,n))



"""
    VectorData <: PointData

A wrapper for a one-dimensional array of two-component vector-valued data. The
resulting wrapper can be indexed as though the first component and second
component are stacked on top of each other.

# Constructors
- `VectorData(d::AbstractVector[,dtype=Float64])` constructs a wrapper for the one-dimensional array of data `d`, splitting `d` into the `u` and `v` components evenly.
- `VectorData(u::AbstractVector,v::AbstractVector)` constructs a wrapper for the vector components data `u` and `v`.
- `VectorData(n::Int)` constructs a wrapper with zeros of length `n` for both components.
- `VectorData(x::PointData)` constructs a wrapper for zero components of the
   same length as that wrapped by `x`.
- `VectorData(n::Int,dtype=ComplexF64)` constructs a wrapper with complex-valued zeros
   of length `n` for both components.

# Example

```jldoctest
julia> f = VectorData(10,dtype=ComplexF64);

julia> f.v[1:5] = 1:5;

julia> f
10 points of vector-valued Complex{Float64} data
10×2 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im
 0.0+0.0im  2.0+0.0im
 0.0+0.0im  3.0+0.0im
 0.0+0.0im  4.0+0.0im
 0.0+0.0im  5.0+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> f[7] = 1im; f[18] = 0.2;

julia> f
10 points of vector-valued Complex{Float64} data
10×2 Array{Complex{Float64},2}:
 0.0+0.0im  1.0+0.0im
 0.0+0.0im  2.0+0.0im
 0.0+0.0im  3.0+0.0im
 0.0+0.0im  4.0+0.0im
 0.0+0.0im  5.0+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+1.0im  0.0+0.0im
 0.0+0.0im  0.2+0.0im
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
```
"""
struct VectorData{N,T,DT<:AbstractVector} <: PointData{N,T}
    data::DT
    u::ScalarData{N,T}
    v::ScalarData{N,T}
end

function VectorData(data::AbstractVector{T}) where {T <: Number}
  nc = NDIM
  @assert mod(length(data),nc) == 0
  numpts = div(length(data),nc)
  u = ScalarData(view(data,1:numpts))
  v = ScalarData(view(data,numpts+1:length(data)))
  VectorData{numpts,T,typeof(data)}(data,u,v)
end

function VectorData(u::AbstractVector{T},v::AbstractVector{T}) where {T <: Number}
  @assert length(u) == length(v)
  VectorData([u;v])
end

VectorData(x::Tuple{AbstractVector{T},AbstractVector{T}}) where {T <: Number} = VectorData(x...)
VectorData(n::Int;dtype=Float64) = VectorData(zeros(dtype,NDIM*n))

"""
    TensorData <: PointData

A wrapper for a one-dimensional array of 2x2 tensor-valued data, with fields
`dudx`, `dudy`, `dvdx`, `dvdy`. The resulting wrapper can be indexed as though these four components are stacked
on top of each other.

# Constructors
- `TensorData(d::AbstractVector[,dtype=Float64])` constructs a wrapper for the one-dimensional array of data `d`, splitting `d` into the four components evenly.
- `TensorData(dudx,dudy,dvdx,dvdy)` constructs a wrapper for the tensor components data, each of type `AbstractVector`
- `TensorData(n::Int)` constructs a wrapper with zeros of length `n` for all components.
- `TensorData(x::PointData[,dtype=Float64])` constructs a wrapper for zero components of the
   same length as that wrapped by `x`.

# Example

"""
struct TensorData{N,T,DT<:AbstractVector} <: PointData{N,T}
    data::DT
    dudx::ScalarData{N,T}
    dudy::ScalarData{N,T}
    dvdx::ScalarData{N,T}
    dvdy::ScalarData{N,T}
end

function TensorData(data::AbstractVector{T}) where {T <: Number}
  nc = NDIM*NDIM
  @assert mod(length(data),nc) == 0
  numpts = div(length(data),nc)
  dudx = ScalarData(view(data,         1:numpts))
  dudy = ScalarData(view(data,  numpts+1:2*numpts))
  dvdx = ScalarData(view(data,2*numpts+1:3*numpts))
  dvdy = ScalarData(view(data,3*numpts+1:4*numpts))
  TensorData{numpts,T,typeof(data)}(data,dudx,dudy,dvdx,dvdy)
end

function TensorData(dudx::AbstractVector{T},dudy::AbstractVector{T},
                dvdx::AbstractVector{T},dvdy::AbstractVector{T}) where {T <: Number}
  @assert length(dudx) == length(dudy) == length(dvdx) == length(dvdy)
  TensorData([dudx;dudy;dvdx;dvdy])
end


TensorData(x::NTuple{NDIM*NDIM,AbstractVector{T}}) where {T <: Number} = TensorData(x...)
TensorData(n::Int;dtype=Float64) = TensorData(zeros(dtype,n*NDIM*NDIM))

for f in (:ScalarData,:VectorData,:TensorData)
  @eval (::Type{$f{N,T,DT}})() where {N,T<:Number,DT<:AbstractVector} = $f(N,dtype=T)
  @eval $f(::PointData{N,T};dtype=T) where {N,T<:Number} = $f(N,dtype=dtype)
  @eval similar(::$f{N,T,DT};element_type=T) where {N,T<:Number,DT<:AbstractVector} = $f(N,dtype=element_type)
end

for f in (:VectorData,:TensorData)
  # if argument u is an AbstractVector, then make this new instance a point to it
  @eval (::Type{$f{N,T,DT}})(u::AbstractVector) where {N,T<:Number,DT<:AbstractVector} = $f(u)
  # if argument u is the same PointData type, then make this new instance
  # a pointer to the `data` field of u
  @eval (::Type{$f{N,T,DT}})(u::$f{N,T,DT}) where {N,T<:Number,DT<:AbstractVector} = $f(u.data)
end


parent(A::PointData) = A.data
parentindices(A::PointData) = parentindices(A.data)

@propagate_inbounds Base.getindex(A::PointData,i::Int) = getindex(A.data,i)
@propagate_inbounds Base.setindex!(A::PointData, v, i::Int) = setindex!(A.data,v,i)

size(A::PointData{N,T},d::Int) where {N,T} = size(A.data,d)
size(A::PointData{N,T}) where {N,T} = size(A.data)


#=
"""
    size(A::PointData,d::Int) -> Int

Return twice the number of vector data points if `d` is 1 (the sum of the length of the `u` and
`v` vectors) and 1 if `d` is 2. This is consistent with the interpretation of VectorData as
a stacked pair of columns, corresponding to the `u` and `v` components, respectively.
"""
Base.size(::VectorData{N,T},d::Int) where {N,T} = d == 1 ? NDIM*N : 1

"""
    size(A::VectorData) -> Tuple

Return a tuple of the number of vector data points by the number of dimensions.
"""
Base.size(A::VectorData) = (size(A,1),)



#@propagate_inbounds Base.getindex(A::VectorData{N,T},i::Int) where {N,T} =
#   i > N ? A.data[i-N] : A.u[i]
#@propagate_inbounds Base.setindex!(A::VectorData{N,T}, v, i::Int) where {N,T} =
#   i > N ? A.v[i-N] = convert(T, v) : A.u[i] = convert(T, v)

"""
    size(A::TensorData,d::Int) -> Int

Return four times the number of tensor data points if `d` is 1 (the sum of the length of the four components)
and 1 if `d` is 2. This is consistent with the interpretation of TensorData as
a stacked set of columns.
"""
Base.size(::TensorData{N,T},d::Int) where {N,T} = d == 1 ? NDIM*NDIM*N : 1

"""
    size(A::TensorData) -> Tuple

Return a tuple of the number of tensor data points by the number of dimensions.
"""
Base.size(A::TensorData) = (size(A,1),)
=#

function show(io::IO, m::MIME"text/plain", pts::ScalarData{N,T}) where {N,T}
  println(io,"$N points of scalar-valued $T data")
  show(io,m,pts.data)
end

function show(io::IO, m::MIME"text/plain", pts::VectorData{N,T}) where {N,T}
  println(io,"$N points of vector-valued $T data")
  show(io,m,pts.data)
end

function show(io::IO, m::MIME"text/plain", pts::TensorData{N,T}) where {N,T}
  println(io,"$N points of tensor-valued $T data")
  show(io,m,pts.data)
end


include("basicpointoperations.jl")
