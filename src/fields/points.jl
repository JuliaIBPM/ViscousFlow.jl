import Base: size, show, summary, similar

abstract type PointData{N,T} <: AbstractVector{T} end

"""
    ScalarData <: PointData

A wrapper for a one-dimensional array of scalar-valued data. The resulting wrapper
can be indexed in the same way as the array itself.

# Constructors
- `ScalarData(d[,dtype=Float64])` constructs a wrapper for the one-dimensional array of data `d`
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
struct ScalarData{N,T} <: PointData{N,T}
    data::Vector{T}
end

@wraparray ScalarData data 1

"""
    VectorData <: GridData

A wrapper for a one-dimensional array of two-component vector-valued data. The
resulting wrapper can be indexed as though the first component and second
component are stacked on top of each other.

# Constructors
- `VectorData(u,v)` constructs a wrapper for the vector components data `u` and `v`.
- `VectorData(n::Int)` constructs a wrapper with zeros of length `n` for both components.
- `VectorData(x::GridData)` constructs a wrapper for zero components of the
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

ScalarData(n::Int;dtype=Float64) = ScalarData(zeros(dtype,n))

VectorData(x::Tuple{Vector{T},Vector{T}}) where {T <: Number} = VectorData(x[1],x[2])
VectorData(n::Int;dtype=Float64) = VectorData(zeros(dtype,n),zeros(dtype,n))

TensorData(x::NTuple{4,Vector{T}}) where {T <: Number} = TensorData(x[1],x[2],x[3],x[4])
TensorData(n::Int;dtype=Float64) = TensorData(zeros(dtype,n),zeros(dtype,n),zeros(dtype,n),zeros(dtype,n))

for f in (:ScalarData,:VectorData,:TensorData)
  @eval (::Type{$f{N,T}})() where {N,T} = $f(N,dtype=T)
  @eval $f(::PointData{N,T};dtype=T) where {N,T} = $f(N,dtype=dtype)
  @eval similar(::$f{N,T};element_type=T) where {N,T} = $f(N,dtype=element_type)
end


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
