import Base: +, -, ∘, real, imag, abs

# Set it to negative of itself
function (-)(p_in::PointData)
  p = similar(p_in)
  @. p.data = -p_in.data
  return p
end

# Add and subtract the same type
function (-)(p1::T,p2::T) where {T<:PointData}
  q = similar(p1,element_type=promote_type(eltype(p1),eltype(p2)))
  @. q.data = p1.data - p2.data
  return q
end

function (+)(p1::T,p2::T) where {T<:PointData}
  q = similar(p1,element_type=promote_type(eltype(p1),eltype(p2)))
  @. q.data = p1.data + p2.data
  return q
end


# Multiply and divide by a constant
function (*)(p::T,c::Number) where {T<:PointData}
  q = similar(p,element_type=promote_type(typeof(c),eltype(p)))
  @. q.data = c*p.data
  return q
end


function (/)(p::T,c::Number) where {T<:PointData}
  q = similar(p,element_type=promote_type(typeof(c),eltype(p)))
  @. q.data = p.data/c
  return q
end

(*)(c::Number,p::T) where {T<:PointData} = *(p,c)


function (*)(p1::T,p2::T) where {T<:PointData}
  q = similar(p1,element_type=promote_type(eltype(p1),eltype(p2)))
  @. q.data = p1.data * p2.data
  return q
end


"""
    (+)(X::VectorData,a::Tuple{T,T}) where {T<:Number} -> VectorData
    (-)(X::VectorData,a::Tuple{T,T}) where {T<:Number} -> VectorData

Adds or subtracts the tuple `a` component by component to each element of `X`. All data in `a` are
converted to Float64. Can also switch the arguments.

# Example

```jldoctest
julia> f = VectorData(5);

julia> f + (2,3)
5 points of vector-valued Float64 data
5×2 Array{Float64,2}:
 2.0  3.0
 2.0  3.0
 2.0  3.0
 2.0  3.0
 2.0  3.0
```
"""
function (+)(A::VectorData,a::Tuple{T,T}) where {T <: Number}
    B = similar(A,element_type=promote_type(T,eltype(A)))
    u, v = a
    B.u .+= u
    B.v .+= v
    return B
end

function (-)(A::VectorData,a::Tuple{T,T}) where {T <: Number}
  B = similar(A,element_type=promote_type(T,eltype(A)))
    u, v = a
    B.u .-= u
    B.v .-= v
    return B
end

### Hadamard products

"""
    product!(out::PointData,p::PointData,q::PointData)

Compute the Hadamard (i.e. element by element) product of point data
data `p` and `q` and return the result in `out`. Note that `p` and `q`
can be of mixed type (scalar, vector, tensor), as long as one of them is
a scalar. Also, `out` must have element type that is consistent with the promoted
type of `p` and `q`.

# Example
```jldoctest
julia> fcs = ScalarData(5,dtype=ComplexF64);

julia> fill!(fcs,2im)
5 points of scalar-valued Complex{Float64} data
5-element Array{Complex{Float64},1}:
 0.0 + 2.0im
 0.0 + 2.0im
 0.0 + 2.0im
 0.0 + 2.0im
 0.0 + 2.0im

julia> frt = TensorData(fcs,dtype=Float64);

julia> fill!(frt,1.0)
5 points of tensor-valued Float64 data dudx, dudy, dvdx, dvdy
5×4 Array{Float64,2}:
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0

julia> out = similar(frt,element_type=ComplexF64);

julia> product!(out,frt,fcs)
5 points of tensor-valued Complex{Float64} data dudx, dudy, dvdx, dvdy
5×4 Array{Complex{Float64},2}:
 0.0+2.0im  0.0+2.0im  0.0+2.0im  0.0+2.0im
 0.0+2.0im  0.0+2.0im  0.0+2.0im  0.0+2.0im
 0.0+2.0im  0.0+2.0im  0.0+2.0im  0.0+2.0im
 0.0+2.0im  0.0+2.0im  0.0+2.0im  0.0+2.0im
 0.0+2.0im  0.0+2.0im  0.0+2.0im  0.0+2.0im
```
"""
function product!(out::PointData{N},
                  p::PointData{N},
                  q::PointData{N}) where {N}

    @inbounds for x in 1:N
        out.data[x] = p.data[x] * q.data[x]
    end
    out
end

# function product!(out::VectorData{N},
#                   p::VectorData{N},
#                   q::VectorData{N}) where {N}
#
#     @inbounds for x in 1:N
#         out.u[x] = p.u[x] * q.u[x]
#         out.v[x] = p.v[x] * q.v[x]
#     end
#     out
# end

# function product!(out::TensorData{N},
#                   p::TensorData{N},
#                   q::TensorData{N}) where {N}
#
#     @inbounds for x in 1:N
#         out.dudx[x] = p.dudx[x] * q.dudx[x]
#         out.dudy[x] = p.dudy[x] * q.dudy[x]
#         out.dvdx[x] = p.dvdx[x] * q.dvdx[x]
#         out.dvdy[x] = p.dvdy[x] * q.dvdy[x]
#     end
#     out
# end

function product!(out::VectorData{N},
                  p::ScalarData{N},
                  q::VectorData{N}) where {N}

    @inbounds for x in 1:N
        out.u[x] = p[x] * q.u[x]
        out.v[x] = p[x] * q.v[x]
    end
    out
end

function product!(out::TensorData{N},
                  p::ScalarData{N},
                  q::TensorData{N}) where {N}

    @inbounds for x in 1:N
        out.dudx[x] = p[x] * q.dudx[x]
        out.dudy[x] = p[x] * q.dudy[x]
        out.dvdx[x] = p[x] * q.dvdx[x]
        out.dvdy[x] = p[x] * q.dvdy[x]
    end
    out
end

"""
    product(p::PointData,q::PointData) -> PointData
    (∘)(p::PointData,q::PointData) -> PointData

Compute the Hadamard (i.e. element by element) product of point data
data `p` and `q`. Works similarly to `product!`.
"""
function product end

for f in (:ScalarData, :VectorData, :TensorData)
    @eval function product(p::$f{N}, q::$f{N}) where {N}
        product!($f(N,dtype=promote_type(eltype(p),eltype(q))), p, q)
    end

    @eval (∘)(p::$f{N}, q::$f{N}) where {N} = product(p,q)
end

for f in (:VectorData, :TensorData)
    @eval product!(out::$f{N},p::$f{N},q::ScalarData{N}) where {N} = product!(out,q,p)

    @eval function product(p::$f{N},q::ScalarData{N}) where {N}
        product!($f(N,dtype=promote_type(eltype(p),eltype(q))), p, q)
    end

    @eval product(q::ScalarData{N},p::$f{N}) where {N} = product(p,q)

    @eval (∘)(p::$f{N}, q::ScalarData{N}) where {N} = product(p,q)

    @eval (∘)(q::ScalarData{N}, p::$f{N}) where {N} = product(p,q)

end

### Operations between tuples and vectors

a::(Tuple{T,T} where {T}) + A::VectorData = A + a
a::(Tuple{T,T} where {T}) - A::VectorData = (B = A - a; @. B.data = -B.data; return B)

"""
    cross(a::Number/ScalarData,A::VectorData) -> VectorData
    ×(a::Number/ScalarData,A::VectorData) -> VectorData

Compute the cross product between the scalar `a` (treated as an out-of-plane component of a vector)
and the planar vector data `A`.
"""
function cross(a::Union{Number,ScalarData},A::VectorData)
    B = similar(A)
    @. B.u = -a*A.v
    @. B.v = a*A.u
    return B
end

### Vector operations

"""
    cross(A::VectorData,B::VectorData) -> ScalarData
    ×(A::VectorData,A::VectorData) -> ScalarData

Compute the cross product between the vector point data `A` and `B`
and return the result as scalar data (treated as an out-of-plane
component of a vector).
"""
function cross(A::VectorData{N},B::VectorData{N}) where N
    C = ScalarData(N,dtype=promote_type(eltype(A),eltype(B)))
    @. C = A.u*B.v - A.v*B.u
    return C
end

"""
    dot(A::Tuple{T,T},B::VectorData) where {T<:Number} -> ScalarData
    ⋅(A::Tuple{T,T},B::VectorData) where {T<:Number} -> ScalarData

Computes the dot product between the tuple `v` and the elements of a tensor `B` on
a set of points and returns scalar data on the same set of points.
"""
function dot(A::Tuple{T,T},B::VectorData) where {T<:Number}
    C = ScalarData(B)
    x, y = A
    @. C.data = x*B.u + y*B.v
    return C
end

"""
    dot(A::Tuple{T,T},B::TensorData) where {T<:Number} -> VectorData
    ⋅(A::Tuple{T,T},B::TensorData) where {T<:Number} -> VectorData

Computes the dot product between the tuple `A` and the elements of a tensor `B` on
a set of points and returns vector data on the same set of points.
"""
function dot(A::Tuple{T,T},B::TensorData) where {T<:Number}
    C = VectorData(B)
    x, y = A
    @. C.u = x*B.dudx + y*B.dudy
    @. C.v = x*B.dvdx + y*B.dvdy
    return C
end

### Operations on complex point data

for f in (:real, :imag, :abs)
  @eval function $f(A::PointData{N,T}) where {N,T <: ComplexF64}
      Acopy = similar(A,element_type=Float64)
      Acopy .= broadcast($f,A)
      return Acopy
  end
end
