# Set it to negative of itself
function (-)(p_in::ScalarData)
  p = deepcopy(p_in)
  @. p.data = -p.data
  return p
end

function (-)(p_in::VectorData)
  p = deepcopy(p_in)
  @. p.u = -p.u
  @. p.v = -p.v
  return p
end

# Add and subtract the same type
function (-)(p1::T,p2::T) where {T<:ScalarData}
  return T(p1.data .- p2.data)
end

function (+)(p1::T,p2::T) where {T<:ScalarData}
  return T(p1.data .+ p2.data)
end

function (-)(p1::T,p2::T) where {T<:VectorData}
  return T(p1.u - p2.u, p1.v - p2.v)
end

function (+)(p1::T,p2::T) where {T<:VectorData}
  return T(p1.u + p2.u, p1.v + p2.v)
end

# Multiply and divide by a constant
function (*)(p::T,c::Number) where {T<:ScalarData}
  return T(c*p.data)
end


function (/)(p::T,c::Number) where {T<:ScalarData}
  return T(p.data ./ c)
end

function (*)(p::T,c::Number) where {T<:VectorData}
  return T(c*p.u,c*p.v)
end

(*)(c::Number,p::T) where {T<:PointData} = *(p,c)

function (/)(p::T,c::Number) where {T<:VectorData}
  return T(p.u / c, p.v / c)
end

function (*)(p1::T,p2::T) where {T<:ScalarData}
  return T(p1.data .* p2.data)
end

function (*)(p1::T,p2::T) where {T<:VectorData}
  return T(p1.u * p2.u, p1.v * p2.v)
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
   B = deepcopy(A)
    u, v = a
    B.u .+= u
    B.v .+= v
    return B
end

function (-)(A::VectorData,a::Tuple{T,T}) where {T <: Number}
   B = deepcopy(A)
    u, v = a
    B.u .-= u
    B.v .-= v
    return B
end

a::(Tuple{T,T} where {T}) + A::VectorData = A + a
a::(Tuple{T,T} where {T}) - A::VectorData = (B = A - a; @. B.u = -B.u; @. B.v = -B.v; return B)

"""
    cross(a::Number/ScalarData,A::VectorData) -> VectorData
    ×(a::Number/ScalarData,A::VectorData) -> VectorData

Compute the cross product between the scalar `a` (treated as an out-of-plane component of a vector)
and the planar vector data `A`.
"""
function cross(a::Union{Number,ScalarData},A::VectorData)
    B = deepcopy(A)
    @. B.u = -a*A.v
    @. B.v = a*A.u
    return B
end

"""
    cross(A::VectorData,B::VectorData) -> ScalarData
    ×(A::VectorData,A::VectorData) -> ScalarData

Compute the cross product between the vector point data `A` and `B`
and return the result as scalar data (treated as an out-of-plane
component of a vector).
"""
function cross(A::VectorData{N},B::VectorData{N}) where N
    C = ScalarData(N)
    @. C = A.u*B.v - A.v*B.u
    return C
end

"""
    dot(v::Tuple{T,T},B::TensorData) where {T<:Number} -> VectorData
    ⋅(v::Tuple{T,T},B::TensorData) where {T<:Number} -> VectorData

Computes the dot product between the tuple `v` and the elements of a tensor `B` on
a set of points and returns vector data on the same set of points.
"""
function dot(A::Tuple{T,T},B::TensorData) where {T<:Number}
    C = VectorData(B)
    x, y = A
    @. C.u = x*B.dudx + y*B.dudy
    @. C.v = x*B.dvdx + y*B.dvdy
    return C
end
