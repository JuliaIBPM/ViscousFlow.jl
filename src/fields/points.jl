import Base: size

abstract type Points end

struct ScalarData{N} <: Points
    data::Vector{Float64}
end

struct VectorData{N} <: Points
    u::Vector{Float64}
    v::Vector{Float64}
end

function Points(data::Vector{T}) where {T <: Real}
  ScalarData{length(data)}(convert.(Float64,data))
end

function Points(u::Vector{T},v::Vector{T}) where {T <: Real}
  @assert length(u) == length(v)
  VectorData{length(u)}(convert.(Float64,u),convert.(Float64,v))
end

Points(x::ScalarData) = Points(zeros(x.data))

Points(x::VectorData) = Points(zeros(x.u),zeros(x.v))

ScalarData(x::VectorData) = Points(zeros(x.u))

VectorData(x::ScalarData) = Points(zeros(x.data),zeros(x.data))

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
