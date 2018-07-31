module TimeMarching

import Whirl:@get

export System, Constrained, Unconstrained, RK, IFRK, IFHERK, r₁, r₂, B₂, B₁ᵀ

"Abstract type for a system of ODEs"
abstract type System{C} end

const Constrained = true
const Unconstrained = false

# These extend functions to tuples of abstract arrays
#=
Base.broadcast!(f,a::Tuple{T1,T2},b::Tuple{T1,T2}) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
    (broadcast!(f,a[1],b[1]); broadcast!(f,a[2],b[2]); a)
Base.broadcast!(f,a::Tuple{T1,T2},c::Number,b::Tuple{T1,T2}) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
    (broadcast!(f,a[1],c,b[1]); broadcast!(f,a[2],c,b[2]); a)
Base.broadcast(f,a::Tuple{T1,T2},b::Tuple{T1,T2}) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
    broadcast(f,a[1],b[1]), broadcast(f,a[2],b[2])
Base.broadcast(f,a::Number,b::Tuple{T1,T2}) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
    broadcast(f,a,b[1]), broadcast(f,a,b[2])
Base.broadcast!(f,a::AbstractArray{Tuple{T1,T2}},c::Number,b::Tuple{T1,T2}...) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
(@inbounds for I in eachindex(a); a[I] = f.(c,b...); end; a)
Base.broadcast!(f,a::AbstractArray{Tuple{T1,T2}},b::Tuple{T1,T2}...) where {T1<:Union{Tuple,AbstractArray},T2<:Union{Tuple,AbstractArray}} =
(@inbounds for I in eachindex(a); a[I] = f.(b...); end; a)
=#

# Functions that get extended by individual systems
function r₁ end
function r₂ end
function B₂ end
function B₁ᵀ end

struct RKParams{N}
  c::Vector{Float64}
  a::Matrix{Float64}
end

using ..SaddlePointSystems

const RK31 = RKParams{3}([0.5, 1.0, 1.0],
                      [1/2        0        0
                       √3/3 (3-√3)/3        0
                       (3+√3)/6    -√3/3 (3+√3)/6])

const Euler = RKParams{1}([1.0],ones(1,1))

include("timemarching/rk.jl")
include("timemarching/ifrk.jl")
include("timemarching/ifherk.jl")

end
