#=
 Methods and types associated with single and double layer potentials
=#

module Layers

using CartesianGrids
using RigidBodyTools

export DoubleLayer, SingleLayer, MaskType, Mask, ComplementaryMask

abstract type LayerType{N,NX,NY} end

struct DoubleLayer{N,NX,NY,G,T,DT,DDT} <: LayerType{N,NX,NY}
    nds :: VectorData{N,Float64,DT}
    H :: RegularizationMatrix{Edges{G,NX,NY,T,DDT},VectorData{N,T,DT}}
end

function DoubleLayer(body::Union{Body,BodyList},H::RegularizationMatrix;weight::Float64 = 1.0)
  ds = ScalarData(dlengthmid(body))*weight
  normals = VectorData(normalmid(body))
  return DoubleLayer(normals∘ds,H)
end

function DoubleLayer(body::Union{Body,BodyList},reg::Regularize{N},w::GridData{NX,NY,T}) where {N,NX,NY,T}

  out = RegularizationMatrix(reg,VectorData(N,dtype=T),Edges(celltype(w),w,dtype=T))
  return DoubleLayer(body,reg._issymmetric ? out[1] : out, weight = sqrt(reg.overdv))
end

(μ::DoubleLayer{N})(p::ScalarData{N}) where {N} = divergence(μ.H*(p∘μ.nds))

function (μ::DoubleLayer{N,NX,NY,G,T,DT,DDT})(p::Number) where {N,NX,NY,G,T,DT,DDT}
  ϕ = ScalarData(N,dtype=T)
  ϕ .= p
  return μ(ϕ)
end

function Base.show(io::IO, H::DoubleLayer{N,NX,NY,G,T,DT,DDT}) where {N,NX,NY,G,T,DT,DDT}
    println(io, "Double-layer potential mapping")
    println(io, "  from $N scalar-valued point data of $T type")
    println(io, "  to a $NX x $NY grid of $G nodal data")
end

struct SingleLayer{N,NX,NY,G,T,DT,DDT} <: LayerType{N,NX,NY}
    ds :: ScalarData{N,Float64,DT}
    H :: RegularizationMatrix{Nodes{G,NX,NY,T,DDT},ScalarData{N,T,DT}}
end

function SingleLayer(body::Union{Body,BodyList},H::RegularizationMatrix)
  ds = ScalarData(dlengthmid(body))
  return SingleLayer(ds,H)
end

function SingleLayer(body::Union{Body,BodyList},reg::Regularize{N},w::GridData{NX,NY,T}) where {N,NX,NY,T}

  out = RegularizationMatrix(reg,ScalarData(N,dtype=T),Nodes(celltype(w),w,dtype=T))
  return SingleLayer(body,reg._issymmetric ? out[1] : out)
end

(μ::SingleLayer{N})(p::ScalarData{N}) where {N} = μ.H*(p∘μ.ds)

function (μ::SingleLayer{N,NX,NY,G,T,DT,DDT})(p::Number) where {N,NX,NY,G,T,DT,DDT}
  ϕ = ScalarData(N,dtype=T)
  ϕ .= p
  return μ(ϕ)
end


function Base.show(io::IO, H::SingleLayer{N,NX,NY,G,T,DT}) where {N,NX,NY,G,T,DT}
    println(io, "Single-layer potential mapping")
    println(io, "  from $N scalar-valued point data of $T type")
    println(io, "  to a $NX x $NY grid of $G nodal data")
end

abstract type MaskType end

struct Mask{N,NX,NY,G} <: MaskType
  data :: Nodes{G,NX,NY}
end

struct ComplementaryMask{N,NX,NY,G} <: MaskType
  mask :: Mask{N,NX,NY,G}
end

"""
    Mask(b::Body,reg::Regularize,w::GridData)

Returns a `Nodes` mask that sets every grid point equal to 1 if it lies inside the body
`b` and 0 outside. Returns grid data of the same type as `w`. Uses regularization
defined by `reg`.
"""
function Mask(dlayer::DoubleLayer{N,NX,NY,G,T,DT,DDT}) where {N,NX,NY,G,T,DT,DDT}
  L = plan_laplacian(NX,NY,with_inverse=true,dtype=T)
  return Mask{N,NX,NY,G}(L\dlayer(1))
end

Mask(body::Union{Body,BodyList},regop::Regularize{N},w::GridData{NX,NY,T}) where {N,NX,NY,T} =
    Mask(DoubleLayer(body,regop,w))

(m::Mask{N,NX,NY,G})(w::Nodes{G,NX,NY,T}) where {N,NX,NY,G,T} = m.data ∘ w

(m::ComplementaryMask{N,NX,NY,G})(w::Nodes{G,NX,NY,T}) where {N,NX,NY,G,T} = w - m.mask(w)

end
