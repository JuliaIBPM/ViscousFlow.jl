#=
 Methods and types associated with single and double layer potentials
=#

module Layers

using ..Fields
using ..Bodies

export DoubleLayer, SingleLayer

abstract type LayerType{N,NX,NY} end

struct DoubleLayer{N,NX,NY,G,DT} <: LayerType{N,NX,NY}
    nds :: VectorData{N,Float64}
    H :: RegularizationMatrix{Edges{G,NX,NY,DT},VectorData{N,DT}}
end

function DoubleLayer(body::Body,H::RegularizationMatrix)
  ds = ScalarData(Bodies.dlength(body))
  normals = VectorData(Bodies.normal(body))
  return DoubleLayer(normals∘ds,H)
end

function DoubleLayer(a,reg::Regularize{N},w::GridData{NX,NY,T}) where {N,NX,NY,T}

  out = RegularizationMatrix(reg,VectorData(N,dtype=T),Edges(celltype(w),w,dtype=T))
  return DoubleLayer(a,reg._issymmetric ? out[1] : out)
end

(μ::DoubleLayer{N})(p::ScalarData{N}) where {N} = divergence(μ.H*(p∘μ.nds))

function Base.show(io::IO, H::DoubleLayer{N,NX,NY,G,DT}) where {N,NX,NY,G,DT}
    println(io, "Double-layer potential mapping")
    println(io, "  from $N scalar-valued point data of $DT type")
    println(io, "  to a $NX x $NY grid of $G nodal data")
end

struct SingleLayer{N,NX,NY,G,DT} <: LayerType{N,NX,NY}
    ds :: ScalarData{N,Float64}
    H :: RegularizationMatrix{Nodes{G,NX,NY,DT},ScalarData{N,DT}}
end

function SingleLayer(body::Body,H::RegularizationMatrix)
  ds = ScalarData(Bodies.dlength(body))
  return SingleLayer(ds,H)
end

function SingleLayer(a,reg::Regularize{N},w::GridData{NX,NY,T}) where {N,NX,NY,T}

  out = RegularizationMatrix(reg,ScalarData(N,dtype=T),Nodes(celltype(w),w,dtype=T))
  return SingleLayer(a,reg._issymmetric ? out[1] : out)
end

(μ::SingleLayer{N})(p::ScalarData{N}) where {N} = μ.H*(p∘μ.ds)

function Base.show(io::IO, H::SingleLayer{N,NX,NY,G,DT}) where {N,NX,NY,G,DT}
    println(io, "Single-layer potential mapping")
    println(io, "  from $N scalar-valued point data of $DT type")
    println(io, "  to a $NX x $NY grid of $G nodal data")
end


end
