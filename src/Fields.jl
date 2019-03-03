#=
Philosophy and convention:
A node is defined as the center of a grid cell. Grid cells may be those
in a primary grid or those in a dual grid.

The definitions "dual" and "primal" are made based on their conventional use
in holding fluid dynamic data. For example, pressure is held in primal nodes,
vorticity and streamfunction at dual nodes.

The definitions here are made with a "dual grid" in mind. That is, the grid is
defined by an integer number of dual cells in each direction. If a "primal grid"
is needed, then all of the defintions can be swapped (primal -> dual, dual -> primal).

Also, note that there might be dual cells that are "ghosts" (i.e. lie outside
the grid), but these are not distinguished in these basic definitions and operators.
=#


module Fields

import Base: @propagate_inbounds, show, summary

using Compat
using FFTW
using SpecialFunctions
using Compat.LinearAlgebra
using Compat.SparseArrays
using Compat: copyto!, reverse

@static if VERSION < v"0.7-"
  import Base: A_mul_B!, A_ldiv_B!, scale!, parentindexes
  import Base.LinAlg: cross, ×, dot, ⋅
  mul!(x,B,y) = A_mul_B!(x,B,y)
  rmul!(B,x) = scale!(B,x)
  ldiv!(x,B,y) = A_ldiv_B!(x,B,y)
  parentindices(A) = parentindexes(A)
  CartesianIndices(R) = CartesianRange(R)
  const GAMMA = γ
  export mul!, ldiv!, rmul!, parentindices
else
  import LinearAlgebra: mul!, ldiv!, cross, ×, dot, ⋅
  import Base: parentindices
  const GAMMA = MathConstants.γ
end

export Primal, Dual, Edges, Nodes,
       EdgeGradient, NodePair,
       Points, ScalarData, VectorData, TensorData,
       curl, curl!, Curl, divergence, divergence!, Divergence,
       grad, grad!, Grad,
       laplacian, laplacian!, plan_laplacian, plan_laplacian!,
       plan_intfact,plan_intfact!,Identity,
       product, product!, ∘,
       cellshift!,
       coordinates,
       DDF, GradDDF,
       Regularize, RegularizationMatrix, InterpolationMatrix,
       CircularConvolution

abstract type CellType end
abstract type Primal <: CellType end
abstract type Dual <: CellType end

macro wraparray(wrapper, field)
    T = supertype(eval(wrapper))
    @assert T <: AbstractArray "Wrapped type must be a subtype of AbstractArray"
    el_type, N = T.parameters

    quote
        Base.parent(A::$wrapper) = A.$field
        Base.size(A::$wrapper) = size(A.$field)
        parentindices(A::$wrapper) = parentindices(A.$field)

        if $N > 1
          function Base.show(io::IO, m::MIME"text/plain", A::$wrapper)
            println(io, "$(typeof(A)) data")
            println(io, "Printing in grid orientation (lower left is (1,1))")
            show(io,m, reverse(transpose(A.$field),dims=1))
          end
          #function Base.summary(io::IO, A::$wrapper)
          #  println(io, "$(typeof(A)) data")
          #  print(io, "Printing in grid orientation (lower left is (1,1))")
          #end
        end

        @propagate_inbounds Base.getindex(A::$wrapper, i::Int) = A.$field[i]
        @propagate_inbounds Base.setindex!(A::$wrapper, v, i::Int) = A.$field[i] = convert($el_type, v)
        if $N > 1
          @propagate_inbounds Base.getindex(A::$wrapper, I::Vararg{Int, $N}) = A.$field[I...]
          @propagate_inbounds Base.setindex!(A::$wrapper, v, I::Vararg{Int, $N}) = A.$field[I...] = convert($el_type, v)
        end
    end
end

function othertype end

macro othertype(celltype, k)
    esc(quote
        Fields.othertype(::$celltype) = $k
        Fields.othertype(::Type{$celltype}) = $k
    end)
end

@othertype Primal Dual
@othertype Dual Primal
@othertype CellType CellType

include("fields/nodes.jl")
include("fields/edges.jl")
include("fields/collections.jl")

scalarlist = ((:(Nodes{Primal,NX,NY}),1,1,0.0,0.0),
              (:(Nodes{Dual,NX,NY}),  0,0,0.5,0.5))

vectorlist = ((:(Edges{Primal,NX,NY}),          0,1,1,0,0.5,0.0,0.0,0.5),
              (:(Edges{Dual,NX,NY}),            1,0,0,1,0.0,0.5,0.5,0.0),
              (:(NodePair{Dual,Dual,NX,NY}),    0,0,0,0,0.5,0.5,0.5,0.5),
              (:(NodePair{Primal,Dual,NX,NY}),  1,1,0,0,0.0,0.0,0.5,0.5),
              (:(NodePair{Dual,Primal,NX,NY}),  0,0,1,1,0.5,0.5,0.0,0.0),
              (:(NodePair{Primal,Primal,NX,NY}),1,1,1,1,0.0,0.0,0.0,0.0))

tensorlist = ((:(EdgeGradient{Dual,Primal,NX,NY}), 0,0,1,1,0.5,0.5,0.0,0.0),
              (:(EdgeGradient{Primal,Dual,NX,NY}), 1,1,0,0,0.0,0.0,0.5,0.5))

include("fields/points.jl")

GridData = Union{Nodes{T,NX,NY},Edges{T,NX,NY}} where {T,NX,NY}

CollectedData = Union{EdgeGradient{T,S,NX,NY},NodePair{T,S,NX,NY}} where {T,S,NX,NY}

"""
    coordinates(w::Nodes/Edges;[dx=1.0],[I0=(1,1)])

Return a tuple of the ranges of the physical coordinates in each direction for grid
data `w`. If `w` is of `Nodes` type, then it returns a tuple of the form
`xg,yg`. If `w` is of `Edges` or `NodePair` type, then it returns a tuple of
the form `xgu,ygu,xgv,ygv`.

The optional keyword argument `dx` sets the grid spacing; its default is `1.0`. The
optional keyword `I0` accepts a tuple of integers to set the index pair of the
primal nodes that coincide with the origin. The default is `(1,1)`.

# Example

```jldoctest
julia> w = Nodes(Dual,(12,22));

julia> xg, yg = coordinates(w,dx=0.1)
(-0.05:0.1:1.05, -0.05:0.1:2.0500000000000003)
```
"""
function coordinates end

for (ctype,dnx,dny,shiftx,shifty) in scalarlist
   @eval coordinates(w::$ctype;dx::Float64=1.0,I0::Tuple{Int,Int}=(1,1)) where {NX,NY} =
    dx.*((1-I0[1]-$shiftx):(NX-$dnx-I0[1]-$shiftx),
         (1-I0[2]-$shifty):(NY-$dny-I0[2]-$shifty))

end

for (ctype,dunx,duny,dvnx,dvny,shiftux,shiftuy,shiftvx,shiftvy) in vectorlist
   @eval coordinates(w::$ctype;dx::Float64=1.0,I0::Tuple{Int,Int}=(1,1)) where {NX,NY} =
    dx.*((1-I0[1]-$shiftux):(NX-$dunx-I0[1]-$shiftux),
         (1-I0[2]-$shiftuy):(NY-$duny-I0[2]-$shiftuy),
         (1-I0[1]-$shiftvx):(NX-$dvnx-I0[1]-$shiftvx),
         (1-I0[2]-$shiftvy):(NY-$dvny-I0[2]-$shiftvy))


end

include("fields/physicalgrid.jl")
include("fields/operators.jl")
include("fields/shift.jl")



end
