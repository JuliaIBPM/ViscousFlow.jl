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

#=
Notes on data types:
All data types can be reduced to one of two different types in each direction.
(1) Lying on dual cell centers (C), numbered 1 through N
(2) Lying on dual cell edges (E), numbered 1 through N-1 (i.e. lying midway between C)

For example, Nodes{Dual} by definition lie on dual cell centers in both directions,
so they would be C x C. Nodes{Primal} are aligned with the corners of the dual
cells, so they lie along dual cell edges in both directions, E x E.

Nodes{Dual}    -- C x C
Nodes{Primal}  -- E x E
XEdges{Dual}   -- E x C
YEdges{Dual}   -- C x E
XEdges{Primal} -- C x E
YEdges{Primal} -- E x C

This is important when considering the interpolation and differentiation operations,
both of which must account for the relative indexing between E and C. The
key to remember is that E[1] is bounded by C[1] and C[2] in our indexing convention,
so that, for either interpolation or differentiation:
      E[i] <- C[i], C[i+1] for i in 1:N-1
and
      C[i] <- E[i-1], E[i] for i in 2:N-1
and
      C[i] <- C[i] for i in 1:N
and
      E[i] <- E[i] for i in 1:N-1
=#


module Fields

import Base: @propagate_inbounds, show, summary, fill!

#using Compat
using FFTW
using SpecialFunctions
using Statistics

using LinearAlgebra
using SparseArrays
using Interpolations

import LinearAlgebra: mul!, ldiv!, cross, ×, dot, ⋅

import Base: parentindices
const GAMMA = MathConstants.γ

export Primal, Dual, ScalarGridData, VectorGridData, GridData,
       Points, ScalarData, VectorData, TensorData,
       celltype, griddatatype,
       diff!,grid_interpolate!,
       curl, curl!, Curl, divergence, divergence!, Divergence,
       grad, grad!, Grad,
       laplacian, laplacian!, laplacian_symm!, plan_laplacian, plan_laplacian!,
       helmholtz, helmholtz!, plan_helmholtz, plan_helmholtz!,
       plan_intfact,plan_intfact!,Identity,
       product, product!, ∘,
       directional_derivative!, directional_derivative_conserve!, curl_cross!,
       convective_derivative!, convective_derivative_rot!,
       DDF, GradDDF,
       Regularize, RegularizationMatrix, InterpolationMatrix,
       CircularConvolution

abstract type CellType end
abstract type Primal <: CellType end
abstract type Dual <: CellType end

abstract type GridData{NX,NY,T} <: AbstractMatrix{T} end

abstract type ScalarGridData{NX,NY,T} <: GridData{NX,NY,T} end

abstract type VectorGridData{NX,NY,T} <: GridData{NX,NY,T} end

# List of scalar grid types. Each pair of numbers specifies
# the number of grid points in each direction for this data type, relative
# to the reference grid. The two pairs of numbers correspond to Primal
# and Dual versions of this grid data type.
const SCALARLIST = [ :Nodes, (-1,-1), (0,0)],
                   [ :XEdges, (0,-1), (-1,0)],
                   [ :YEdges, (-1,0), (0,-1)]

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

#@wraparray ScalarGridData data 2

include("fields/fieldmacros.jl")
include("fields/scalargrid.jl")

# Generate the scalar grid field types and associated functions
for (wrapper,primaldn,dualdn) in SCALARLIST
    @eval @scalarfield $wrapper $primaldn $dualdn
end

include("fields/collections.jl")


# The information in this list (after the cell types) follows directly from
# the known component types.
# (gtype,ctype(s),dunx,duny,dvnx,dvny,shiftux,shiftuy,shiftvx,shiftvy)
vectorlist = ((:Edges, [:Primal],      0,1,1,0,0.5,0.0,0.0,0.5),
              (:Edges, [:Dual],            1,0,0,1,0.0,0.5,0.5,0.0),
              (:NodePair, [:Primal,:Dual],  1,1,0,0,0.0,0.0,0.5,0.5),
              (:NodePair, [:Dual,:Primal],  0,0,1,1,0.5,0.5,0.0,0.0))

tensorlist = ((:EdgeGradient, [:Dual,:Primal], 0,0,1,1,0.5,0.5,0.0,0.0),
              (:EdgeGradient, [:Primal,:Dual], 1,1,0,0,0.0,0.0,0.5,0.5))

include("fields/basicoperations.jl")
include("fields/points.jl")

#CollectedData = Union{EdgeGradient{R,S,NX,NY,T},NodePair{R,S,NX,NY,T}} where {R,S,NX,NY,T}

include("fields/physicalgrid.jl")
include("fields/operators.jl")

end
