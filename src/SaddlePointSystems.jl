module SaddlePointSystems

#=
Still to do:
* Need tests of vector force data
* Recursive saddle point systems
* Iterative (non-stored) Schur complement solution
=#

using LinearMaps
using RecursiveArrayTools

using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \

import Base: size, eltype

export SaddleSystem, SaddleVector, state, constraint

struct SaddleSystem{T,Ns,Nc,TF,TU}
    A :: LinearMap{T}
    B₂ :: LinearMap{T}
    B₁ᵀ :: LinearMap{T}
    C :: LinearMap{T}
    A⁻¹ :: LinearMap{T}
    A⁻¹B₁ᵀf :: Vector{T}
    B₂A⁻¹r₁ :: Vector{T}
    _u_buf :: Vector{T}
    _f_buf :: Vector{T}
    S :: LinearMap{T}
    S⁻¹ :: LinearMap{T}
end


### CONSTRUCTORS ###

"""
    SaddleSystem

Construct a saddle-point system operator from the constituent operator blocks. The resulting object can be used
with `*` and `\` to multiply and solve. The saddle-point problem has the form

``
\\begin{bmatrix}A & B_1^T \\\\ B_2 & C \\end{bmatrix} \\begin{pmatrix} u \\\\ f \\end{pmatrix} = \\begin{pmatrix} r_1 \\\\ r_2 \\end{pmatrix}
``

### Constructors
`SaddleSystem(A::AbstractMatrix,B₂::AbstractMatrix,B₁ᵀ::AbstractMatrix,C::AbstractMatrix[,eltype=Float64])`.
Blocks are given as matrices. Must have consistent sizes to stack appropriately. If
this is called with `SaddleSystem(A,B₂,B₁ᵀ)`, it sets `C` to zero automatically.

`SaddleSystem(A,B₂,B₁ᵀ,C,u,f[,eltype=Float64])`.
Operators `A`, `B₂`, `B₁ᵀ`, `C` are given in various forms, including matrices, functions, and function-like objects.
`u` and `f` are examples of the data types in the corresponding solution and right-hand side vectors.
Guidelines:

* The entries `A` and `B₂` must be able to act upon `u` (either by multiplication or as a function) and `B₁ᵀ` and `C` must be able to act on `f` (also, either by multiplication or as a function).
* `A` and `B₁ᵀ` should return data of type `u`, and `B₂` and `C` should return data of type `f`.
* `A` must be invertible and be outfitted with operators `\` and `ldiv!`.
* Both `u` and `f` must be subtypes of `AbstractArray`: they must be equipped with `size`
  and `vec` functions and with a constructor of the form `T(data)` where `T` is the data type of
  `u` or `f` and `data` is the wrapped data array.

If called as `SaddleSystem(A,B₂,B₁ᵀ,u,f)`, the `C` block is omitted and assumed to be zero.

If called with `SaddleSystem(A,u)`, this is equivalent to calling `SaddleSystem(A,nothing,nothing,u,[])`, then this reverts
to the unconstrained system described by operator `A`.

The list of vectors in any of these constructors can be replaced by a `SaddleVector`,
e.g. `SaddleSystem(A,B₂,B₁ᵀ,SaddleVector(u,f))`.
"""
function SaddleSystem(A::LinearMap{T},B₂::LinearMap{T},B₁ᵀ::LinearMap{T},C::LinearMap{T},
                      A⁻¹::LinearMap{T},TU,TF) where {T}

    ns, nc = _check_sizes(A,B₂,B₁ᵀ,C)

    S = C - B₂*A⁻¹*B₁ᵀ

    Sfact = factorize(Matrix(S))
    S⁻¹ = LinearMap{T}(x -> Sfact\x,nc)

    return SaddleSystem{T,ns,nc,TU,TF}(A,B₂,B₁ᵀ,C,A⁻¹,zeros(T,ns),zeros(T,nc),zeros(T,ns),zeros(T,nc),S,S⁻¹)
end


function SaddleSystem(A::AbstractMatrix{T},B₂::AbstractMatrix{T},B₁ᵀ::AbstractMatrix{T},
                      C::AbstractMatrix{T}) where {T}

    Afact = factorize(A)
    Ainv = LinearMap{T}(x -> Afact\x,size(A,1))

    return SaddleSystem(LinearMap{T}(A),LinearMap{T}(B₂),LinearMap{T}(B₁ᵀ),
                        LinearMap{T}(C),Ainv,Vector{T},Vector{T})
end

# For cases in which C is zero, no need to pass along the argument
SaddleSystem(A::AbstractMatrix{T},B₂::AbstractMatrix{T},B₁ᵀ::AbstractMatrix{T}) where {T} =
        SaddleSystem(A,B₂,B₁ᵀ,zeros(T,size(B₂,1),size(B₁ᵀ,2)))

# This version should take in functions or function-like objects that act upon given
# data types u and f. Should transform them into operators that act on abstract vectors
# of the same size
# There should already be an \ operator associated with A
function SaddleSystem(A,B₂,B₁ᵀ,C,u::TU,f::TF;eltype=Float64) where {TU,TF}

    return SaddleSystem(linear_map(A,u,eltype=eltype),linear_map(B₂,u,f,eltype=eltype),
                        linear_map(B₁ᵀ,f,u,eltype=eltype),
                        linear_map(C,f,eltype=eltype),
                        linear_inverse_map(A,u,eltype=eltype),TU,TF)
end

SaddleSystem(A,B₂,B₁ᵀ,u::TU,f::TF;eltype=Float64) where {TU,TF} =
    SaddleSystem(A,B₂,B₁ᵀ,C_zero(f,eltype),u,f,eltype=eltype)


SaddleSystem(A,u::TU;eltype=Float64) where {TU} = SaddleSystem(A,nothing,nothing,u,Type{eltype}[])

SaddleSystem(A,B₂,B₁ᵀ,C,v::ArrayPartition;eltype=Float64) = SaddleSystem(A,B₂,B₁ᵀ,C,v.x[1],v.x[2],eltype=eltype)
SaddleSystem(A,B₂,B₁ᵀ,v::ArrayPartition;eltype=Float64) = SaddleSystem(A,B₂,B₁ᵀ,v.x[1],v.x[2],eltype=eltype)
SaddleSystem(A,v::ArrayPartition;eltype=Float64) = SaddleSystem(A,v.x[1],eltype=eltype)


function Base.show(io::IO, S::SaddleSystem{T,Ns,Nc,TU,TF}) where {T,Ns,Nc,TU,TF}
     println(io, "Saddle system with $Ns states and $Nc constraints and")
     println(io, "   State vector of type $TU")
     println(io, "   Constraint vector of type $TF")
     println(io, "   Elements of type $T")
end

### AUXILIARY ROUTINES

size(::SaddleSystem{T,Ns,Nc}) where {T,Ns,Nc} = (Ns+Nc,Ns+Nc)

eltype(::SaddleSystem{T,Ns,Nc}) where {T,Ns,Nc} = T

C_zero(f,eltype) = zeros(eltype,length(f),length(f))

function _check_sizes(A,B₂,B₁ᵀ,C)
    mA, nA = size(A)
    mB1, nB1 = size(B₁ᵀ)
    mB2, nB2 = size(B₂)
    mC, nC = size(C)

    # check compatibility of sizes
    mA == nA  || error("A is not square")
    mA == mB1 || error("Incompatible number of rows in A and B₁ᵀ")
    nA == nB2 || error("Incompatible number of columns in A and B₂")
    mC == mB2 || error("Incompatible number of rows in C and B₂")
    nC == nB1 || error("Incompatible number of columns in C and B₁ᵀ")

    ns = nA
    nc = nB1

    return ns, nc
end

include("saddlepoint/vectors.jl")
include("saddlepoint/linearmaps.jl")
include("saddlepoint/arithmetic.jl")



#
#
# """
#     ldiv!(state,sys::SaddleSystem,rhs)
#
# Solve a saddle-point system. `rhs` is a tuple of the right-hand side `(ru,rf)`.
# Output `state`, a tuple (u,f), is updated. Note that `sys` is also mutated:
# its scratch space `sys.B₂A⁻¹r₁` and `sys.A⁻¹B₁ᵀf` hold the intermediate results
# of the solution.
#
# A shorthand can be used for this operation: `state = sys\\rhs`
# """
# function ldiv!(state::Tuple{TU,TF},
#                     sys::SaddleSystem{FA,FAB,FBA,FP,TU,TF,T,N,false},
#                     rhs::Tuple{TU,TF}) where {TU,TF,T,FA,FAB,FBA,FP,N}
#   # non-stored matrix
#
#   ru, rf = rhs
#   u, f = state
#   sys.B₂A⁻¹r₁ .= sys.B₂A⁻¹(ru)
#   sys.tmpvec .= rf
#   sys.tmpvec .-= sys.B₂A⁻¹r₁
#   if N > 0
#     sys.tmpvecout .= f
#     cg!(sys.tmpvecout,sys.S,sys.tmpvec,tol=sys.tol)
#     f .= sys.tmpvecout
#     f .= sys.P(f)
#   end
#   u .= sys.A⁻¹(ru)
#   sys.A⁻¹B₁ᵀf .= sys.A⁻¹B₁ᵀ(f)
#   u .-= sys.A⁻¹B₁ᵀf
#   state = u, f
# end
#



end
