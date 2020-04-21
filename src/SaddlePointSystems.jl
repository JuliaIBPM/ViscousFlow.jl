module SaddlePointSystems

using LinearMaps
using RecursiveArrayTools

using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \

export SaddleSystem

struct SaddleSystem{T,Ns,Nc}
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
- `SaddleSystem(A::AbstractMatrix,B₂::AbstractMatrix,B₁ᵀ::AbstractMatrix,C::AbstractMatrix)`.
Blocks are given as matrices. Must have consistent sizes to stack appropriately.
- `SaddleSystem(A,B₂,B₁ᵀ,C,u,f)`.
Operators `A`, `B₂`, `B₁ᵀ`, `C` are given in various forms, including matrices, functions, and function-like objects.
`u` and `f` are examples of the data types in the corresponding solution and right-hand side vectors.
Guidelines:
* The entries `A` and `B₂` must be able to act upon `u` (either by multiplication or as a function) and
 `B₁ᵀ` and `C` must be able to act on `f` (also, either by multiplication or as a function).
* `A` and `B₁ᵀ` should return data of type `u`, and `B₂` and `C` should
 return data of type `f`.
* Both `u` and `f` must be mappable to/from `AbstractVector` type.
- `SaddleSystem(A,B₂,B₁ᵀ)` or `SaddleSystem(A,B₂,B₁ᵀ,u,f)`.
The `C` block is omitted and assumed to be zero.
"""
function SaddleSystem(A::LinearMap{T},B₂::LinearMap{T},B₁ᵀ::LinearMap{T},C::LinearMap{T},
                      A⁻¹::LinearMap{T}) where {T}

    ns, nc = _check_sizes(A,B₂,B₁ᵀ,C)

    S = C - B₂*A⁻¹*B₁ᵀ

    Sfact = factorize(Matrix(S))
    S⁻¹ = LinearMap{T}(x -> Sfact\x,nc)

    return SaddleSystem{T,ns,nc}(A,B₂,B₁ᵀ,C,A⁻¹,zeros(T,ns),zeros(T,nc),zeros(T,ns),zeros(T,nc),S,S⁻¹)
end


function SaddleSystem(A::AbstractMatrix{T},B₂::AbstractMatrix{T},B₁ᵀ::AbstractMatrix{T},
                      C::AbstractMatrix{T}) where {T}

    Afact = factorize(A)
    Ainv = LinearMap{T}(x -> Afact\x,size(A,1))

    return SaddleSystem(LinearMap{T}(A),LinearMap{T}(B₂),LinearMap{T}(B₁ᵀ),
                        LinearMap{T}(C),Ainv)
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
                        linear_inverse_map(A,u,eltype=eltype))
end

SaddleSystem(A,B₂,B₁ᵀ,u::TU,f::TF;eltype=Float64) where {TU,TF} =
    SaddleSystem(A,B₂,B₁ᵀ,zeros(eltype,length(f),length(f)),u,f,eltype=eltype)

### AUXILIARY ROUTINES

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


include("saddlepoint/linearmaps.jl")
include("saddlepoint/arithmetic.jl")


# """
#     SaddleSystem((u,f),(A⁻¹,B₁ᵀ,B₂);[tol=1e-3],[issymmetric=false],[isposdef=false],[conditioner=Identity],[store=false])
#
# Construct the computational operators for a saddle-point system of the form
# \$[A B₁ᵀ; B₂ 0][u;f]\$. Note that the constituent operators are passed in as a
# tuple in the order seen here. Each of these operators could act on its corresponding
# data type in a function-like way, e.g. `A⁻¹(u)`, or in a matrix-like way, e.g.,
# `A⁻¹*u`.
#
# The optional argument `tol` sets the tolerance for iterative solution (if
#   applicable). Its default is 1e-3.
#
# The optional argument `conditioner` can be used to supply a function that acts
# upon the result `f` to 'condition' it (e.g. filter it). It is, by default, set
# to the identity.
#
# The optional Boolean argument `store` will compute and store the Schur
# complement matrix's factorization. This makes the inversion faster, though
# it comes at the expense of memory and overhead time for pre-computing it. The
# resulting solution is somewhat noiser, too.
#
# # Arguments
#
# - `u` : example of state vector data.
# - `f` : example of constraint force vector data. This data must be of
#         AbstractVector supertype.
# - `A⁻¹` : operator evaluating the inverse of `A` on data of type `u`, return type `u`
# - `B₁ᵀ` : operator evaluating the influence of constraint force,
#             acting on `f` and returning type `u`
# - `B₂` : operator evaluating the influence of state vector on constraints,
#             acting on `u` and returning type `f`
# """
# struct SaddleSystem{FA,FAB,FBA,FP,TU,TF,T,N,Storage}
#     A⁻¹B₁ᵀf :: TU
#     B₂A⁻¹r₁ :: TF
#     tmpvec :: Vector{T}
#     tmpvecout :: Vector{T}
#     A⁻¹ :: FA
#     A⁻¹B₁ᵀ :: FAB
#     B₂A⁻¹ :: FBA
#     P :: FP
#     S  :: LinearMap
#     S⁻¹ :: Union{Factorization{T},Nothing}
#     tol :: Float64
#     _issymmetric :: Bool
#     _isposdef :: Bool
# end
#
# function (::Type{SaddleSystem})(state::Tuple{TU,TF},sys::Tuple{FA,FB1,FB2};
#                                 tol::Float64=1e-3,
#                                 conditioner::FP=x->x,
#                                 issymmetric::Bool=false,
#                                 isposdef::Bool=false,
#                                 store::Bool=false,
#                                 precompile::Bool=true) where {TU,TF,FA,FB1,FB2,FP}
#     u,f = state
#
#     optypes = (TU,TF,TU)
#     opnames = ("A⁻¹","B₁ᵀ","B₂")
#     ops = []
#     T = eltype(TU)
#
#     # check for methods
#     for (i,typ) in enumerate(optypes)
#       if hasmethod(sys[i],Tuple{typ})
#         push!(ops,sys[i])
#       elseif hasmethod(*,Tuple{typeof(sys[i]),typ})
#         # generate a method that acts on TU
#         push!(ops,x->sys[i]*x)
#       else
#         error("No valid operator for $(opnames[i]) supplied")
#       end
#     end
#
#     A⁻¹, B₁ᵀ, B₂ = ops
#     ubuffer = deepcopy(u)
#     fbuffer = deepcopy(f)
#     N = length(f)
#     tmpvec = zeros(T,N)
#     tmpvecout = zeros(T,N)
#
#     # Schur complement
#     function Schur!(fout::AbstractVector{T},fin::AbstractVector{T}) where {T<:Number}
#        fbuffer .= fin
#        fbuffer .= (B₂∘A⁻¹∘B₁ᵀ)(fbuffer)
#        fout .= -fbuffer
#        return fout
#     end
#     S = LinearMap{T}(Schur!,N;ismutating=true,issymmetric=issymmetric,isposdef=isposdef)
#
#     if store && N > 0
#       Smat = zeros(T,N,N)
#       fill!(f,0.0)
#       for i = 1:N
#         f[i] = 1.0
#         fbuffer .= S*f
#         Smat[1:N,i] .= fbuffer
#         f[i] = 0.0
#       end
#       #S⁻¹ = Nullable(factorize(Smat))
#       S⁻¹ = factorize(Smat)
#     else
#       #S⁻¹ = Nullable()
#       S⁻¹ = nothing
#     end
#
#     A⁻¹B₁ᵀ(f::TF) = (A⁻¹∘B₁ᵀ)(f)
#
#     B₂A⁻¹(w::TU) = (B₂∘A⁻¹)(w)
#
#     saddlesys = SaddleSystem{typeof(A⁻¹),typeof(A⁻¹B₁ᵀ),typeof(B₂A⁻¹),typeof(conditioner),TU,TF,T,N,store}(
#                                 ubuffer,fbuffer,tmpvec,tmpvecout,
#                                 A⁻¹,A⁻¹B₁ᵀ,B₂A⁻¹,conditioner,S,S⁻¹,
#                                 tol,issymmetric,isposdef)
#     # run once in order to precompile it
#     if precompile
#       saddlesys\(u,f)
#     end
#
#     return saddlesys
#
# end
#
# (::Type{SaddleSystem})(state::Tuple{TU,TF},sys::Array{Any,2};kwarg...) where {TU,TF} =
#             SaddleSystem(state,(sys[1,1],sys[1,2],sys[2,1]);kwarg...)
#
# function Base.show(io::IO, S::SaddleSystem{FA,FAB,FBA,FP,TU,TF,T,N,Storage}) where {FA,FAB,FBA,FP,TU,TF,T,N,Storage}
#     println(io, "Saddle system with $N constraints and")
#     println(io, "   State of type $TU")
#     println(io, "   Force of type $TF")
# end
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
# # stored matrix
# function ldiv!(state::Tuple{TU,TF},
#                     sys::SaddleSystem{FA,FAB,FBA,FP,TU,TF,T,N,true},
#                     rhs::Tuple{TU,TF}) where {TU,TF,T,FA,FAB,FBA,FP,N}
#
#   ru, rf = rhs
#   u, f = state
#   sys.B₂A⁻¹r₁ .= sys.B₂A⁻¹(ru)
#   sys.tmpvec .= rf
#   sys.tmpvec .-= sys.B₂A⁻¹r₁
#   if N > 0
#     ldiv!(sys.S⁻¹,sys.tmpvec)
#     f .= sys.tmpvec
#     f .= sys.P(f)
#   end
#   u .= sys.A⁻¹(ru)
#   sys.A⁻¹B₁ᵀf .= sys.A⁻¹B₁ᵀ(f)
#   u .-= sys.A⁻¹B₁ᵀf
#   state = u, f
# end
#
#
# \(sys::SaddleSystem{FA,FAB,FBA,FP,TU,TF,T,N,Storage},rhs::Tuple{TU,TF}) where {TU,TF,T,FA,FAB,FBA,FP,N,Storage} =
#       ldiv!(similar.(rhs),sys,rhs)
#
#
# # solving tuples of systems
# function ldiv!(state,sys::NTuple{M,SaddleSystem},rhs) where {M}
#   for (i,sysi) in enumerate(sys)
#     ldiv!(state[i],sysi,rhs[i])
#   end
#   state
# end
#
# \(sys::NTuple{M,SaddleSystem},rhs) where {M} =
#       ldiv!(deepcopy.(rhs),sys,rhs)

end
