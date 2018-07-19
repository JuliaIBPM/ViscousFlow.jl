module SaddlePointSystems

using LinearMaps
using IterativeSolvers
using ..Fields

import Base: *, \, A_mul_B!, A_ldiv_B!

export SaddleSystem

struct SaddleSystem{FA,FAB,FBA,FP,TU,TF,N}
    A⁻¹B₁ᵀf :: TU
    B₂A⁻¹r₁ :: TF
    A⁻¹ :: FA
    A⁻¹B₁ᵀ :: FAB
    B₂A⁻¹ :: FBA
    P :: FP
    S  :: LinearMap
    _issymmetric :: Bool
    _isposdef :: Bool
end

"""
    SaddleSystem((u,f),(A⁻¹,B₁ᵀ,B₂);[issymmetric=false],[isposdef=false])

Construct the computational operators for a saddle-point system of the form
\$[A B₁ᵀ; B₂ 0][u;f]\$. Note that the constituent operators are passed in as a
tuple in the order seen here. Each of these operators could act on its corresponding
data type in a function-like way, e.g. A⁻¹(u), or in a matrix-like way, e.g.,
A⁻¹*u.

# Arguments

- `u` : example of state vector data.
- `f` : example of constraint force vector data. This data must be of
        AbstractVector supertype.
- `A⁻¹` : operator evaluating the inverse of `A` on data of type `u`, return type `u`
- `B₁ᵀ` : operator evaluating the influence of constraint force,
            acting on `f` and returning type `u`
- `B₂` : operator evaluating the influence of state vector on constraints,
            acting on `u` and returning type `f`
"""
function (::Type{SaddleSystem})(state::Tuple{TU,TF},sys::Tuple{FA,FB1,FB2};
                                conditioner::FP=x->x,
                                issymmetric::Bool=false,
                                isposdef::Bool=false) where {TU,TF,FA,FB1,FB2,FP}
    u,f = state

    optypes = (TU,TF,TU)
    opnames = ("A⁻¹","B₁ᵀ","B₂")
    ops = []

    # check for methods
    for (i,typ) in enumerate(optypes)
      if method_exists(sys[i],Tuple{typ})
        push!(ops,sys[i])
      elseif method_exists(*,Tuple{typeof(sys[i]),typ})
        # generate a method that acts on TU
        push!(ops,x->sys[i]*x)
      else
        error("No valid operator for $(opnames[i]) supplied")
      end
    end

    A⁻¹, B₁ᵀ, B₂ = ops
    ubuffer = deepcopy(u)
    fbuffer = deepcopy(f)
    N = length(f)

    # Schur complement
    function Schur!(fout::AbstractVector{Float64},fin::AbstractVector{Float64})
       fbuffer .= fin
       fbuffer .= (B₂∘A⁻¹∘B₁ᵀ)(fbuffer)
       fout .= -fbuffer
       return fout
    end
    S = LinearMap(Schur!,N;ismutating=true,issymmetric=issymmetric,isposdef=isposdef)

    A⁻¹B₁ᵀ(f::TF) = (A⁻¹∘B₁ᵀ)(f)

    B₂A⁻¹(w::TU) = (B₂∘A⁻¹)(w)

    saddlesys = SaddleSystem{typeof(A⁻¹),typeof(A⁻¹B₁ᵀ),typeof(B₂A⁻¹),typeof(conditioner),TU,TF,N}(ubuffer,fbuffer,
                                A⁻¹,A⁻¹B₁ᵀ,B₂A⁻¹,conditioner,S,
                                issymmetric,isposdef)
    # run once in order to precompile it
    saddlesys\(u,f)

    return saddlesys

end

(::Type{SaddleSystem})(state::Tuple{TU,TF},sys::Array{Any,2};kwarg...) where {TU,TF} =
            SaddleSystem(state,(sys[1,1],sys[1,2],sys[2,1]);kwarg...)

function Base.show(io::IO, S::SaddleSystem{FA,FAB,FBA,FP,TU,TF,N}) where {FA,FAB,FBA,FP,TU,TF,N}
    println(io, "Saddle system with")
    println(io, "   State of type $TU")
    println(io, "   Force of type $TF")
end


"""
    A_ldiv_B!(state,sys::SaddleSystem,rhs)

Solve a saddle-point system. `rhs` is a tuple of the right-hand side `(ru,rf)`.
Output `state`, a tuple (u,f), is updated. Note that `sys` is also mutated:
its scratch space `sys.B₂A⁻¹r₁` and `sys.A⁻¹B₁ᵀf` hold the intermediate results
of the solution.

A shorthand can be used for this operation: `state = sys\rhs`
"""
function A_ldiv_B!(state::Tuple{TU,TF},
                    sys::SaddleSystem{FA,FAB,FBA,FP,TU,TF,N},
                    rhs::Tuple{TU,TF}) where {TU,TF,FA,FAB,FBA,FP,N}

  ru, rf = rhs
  u, f = state
  sys.B₂A⁻¹r₁ .= sys.B₂A⁻¹(ru)
  rf .-= sys.B₂A⁻¹r₁
  cg!(f,sys.S,rf,tol=1e-3)
  f .= sys.P(f)
  u .= sys.A⁻¹(ru)
  sys.A⁻¹B₁ᵀf .= sys.A⁻¹B₁ᵀ(f)
  u .-= sys.A⁻¹B₁ᵀf
  state = u, f
end

\(sys::SaddleSystem{FA,FAB,FBA,FP,TU,TF,N},rhs::Tuple{TU,TF}) where {TU,TF,FA,FAB,FBA,FP,N} =
      A_ldiv_B!((TU(),TF()),sys,rhs)

end
