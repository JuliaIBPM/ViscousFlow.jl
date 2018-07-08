module SaddlePointSystems

using LinearMaps
using IterativeSolvers
using ..Fields

import Base: *, \, A_mul_B!, A_ldiv_B!

export SaddleSystem

struct SaddleSystem{FA,FAB,FBA,TU,TF,N}
    A⁻¹B₁ᵀf :: TU
    B₂A⁻¹r₁ :: TF
    A⁻¹ :: FA
    A⁻¹B₁ᵀ :: FAB
    B₂A⁻¹ :: FBA
    S  :: LinearMap
    _issymmetric :: Bool
    _isposdef :: Bool
end

"""
    SaddleSystem(u,f,A⁻¹,B₁ᵀ,B₂;[issymmetric=false],[isposdef=false])

Construct a saddle-point system.

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
function (::Type{SaddleSystem})(u::TU,f::TF,A⁻¹::FA,B₁ᵀ::FB1,B₂::FB2;
                                issymmetric::Bool=false,
                                isposdef::Bool=false) where {TU,TF,FA,FB1,FB2}
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

    SaddleSystem{FA,typeof(A⁻¹B₁ᵀ),typeof(B₂A⁻¹),TU,TF,N}(ubuffer,fbuffer,
                                A⁻¹,A⁻¹B₁ᵀ,B₂A⁻¹,S,
                                issymmetric,isposdef)
end

function Base.show(io::IO, S::SaddleSystem{FA,FAB,FBA,TU,TF,N}) where {FA,FAB,FBA,TU,TF,N}
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

A shorthand can be used for this operation: state = sys\rhs
"""
function A_ldiv_B!(state::Tuple{TU,TF},
                    sys::SaddleSystem{FA,FAB,FBA,TU,TF,N},
                    rhs::Tuple{TU,TF}) where {TU,TF,FA,FAB,FBA,N}

  ru, rf = rhs
  u, f = state
  sys.B₂A⁻¹r₁ .= sys.B₂A⁻¹(ru)
  rf .-= sys.B₂A⁻¹r₁
  cg!(f,sys.S,rf,tol=1e-3)
  u .= sys.A⁻¹(ru)
  sys.A⁻¹B₁ᵀf .= sys.A⁻¹B₁ᵀ(f)
  u .-= sys.A⁻¹B₁ᵀf
  state = u, f
end

\(sys::SaddleSystem{FA,FAB,FBA,TU,TF,N},rhs::Tuple{TU,TF}) where {TU,TF,FA,FAB,FBA,N} =
      A_ldiv_B!((TU(),TF()),sys,rhs)

end
