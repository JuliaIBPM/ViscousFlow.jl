module SaddlePointSystems

using LinearMaps
using IterativeSolvers
using ..Fields

import Base: *, \, A_mul_B!, A_ldiv_B!

export SaddleSystem

struct SaddleSystem{FA,FAB,FBA,TU,TF,N}
    A⁻¹B₁ᵀf :: TU
    B₂A⁻¹r₁ :: TF
    A⁻¹! :: FA
    A⁻¹B₁ᵀ! :: FAB
    B₂A⁻¹! :: FBA
    S  :: LinearMap
    _issymmetric :: Bool
    _isposdef :: Bool
end

"""
    SaddleSystem(u,f,A⁻¹!,B₁ᵀ!,B₂!;[issymmetric=false],[isposdef=false])

Construct a saddle-point system.

# Arguments

- `u` : example of state vector data
- `f` : example of constraint force vector data
- `A⁻¹!` : operator evaluating the inverse of `A`
- `B₁ᵀ!` : operator evaluating the influence of constraint force on system
- `B₂!` : operator evaluating the influence of state vector on constraints
"""
function (::Type{SaddleSystem})(u::TU,f::TF,A⁻¹!::FA,B₁ᵀ!::FB1,B₂!::FB2;
                                issymmetric::Bool=false,
                                isposdef::Bool=false) where {TU,TF,FA,FB1,FB2}
    ubuffer = deepcopy(u)
    fbuffer = deepcopy(f)
    N = length(f)

    # Schur complement
    function Schur!(fout::AbstractVector{Float64},fin::AbstractVector{Float64})
       fbuffer .= fin
       B₂!(fbuffer,A⁻¹!(ubuffer,B₁ᵀ!(ubuffer,fbuffer)))
       fout .= -fbuffer
       return fout
    end
    S = LinearMap(Schur!,N;ismutating=true,issymmetric=issymmetric,isposdef=isposdef)

    function A⁻¹B₁ᵀ!(w::TU,f::TF)
        A⁻¹!(w,B₁ᵀ!(ubuffer,f))
        return w
    end

    function B₂A⁻¹!(f::TF,w::TU)
        B₂!(f,A⁻¹!(ubuffer,w))
        return f
    end

    SaddleSystem{FA,typeof(A⁻¹B₁ᵀ!),typeof(B₂A⁻¹!),TU,TF,N}(ubuffer,fbuffer,
                                A⁻¹!,A⁻¹B₁ᵀ!,B₂A⁻¹!,S,
                                issymmetric,isposdef)
end

function Base.show(io::IO, S::SaddleSystem{FA,FAB,FBA,TU,TF,N}) where {FA,FAB,FBA,TU,TF,N}
    print(io, "Saddle system with state of type $TU and force of type $TF")
end


"""
    A_ldiv_B!(state,sys::SaddleSystem,rhs)

Solve a saddle-point system. `rhs` is a tuple of the right-hand side (ru,rf).
Output `state`, a tuple (u,f), is updated. Note that `sys` is also mutated:
its scratch space `sys.B₂A⁻¹r₁` and `sys.A⁻¹B₁ᵀf` hold the intermediate results
of the solution.
"""
function A_ldiv_B!(state::Tuple{TU,TF},
                    sys::SaddleSystem{FA,FAB,FBA,TU,TF,N},
                    rhs::Tuple{TU,TF}) where {TU,TF,FA,FAB,FBA,N}

  ru, rf = rhs
  u, f = state
  rf .-= sys.B₂A⁻¹!(sys.B₂A⁻¹r₁,ru)
  cg!(f,sys.S,rf,tol=1e-3)
  sys.A⁻¹!(u,ru)
  u .-= sys.A⁻¹B₁ᵀ!(sys.A⁻¹B₁ᵀf,f)
  state = u, f
end

end
