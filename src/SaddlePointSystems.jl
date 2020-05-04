module SaddlePointSystems

#=
Still to do:
* Recursive saddle point systems
* Iterative (non-stored) Schur complement solution
=#

using LinearMaps
using RecursiveArrayTools
using IterativeSolvers
using UnPack

using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \, I

import Base: size, eltype

export SaddleSystem, SaddleVector, state, constraint
export SchurSolverType, Direct, CG, BiCG, GMRES

include("saddlepoint/saddlesystems.jl")
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
