module SaddlePointSystems

#=
Still to do:
* Recursive saddle point systems
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


end
