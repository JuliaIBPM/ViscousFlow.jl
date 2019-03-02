module Systems

using DocStringExtensions

export NavierStokes, PointForce

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions
using ..Bodies
#import ViscousFlow: r₁
#import ViscousFlow: plan_intfact

using Compat
using Compat.LinearAlgebra
using Compat.SparseArrays


include("systems/navier_stokes.jl")

end
