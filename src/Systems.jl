module Systems

using DocStringExtensions
using SparseArrays

export NavierStokes, PointForce

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions
using ..Bodies
#import ViscousFlow: r₁
#import ViscousFlow: plan_intfact

using Compat
using Compat.LinearAlgebra

include("systems/navier_stokes.jl")

end
