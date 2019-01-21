module Systems

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions
using ..Bodies
#import ViscousFlow: r₁
#import ViscousFlow: plan_intfact

using Compat
using Compat.LinearAlgebra

export NavierStokes, PointForce

include("systems/navier_stokes.jl")

end
