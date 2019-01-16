module Systems

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions
using ..Bodies
#import Whirl: r₁
#import Whirl: plan_intfact

using Compat
using Compat.LinearAlgebra



include("systems/navier_stokes.jl")

end
