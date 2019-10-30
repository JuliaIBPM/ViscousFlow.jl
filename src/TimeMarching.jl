module TimeMarching

import ViscousFlow:@get

#using Compat
#using Compat: hasmethod

using LinearAlgebra

import LinearAlgebra: ldiv!

export System, Constrained, Unconstrained, RK, IFRK, IFHERK, r₁, r₂, B₂, B₁ᵀ,
          plan_constraints

"Abstract type for a system of ODEs"
abstract type System{C} end

const Constrained = true
const Unconstrained = false

# Functions that get extended by individual systems
function r₁ end
function r₂ end
function B₂ end
function B₁ᵀ end
function plan_constraints end

struct RKParams{N}
  c::Vector{Float64}
  a::Matrix{Float64}
end

using ..SaddlePointSystems

const RK31 = RKParams{3}([0.5, 1.0, 1.0],
                      [1/2        0        0
                       √3/3 (3-√3)/3        0
                       (3+√3)/6    -√3/3 (3+√3)/6])

const Euler = RKParams{1}([1.0],ones(1,1))

include("timemarching/rk.jl")
include("timemarching/ifrk.jl")
include("timemarching/ifherk.jl")

end
