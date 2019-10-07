module Systems

using DocStringExtensions

export NavierStokes, PointForce, vorticity, velocity, streamfunction, nl

using ..Fields
using ..TimeMarching
using ..RigidBodyMotions
using ..Bodies
#import ViscousFlow: r₁
#import ViscousFlow: plan_intfact

using Compat
using Compat.LinearAlgebra
using Compat.SparseArrays

"""
    assign_velocity!(V::VectorData,X::VectorData,
                     xc::Real,yc::Real,α::Real,
                     mlist::Vector{RigidBodyMotion},t::Real)

Assign the components of rigid body velocity for every body (in inertial coordinate system)
at time `t` in the overall data structure `V`, using coordinates described by `X` (also in inertial
coordinate system), based on array of supplied motion `mlist` for each body.
"""
function RigidBodyMotions.assign_velocity!(V::VectorData{N},X::VectorData{N},
                                           bl::BodyList,tlist::Vector{RigidTransform},
                                           mlist::Vector{RigidBodyMotion},t::Real) where {N}
    N == numpts(bl) || error("Inconsistent size of data structures")
    for i in 1:length(bl)
        ui = view(V.u,bl,i)
        vi = view(V.v,bl,i)
        xi = view(X.u,bl,i)
        yi = view(X.v,bl,i)
        Ti = tlist[i]
        assign_velocity!(ui,vi,xi,yi,Ti.trans[1],Ti.trans[2],Ti.α,mlist[i],t)
    end
end

include("systems/rigidbodies.jl")
include("systems/navier_stokes.jl")
include("systems/utils.jl")

end
