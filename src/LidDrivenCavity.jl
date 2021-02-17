import RigidBodyTools: assign_velocity!

abstract type SpecifiedMotion end

struct LidDrivenCavity <: SpecifiedMotion
    u :: Float64
end

function assign_velocity!(u::AbstractVector{Float64},
                          v::AbstractVector{Float64},
                          b::Body,
                          m::LidDrivenCavity,
                          t::Real)

fill!(v,0.0)
fill!(u,0.0)

#get normal coordinates of the body surface
nx,ny = normalmid(b)

#find the top wall
top_wall = findall(x->(x=1.0),ny)

#set velocity to the top boundary
u[top_wall] .= m.u

return u,v

end
