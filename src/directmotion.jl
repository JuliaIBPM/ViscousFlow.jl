
import RigidBodyTools: assign_velocity!

struct LidDrivenCavity <: DirectlySpecifiedMotion
    U :: Float64
end


"""
    assign_velocity!(u,v,x,y,m::DirectlySpecifiedMotion,t)

Given a `DirectlySpecifiedMotion`, vectors of `x` and `y` coordinates,
and time `t`, return the velocities `u` and `v` at those points.
"""
function assign_velocity!(u::AbstractVector{Float64},v::AbstractVector{Float64},
                 b::Body,
                 m::LidDrivenCavity,t::Real)

   fill!(v,0.0)
   fill!(u,0.0)
   # find the elements of the x and y vectors that lie on the
   # top wall, and assign velocity m.U to their u values.

   nx, ny = normalmid(b)
   topindices = findall(x -> (x ≈ 1.0), ny)
   u[topindices] .= m.U

   return u, v

end
