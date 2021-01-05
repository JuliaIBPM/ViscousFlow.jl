#### Operators for a system with moving body ####

export rigid_body_rhs!, rigid_body_rhs, zero_body_state

zero_body_state(b::Body) = zeros(Float64,3*(NDIM-1))
zero_body_state(bl::BodyList) = zeros(Float64,3*length(bl)*(NDIM-1))
zero_body_state(sys::NavierStokes) = zero_body_state(sys.bodies)


## Right-hand sides for rigid-body motion equations ##

function rigid_body_rhs!(dx::Vector{T},x::Vector{T},sys::NavierStokes,t::Real) where {T<:Real}
  # note that x itself is not used here, because sys is presumed to contain
  # updated body data
  length(dx) == 3*length(sys.bodies)*(NDIM-1) || error("Wrong length for vector")
  dx .= rigidbodyvelocity(sys.motions,t)
  return dx
end


rigid_body_rhs(x::Vector{T},sys::NavierStokes,t::Real) where {T<:Real} =
              rigid_body_rhs!(zero_body_state(sys),sys,t)
