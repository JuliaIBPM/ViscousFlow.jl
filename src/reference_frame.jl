# Evaluating the free stream velocity


"""
    evaluate_freestream(t,x,sys; [subtract_translation = true])

Provide the components of the free stream velocity at time `t`.
If the problem is set up in the rotational frame of reference,
then this transforms the specified velocity and also subtracts the
translational velocity of the center of rotation. (If `subtract_translation = false`,
then it does not subtract this translation velocity.)
"""
function evaluate_freestream(t, x, sys; subtract_translation = true)
  @unpack phys_params, motions = sys
  evaluate_freestream(t,x,phys_params,motions; subtract_translation = subtract_translation)
end

evaluate_freestream(integ::ImmersedLayers.ConstrainedSystems.OrdinaryDiffEq.ODEIntegrator;kwargs...) =
    evaluate_freestream(integ.t,aux_state(integ.u),integ.p;kwargs...)

evaluate_freestream(sol::ImmersedLayers.ConstrainedSystems.OrdinaryDiffEq.ODESolution,sys::ILMSystem,t;kwargs...) =
    evaluate_freestream(t,aux_state(sol(t)),sys;kwargs...)

function evaluate_freestream(t, x, phys_params, motions; subtract_translation = true)
  @unpack reference_body, m, vl, Xl = motions

  # get the specified freestream, if any
  freestream_func = get_freestream_func(phys_params)
  Vinf = [freestream_func(t,phys_params)...]

  # if the problem is set in the rotational frame, then...
  if reference_body != 0

    # update the motions cache in sys
    evaluate_motion!(motions,x,t)

    # transform the specified freestream to the co-rotating coordinates
    XRref = rotation_transform(motions.Xl[reference_body])
    Vinf_plucker = XRref*PluckerMotion{2}(linear=Vinf)
    Vinf .= Vinf_plucker.linear

    vref = motions.vl[reference_body]

    if subtract_translation
        Vinf .-= vref.linear
    end
  end
  return tuple(Vinf...)
end




#= Changes of reference frame =#

"""
    velocity_rel_to_rotating_frame!(u_prime,x,t,base_cache,phys_params,motions)

Changes the input velocity `u_prime` (u', which is measured relative to the translating
frame) to û (measured relative to the translating/rotating frame)
"""
function velocity_rel_to_rotating_frame!(u_prime::Edges{Primal},x,t,base_cache,phys_params,motions)
    @unpack xg, yg = base_cache
    _velocity_rel_to_rotating_frame!(u_prime.u,u_prime.v,xg.v,yg.u,x,t,base_cache,motions)
    return u_prime
end

function velocity_rel_to_rotating_frame!(u_prime::VectorData,x,t,base_cache,phys_params,motions)
    @unpack pts = base_cache
    _velocity_rel_to_rotating_frame!(u_prime.u,u_prime.v,pts.u,pts.v,x,t,base_cache,motions)
    return u_prime
end

# Calculate cross product of Ω × (x-xr), where xr is the center of the
# reference body,
# and subtract this from u and v (velocity relative to inertial frame) to get
# û and v̂ (velocity relative to the rotating frame), returned in place
function _velocity_rel_to_rotating_frame!(u,v,x,y,xvec,t,base_cache,motions)
    @unpack reference_body, m, vl, Xl = motions
    @unpack bl = base_cache

    evaluate_motion!(motions,xvec,t)

    bref = bl[reference_body]
    xr, yr = bref.cent

    vref = vl[reference_body]

    Ω = vref.angular
    u .+= Ω.*(y .- yr)
    v .-= Ω.*(x .- xr)

    return u, v

end

"""
    velocity_rel_to_inertial_frame!(u_prime,x,t,base_cache,phys_params,motions)

Changes the input velocity `u_prime` (u', which is measured relative to the translating
frame) to u (measured relative to the inertial frame, but expressed in rotating coordinates).
Note that the inertial frame has zero velocity at infinity, so any free stream is
also removed here.
"""
function velocity_rel_to_inertial_frame!(u_prime::Edges{Primal},x,t,base_cache,phys_params,motions)
    Vinf = evaluate_freestream(t,x,phys_params,motions)
    Uxinf, Uyinf = Vinf
    u_prime.u .-= Uxinf
    u_prime.v .-= Uyinf
    return u_prime
end

function velocity_rel_to_inertial_frame!(u_prime::VectorData,x,t,base_cache,phys_params,motions)
    Vinf = evaluate_freestream(t,x,phys_params,motions)
    Uxinf, Uyinf = Vinf
    u_prime.u .-= Uxinf
    u_prime.v .-= Uyinf
    return u_prime
end






function _reference_body_acceleration(x,t,motions::ImmersedLayers.ILMMotion)
    @unpack reference_body, m = motions
    al = _body_accelerations(x,t,m)
    return al[reference_body]
end

# This is very kludgy and does not respect the rate of change of x
# There needs to be a proper calculation of accleration from RigidBodyTools
function _body_accelerations(x,t,m::RigidBodyMotion)
    dt = 1.0e-8
    vl_dt = body_velocities(x,t+dt,m)
    vl = body_velocities(x,t,m)

    return PluckerMotionList([PluckerMotion{2}((vl_dt[jb].data - vl[jb].data)/dt) for jb in 1:m.nbody])

end
