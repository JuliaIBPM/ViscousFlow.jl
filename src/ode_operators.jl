#=
Defining the extra cache and extending prob_cache
=#

@ilmproblem ViscousIncompressibleFlow vector

struct ViscousIncompressibleFlowCache{CDT,FRT,DVT,VFT,VORT,DILT,VELT,FCT} <: AbstractExtraILMCache
   cdcache :: CDT
   fcache :: FRT
   dvb :: DVT
   vb_tmp :: DVT
   v_tmp :: VFT
   dv :: VFT
   v_rot ::VFT
   dv_tmp :: VFT
   w_tmp :: VORT
   divv_tmp :: DILT
   velcache :: VELT
   f :: FCT
end

function ImmersedLayers.prob_cache(prob::ViscousIncompressibleFlowProblem,
                                   base_cache::BasicILMCache{N,scaling}) where {N,scaling}
    @unpack phys_params, forcing, motions = prob
    @unpack gdata_cache, gcurl_cache, g = base_cache
    @unpack reference_body = motions

    dvb = zeros_surface(base_cache)
    vb_tmp = zeros_surface(base_cache)

    v_tmp = zeros_grid(base_cache)
    dv = zeros_grid(base_cache)
    v_rot = zeros_grid(base_cache)
    dv_tmp = zeros_grid(base_cache)
    w_tmp = zeros_gridcurl(base_cache)
    divv_tmp = zeros_griddiv(base_cache)

    velcache = VectorFieldCache(base_cache)

    # Construct a Lapacian outfitted with the viscosity
    Re = get_Reynolds_number(phys_params)
    over_Re = isinf(Re) ? 0.0 : 1.0/Re
    viscous_L = Laplacian(base_cache,over_Re)

    # Create cache for the convective derivative
    cdcache = reference_body > 0 ? RotConvectiveDerivativeCache(base_cache) : ConvectiveDerivativeCache(base_cache)

    fcache = nothing

    # Create cache for the forcing regions
    fmods = get_forcing_models(forcing)
    fcache = ForcingModelAndRegion(fmods,base_cache)

    # The state here is vorticity, the constraint is the surface traction
    f = _get_ode_function_list(viscous_L,base_cache)

    ViscousIncompressibleFlowCache(cdcache,fcache,dvb,vb_tmp,v_tmp,dv,v_rot,dv_tmp,w_tmp,divv_tmp,velcache,f)
end

_get_ode_function_list(viscous_L,base_cache::BasicILMCache{N}) where {N} =
                ODEFunctionList(state = zeros_gridcurl(base_cache),
                                constraint = zeros_surface(base_cache),
                                ode_rhs=viscousflow_vorticity_ode_rhs!,
                                lin_op=viscous_L,
                                bc_rhs=viscousflow_vorticity_bc_rhs!,
                                constraint_force = viscousflow_vorticity_constraint_force!,
                                bc_op = viscousflow_vorticity_bc_op!,
                                ode_implicit_rhs = viscousflow_vorticity_ode_implicit_rhs!)

                                #bc_regulator = viscousflow_vorticity_bc_reg!

_get_ode_function_list(viscous_L,base_cache::BasicILMCache{0}) =
                 ODEFunctionList(state = zeros_gridcurl(base_cache),
                                 ode_rhs=viscousflow_vorticity_ode_rhs!,
                                 lin_op=viscous_L)



#= ODE functions =#


function viscousflow_vorticity_ode_rhs!(dw,w,x,sys::ILMSystem,t)
  @unpack extra_cache, base_cache = sys
  @unpack v_tmp, dv = extra_cache

  velocity!(v_tmp,w,x,sys,t)
  viscousflow_velocity_ode_rhs!(dv,v_tmp,x,sys,t)
  curl!(dw,dv,base_cache)

  return dw
end

function viscousflow_velocity_ode_rhs!(dv,v,x,sys::ILMSystem,t)
    @unpack bc, phys_params, extra_cache, base_cache, motions = sys
    @unpack dvb, dv_tmp, cdcache, fcache, w_tmp = extra_cache

    Re = get_Reynolds_number(phys_params)
    over_Re = 1.0/Re

    fill!(dv,0.0)

    # Calculate the convective derivative
    convective_term!(dv_tmp,v,x,t,base_cache,extra_cache,phys_params,motions,cdcache)
    dv .-= dv_tmp

    # Calculate the double-layer term
    #prescribed_surface_jump!(dvb,x,t,sys)
    #surface_divergence_symm!(dv_tmp,over_Re*dvb,base_cache)
    #dv .-= dv_tmp

    # Apply forcing
    apply_forcing!(dv_tmp,v,t,fcache,phys_params)
    dv .+= dv_tmp

    return dv
end

 convective_term!(dv,v,x,t,base_cache,extra_cache,phys_params,motions,cdcache::ConvectiveDerivativeCache) =
            convective_derivative!(dv,v,base_cache,cdcache)

function convective_term!(dv,v,x,t,base_cache,extra_cache,phys_params,motions,cdcache::RotConvectiveDerivativeCache)
    @unpack v_rot, w_tmp = extra_cache

    v_rot .= v

    # to avoid having to deliver w as input here
    curl!(w_tmp,v_rot,base_cache) # v_rot is v' (relative to inertial frame)
    velocity_rel_to_rotating_frame!(v_rot,x,t,base_cache,phys_params,motions)  # Now v_rot is v̂ (rel. to rotating frame)
    w_cross_v!(dv,w_tmp,v_rot,base_cache,cdcache)
end

function viscousflow_vorticity_ode_implicit_rhs!(dw,x,sys::ILMSystem,t)
  @unpack extra_cache, base_cache = sys
  @unpack v_tmp, dv = extra_cache

  viscousflow_velocity_ode_implicit_rhs!(dv,x,sys,t)
  curl!(dw,dv,base_cache)

  return dw
end

function viscousflow_velocity_ode_implicit_rhs!(dv,x,sys::ILMSystem,t)
    @unpack bc, phys_params, extra_cache, base_cache, motions = sys
    @unpack dvb, dv_tmp, cdcache, fcache, w_tmp = extra_cache

    Re = get_Reynolds_number(phys_params)
    over_Re = 1.0/Re

    fill!(dv,0.0)

    # Calculate the double-layer term
    prescribed_surface_jump!(dvb,x,t,sys)
    surface_divergence_symm!(dv_tmp,over_Re*dvb,base_cache)
    dv .-= dv_tmp

    return dv
end


function viscousflow_vorticity_bc_rhs!(vb,x,sys::ILMSystem,t)
    @unpack bc, extra_cache, base_cache, phys_params, motions = sys
    @unpack dvb, vb_tmp, v_tmp, velcache, divv_tmp = extra_cache
    @unpack dcache, ϕtemp = velcache

    viscousflow_velocity_bc_rhs!(vb,x,sys,t)

    # Subtract influence of scalar potential field
    fill!(divv_tmp,0.0)
    prescribed_surface_jump!(dvb,x,t,sys)
    scalarpotential_from_masked_divv!(ϕtemp,divv_tmp,dvb,base_cache,dcache)
    vecfield_from_scalarpotential!(v_tmp,ϕtemp,base_cache)
    interpolate!(vb_tmp,v_tmp,base_cache)
    vb .-= vb_tmp

    # Subtract influence of free stream
    Uinf, Vinf = evaluate_freestream(t,x,sys)

    vb.u .-= Uinf
    vb.v .-= Vinf

    return vb
end



function viscousflow_velocity_bc_rhs!(vb,x,sys::ILMSystem,t)
    prescribed_surface_average!(vb,x,t,sys)
    return vb
end

function viscousflow_vorticity_constraint_force!(dw,τ,x,sys)
    @unpack extra_cache, base_cache = sys
    @unpack dv = extra_cache

    fill!(dv,0.0)
    viscousflow_velocity_constraint_force!(dv,τ,x,sys)
    curl!(dw,dv,base_cache)

    return dw
end

function viscousflow_velocity_constraint_force!(dv,τ,x,sys::ILMSystem{S,P,N}) where {S,P,N}
    @unpack base_cache = sys
    regularize!(dv,τ,base_cache)
    return dv
end

function viscousflow_velocity_constraint_force!(dv,τ,x,sys::ILMSystem{S,P,0}) where {S,P}
    return dv
end

function viscousflow_vorticity_bc_op!(vb,w,x,sys::ILMSystem)
    @unpack extra_cache, base_cache = sys
    @unpack velcache, v_tmp = extra_cache
    @unpack ψtemp = velcache

    vectorpotential_from_curlv!(ψtemp,w,base_cache)
    vecfield_from_vectorpotential!(v_tmp,ψtemp,base_cache)
    viscousflow_velocity_bc_op!(vb,v_tmp,x,sys)

    return vb
end

function viscousflow_velocity_bc_op!(vb,v,x,sys::ILMSystem)
    @unpack base_cache = sys
    interpolate!(vb,v,base_cache)
    return vb
end


function viscousflow_vorticity_bc_reg!(vb,τ,x,sys::ILMSystem{S,P,N}) where {S,P,N}
    @unpack phys_params = sys
    Re = get_Reynolds_number(phys_params)
    fill!(vb,0.0)
    #vb .= -cellsize(sys)*τ
    return vb

end

function viscousflow_vorticity_bc_reg!(vb,τ,x,sys::ILMSystem{S,P,0}) where {S,P}
    fill!(vb,0.0)
    return vb

end
