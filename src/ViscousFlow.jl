module ViscousFlow

#using DocStringExtensions
using Reexport
using UnPack
@reexport using ImmersedLayers
@reexport using GridUtilities


export ViscousIncompressibleFlowProblem
export setup_grid, viscousflow_system, setup_problem, surface_point_spacing,
        surface_velocity_in_translating_frame!

#= Supporting functions =#

setup_problem(g;kwargs...) =
    ViscousIncompressibleFlowProblem(g;timestep_func=DEFAULT_TIMESTEP_FUNC,
                                       kwargs...)

setup_problem(g,bl;bc=nothing,kwargs...) =
    ViscousIncompressibleFlowProblem(g,bl;timestep_func=DEFAULT_TIMESTEP_FUNC,
                                       bc=get_bc_func(bc),
                                       kwargs...)


function grid_spacing(phys_params)
    gridRe = get(phys_params,"grid Re",DEFAULT_GRID_RE)
    Re = get_Reynolds_number(phys_params)
    Δx = gridRe/Re
end

grid_spacing(Δx::Float64) = Δx

function get_Reynolds_number(phys_params)
  #haskey(phys_params,"Re") || error("No Reynolds number set")
  #return phys_params["Re"]
  return get(phys_params,"Re",Inf)
end

function default_timestep(sys)
    @unpack phys_params = sys
    g = get_grid(sys)
    Fo = get(phys_params,"Fourier",DEFAULT_FOURIER_NUMBER)
    Co = get(phys_params,"CFL",DEFAULT_CFL_NUMBER)
    Re = get_Reynolds_number(phys_params)

    Umax = get_max_velocity(sys)

    Δt = min(Fo*Re*cellsize(g)^2,Co*cellsize(g)/Umax)
    return Δt
end

function get_max_velocity(sys)
    @unpack base_cache, phys_params, motions = sys

    Umax = 0.0
    if !isnothing(motions)
      Umax,i,tmax,bi = maxlistvelocity(sys)
    end

    Uinf, Vinf = evaluate_freestream(0.0,phys_params)
    Umax = max(Umax,sqrt(Uinf^2+Vinf^2))

    mot = get_rotation_func(phys_params)
    Umot,i,tmax,bi = maxlistvelocity(surfaces(sys),mot)
    Umax = max(Umax,Umot)

    # If no velocity has been set yet, just set it to unity
    Umax = Umax == 0.0 ? 1.0 : Umax

    return Umax
end


function default_freestream(t,phys_params)
    Vinfmag = get(phys_params,"freestream speed",0.0)
    Vinf_angle = get(phys_params,"freestream angle",0.0)

    Uinf = Vinfmag*cos(Vinf_angle)
    Vinf = Vinfmag*sin(Vinf_angle)

    return Uinf, Vinf
end


function default_vsplus(t,base_cache,phys_params,motions)
  vsplus = zeros_surface(base_cache)
  return vsplus
end

function default_vsminus(t,base_cache,phys_params,motions)
    vsminus = zeros_surface(base_cache)
    return vsminus
end



const DEFAULT_GRID_RE = 2.0
const DEFAULT_FOURIER_NUMBER = 1.0
const DEFAULT_CFL_NUMBER = 0.5
const DEFAULT_DS_TO_DX_RATIO = 1.4
const DEFAULT_FREESTREAM_FUNC = default_freestream
const DEFAULT_TIMESTEP_FUNC = default_timestep
const DEFAULT_VSPLUS_FUNC = default_vsplus
const DEFAULT_VSMINUS_FUNC = default_vsminus
const DEFAULT_CENTER_OF_ROTATION = (0.0,0.0)

#=
Process keywords
=#

function get_freestream_func(phys_params::Dict)
    return get(phys_params,"freestream",DEFAULT_FREESTREAM_FUNC)
end

get_freestream_func(::Nothing) = get_freestream_func(Dict())


function get_rotation_func(phys_params::Dict)
    omega = get(phys_params,"angular velocity",0.0)
    Xp = get_center_of_rotation(phys_params)
    return get(phys_params,"reference frame",
                  RigidBodyTools.RigidBodyMotion((0.0,0.0),omega;pivot=Xp))
end

get_rotation_func(::Nothing) = get_rotation_func(Dict())

in_rotational_frame(phys_params::Dict) = haskey(phys_params,"angular velocity") || haskey(phys_params,"reference frame")

function get_center_of_rotation(phys_params::Dict)
    return get(phys_params,"center of rotation",DEFAULT_CENTER_OF_ROTATION)
end

function get_forcing_models(forcing::Dict)
    return get(forcing,"forcing models",nothing)
end

get_forcing_models(::Nothing) = get_forcing_models(Dict())

function get_bc_func(bc_in::Dict)
    bc = Dict()
    bc["exterior"] = haskey(bc_in,"exterior") ? bc_in["exterior"] : DEFAULT_VSPLUS_FUNC
    bc["interior"] = haskey(bc_in,"interior") ? bc_in["interior"] : DEFAULT_VSMINUS_FUNC
    return bc
end

get_bc_func(::Nothing) = get_bc_func(Dict())


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
    @unpack phys_params, forcing = prob
    @unpack gdata_cache, gcurl_cache, g = base_cache

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
    cdcache = in_rotational_frame(phys_params) ? RotConvectiveDerivativeCache(base_cache) : ConvectiveDerivativeCache(base_cache)

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
                                bc_op = viscousflow_vorticity_bc_op!)

_get_ode_function_list(viscous_L,base_cache::BasicILMCache{0}) =
                 ODEFunctionList(state = zeros_gridcurl(base_cache),
                                 ode_rhs=viscousflow_vorticity_ode_rhs!,
                                 lin_op=viscous_L)

#= setup API =#

"""
    setup_grid(xlim::Tuple,ylim::Tuple,phys_params::Dict[;nthreads_max=length(Sys.cpu_info())])

Construct a Cartesian grid with limits `xlim` and `ylim`
and spacing determined by the Reynolds number in the `phys_params`.
The maximum number of threads can be optionally set; it defaults
to the number of processor cores.
"""
function setup_grid(xlim::Tuple,ylim::Tuple,phys_params;kwargs...)
    PhysicalGrid(xlim,ylim,grid_spacing(phys_params);kwargs...)
end

"""
    surface_point_spacing(g::PhysicalGrid,phys_params)

Calculate the surface point spacing for a given grid, using
the specified parameter "point spacing ratio" in the physical parameters,
or the default value (1.4) if not specified.
"""
function surface_point_spacing(g::PhysicalGrid,phys_params)
    ds_to_dx = get(phys_params,"point spacing ratio",DEFAULT_DS_TO_DX_RATIO)
    return ds_to_dx*cellsize(g)
end

"""
    viscousflow_system(grid,[bodies];kwargs...)

Construct the operators and cache variables for a viscous flow problem.

The `kwargs` are the optional keyword aguments. There are several, some of
which are crucial for certain types of problems.

- `ddftype = ` to set the DDF type. The default is `CartesianGrids.Yang3`.
- `scaling = ` to set the scaling type, `GridScaling` (default) or `IndexScaling`.
- `dtype = ` to set the element type to `Float64` (default) or `ComplexF64`.
- `phys_params = ` A dictionary to pass in physical parameters, or to pass in
                  alternative models for the freestream velocity (with the "freestream" key)
                  or overall rotational motion (with the "rotation" key)
- `bc = ` A dictionary to pass in boundary condition data or functions, using "external"
          and "internal" keys to pass in functions that provide the
          corresponding surface data outside and inside the surface(s).
- `forcing = ` A dictionary to pass in forcing models (via the "forcing models" key)
- `motions = ` to provide function(s) that specify the velocity of the immersed surface(s). Note: if this keyword is used, it is assumed that surfaces will move.
- `timestep_func =` to pass in a function for time-dependent problems that provides the time-step size.
                  It is expected that this function takes in two arguments,
                  the `grid::PhysicalGrid` and `phys_params`, and returns the time step.
                  It defaults to the basic Fourier/CFL type function `default_timestep`
"""
function viscousflow_system(args...;phys_params=Dict(), kwargs...)
  prob = setup_problem(args...;phys_params=phys_params,kwargs...)
  return construct_system(prob)
end

# Evaluating the free stream velocity

"""
    evaluate_freestream(t,phys_params)

Provide the components of the free stream velocity at time `t`.
If the problem is set up in the rotational frame of reference,
then this transforms the specified velocity and also subtracts the
translational motion of the center of rotation.
"""
function evaluate_freestream(t,phys_params)

  # get the specified freestream, if any
  freestream_func = get_freestream_func(phys_params)
  Vinf = freestream_func(t,phys_params)

  # if the problem is set in the rotational frame, then...
  if in_rotational_frame(phys_params)

    # transform the specified freestream to the co-rotating coordinates
    Vinf_i = Vinf
    Vinf = transform_vector_to_rotating_coordinates(Vinf_i...,t,phys_params)

    # subtract the motion of the center of rotation
    mot = get_rotation_func(phys_params)
    Xp = get_center_of_rotation(phys_params)
    Vp = translational_velocity(mot(t,Xp);inertial=false)
    Vinf = Vinf .- Vp
  end
  return Vinf
end

"""
    surface_velocity_in_translating_frame!(vel,t,base_cache,phys_params)

Evaluates the surface velocity when the problem is set up in a moving
frame of reference (specified with the `reference frame` keyword).
It removes the translational velocity of the center of rotation, since
this is applied as a free stream velocity (with change of sign).
"""
function surface_velocity_in_translating_frame!(vel,t,base_cache,phys_params)
    mot = get_rotation_func(phys_params)

    surface_velocity!(vel,base_cache,mot,t;inertial=false)
    Xp = get_center_of_rotation(phys_params)
    Ũp, Ṽp = translational_velocity(mot(t,Xp);inertial=false)
    vel.u .-= Ũp
    vel.v .-= Ṽp
    return vel
end


#= Changes of reference frame =#

"""
    velocity_rel_to_rotating_frame!(u_prime,t,sys::ILMSystem)

Changes the input velocity `u_prime` (u', which is measured relative to the translating
frame) to û (measured relative to the translating/rotating frame)
"""
function velocity_rel_to_rotating_frame!(u_prime::Edges{Primal},t,base_cache,phys_params)
    xg, yg = x_grid(base_cache), y_grid(base_cache)
    _velocity_rel_to_rotating_frame!(u_prime.u,u_prime.v,xg.v,yg.u,t,phys_params)
    return u_prime
end

function velocity_rel_to_rotating_frame!(u_prime::VectorData,t,base_cache,phys_params)
    pts = points(base_cache)
    _velocity_rel_to_rotating_frame!(u_prime.u,u_prime.v,pts.u,pts.v,t,phys_params)
    return u_prime
end

# Calculate cross product of Ω × (x-xr).
# and subtract this from u and v (velocity relative to inertial frame) to get
# û and v̂ (velocity relative to the rotating frame), returned in place
function _velocity_rel_to_rotating_frame!(u,v,x,y,t,phys_params)
    rot_func = get_rotation_func(phys_params)
    kindata = rot_func(t)
    omega = angular_velocity(kindata)

    xr, yr = get_center_of_rotation(phys_params)

    u .+= omega.*(y .- yr)
    v .-= omega.*(x .- xr)
    return u, v

end

"""
    transform_vector_to_rotating_coordinates(u,v,t,phys_params) -> Tuple

Transform a vector expressed in inertial coordinates with components `u`
and `v` to the rotating coordinate system at time `t`.
"""
function transform_vector_to_rotating_coordinates(u,v,t,phys_params)
  # This computes R^T*V to provide vector components in the corotating
  # coordinate system
  mot = get_rotation_func(phys_params)
  k = mot(t)
  α = angular_position(k)

  up =  u*cos(α) + v*sin(α)
  vp = -u*sin(α) + v*cos(α)
  return up, vp
end

"""
    transform_vector_to_rotating_coordinates(u,v,t,phys_params) -> Tuple

Transform a vector expressed in rotating coordinates with components `u`
and `v` to the inertial coordinate system at time `t`.
"""
function transform_vector_to_inertial_coordinates(u,v,t,phys_params)
  # This computes R^T*Vinf to provide freestream
  # velocity in the corotating coordinate system, if appropriate
  mot = get_rotation_func(phys_params)
  k = mot(t)
  α = angular_position(k)

  up =  u*cos(α) - v*sin(α)
  vp =  u*sin(α) + v*cos(α)
  return up, vp
end

#= ODE functions =#


function viscousflow_vorticity_ode_rhs!(dw,w,sys::ILMSystem,t)
  @unpack extra_cache, base_cache = sys
  @unpack v_tmp, dv = extra_cache

  velocity!(v_tmp,w,sys,t)
  viscousflow_velocity_ode_rhs!(dv,v_tmp,sys,t)
  curl!(dw,dv,base_cache)

  return dw
end

function viscousflow_velocity_ode_rhs!(dv,v,sys::ILMSystem,t)
    @unpack bc, phys_params, extra_cache, base_cache = sys
    @unpack dvb, dv_tmp, cdcache, fcache, w_tmp = extra_cache

    Re = get_Reynolds_number(phys_params)
    over_Re = 1.0/Re

    fill!(dv,0.0)

    # Calculate the convective derivative
    convective_term!(dv_tmp,v,t,base_cache,extra_cache,phys_params,cdcache)
    dv .-= dv_tmp

    # Calculate the double-layer term
    prescribed_surface_jump!(dvb,t,sys)
    surface_divergence_symm!(dv_tmp,over_Re*dvb,base_cache)
    dv .-= dv_tmp

    # Apply forcing
    apply_forcing!(dv_tmp,v,t,fcache,phys_params)
    dv .+= dv_tmp

    return dv
end

 convective_term!(dv,v,t,base_cache,extra_cache,phys_params,cdcache::ConvectiveDerivativeCache) =
            convective_derivative!(dv,v,base_cache,cdcache)

function convective_term!(dv,v,t,base_cache,extra_cache,phys_params,cdcache::RotConvectiveDerivativeCache)
    @unpack v_rot, w_tmp = extra_cache

    v_rot .= v

    # to avoid having to deliver w as input here
    curl!(w_tmp,v_rot,base_cache) # v_rot is v' (relative to inertial frame)
    velocity_rel_to_rotating_frame!(v_rot,t,base_cache,phys_params)  # Now v_rot is v̂ (rel. to rotating frame)
    w_cross_v!(dv,w_tmp,v_rot,base_cache,cdcache)
end

function viscousflow_vorticity_bc_rhs!(vb,sys::ILMSystem,t)
    @unpack bc, extra_cache, base_cache, phys_params = sys
    @unpack dvb, vb_tmp, v_tmp, velcache, divv_tmp = extra_cache
    @unpack dcache, ϕtemp = velcache

    viscousflow_velocity_bc_rhs!(vb,sys,t)

    # Subtract influence of scalar potential field
    fill!(divv_tmp,0.0)
    prescribed_surface_jump!(dvb,t,sys)
    scalarpotential_from_masked_divv!(ϕtemp,divv_tmp,dvb,base_cache,dcache)
    vecfield_from_scalarpotential!(v_tmp,ϕtemp,base_cache)
    interpolate!(vb_tmp,v_tmp,base_cache)
    vb .-= vb_tmp

    # Subtract influence of free stream
    Uinf, Vinf = evaluate_freestream(t,phys_params)

    vb.u .-= Uinf
    vb.v .-= Vinf

    return vb
end



function viscousflow_velocity_bc_rhs!(vb,sys::ILMSystem,t)
    prescribed_surface_average!(vb,t,sys)
    return vb
end

function viscousflow_vorticity_constraint_force!(dw,τ,sys)
    @unpack extra_cache, base_cache = sys
    @unpack dv = extra_cache

    viscousflow_velocity_constraint_force!(dv,τ,sys)
    curl!(dw,dv,base_cache)

    return dw
end

function viscousflow_velocity_constraint_force!(dv,τ,sys::ILMSystem{S,P,N}) where {S,P,N}
    @unpack base_cache = sys
    regularize!(dv,τ,base_cache)
    return dv
end

function viscousflow_velocity_constraint_force!(dv,τ,sys::ILMSystem{S,P,0}) where {S,P}
    return dv
end

function viscousflow_vorticity_bc_op!(vb,w,sys::ILMSystem)
    @unpack extra_cache, base_cache = sys
    @unpack velcache, v_tmp = extra_cache
    @unpack ψtemp = velcache

    vectorpotential_from_curlv!(ψtemp,w,base_cache)
    vecfield_from_vectorpotential!(v_tmp,ψtemp,base_cache)
    viscousflow_velocity_bc_op!(vb,v_tmp,sys)

    return vb
end

function viscousflow_velocity_bc_op!(vb,v,sys::ILMSystem)
    @unpack base_cache = sys
    interpolate!(vb,v,base_cache)
    return vb
end



#= Fields =#

vorticity!(masked_w::Nodes{Dual},w::Nodes{Dual},dv::VectorData,base_cache::BasicILMCache,wcache::VectorPotentialCache) =
            masked_curlv_from_curlv_masked!(masked_w,w,dv,base_cache,wcache)

total_vorticity!(w::Nodes{Dual},masked_w::Nodes{Dual},dv::VectorData,base_cache::BasicILMCache,wcache::VectorPotentialCache) =
            curlv_masked_from_masked_curlv!(w,masked_w,dv,base_cache,wcache)

for f in [:vorticity,:total_vorticity]

    f! = Symbol(string(f)*"!")

    @eval function $f!(masked_w::Nodes{Dual},w::Nodes{Dual},τ,sys::ILMSystem,t)
        @unpack extra_cache, base_cache = sys
        @unpack dvb, velcache = extra_cache
        @unpack wcache = velcache

        prescribed_surface_jump!(dvb,t,sys)
        $f!(masked_w,w,dvb,base_cache,wcache)
    end

    @eval $f(w::Nodes{Dual},τ,sys::ILMSystem,t) = $f!(zeros_gridcurl(sys),w,τ,sys,t)

    @eval @snapshotoutput $f
end

function velocity!(v::Edges{Primal},curl_vmasked::Nodes{Dual},masked_divv::Nodes{Primal},dv::VectorData,vp,base_cache::BasicILMCache,velcache::VectorFieldCache,w_tmp::Nodes{Dual})
    @unpack wcache  = velcache

    masked_curlv = w_tmp

    masked_curlv_from_curlv_masked!(masked_curlv,curl_vmasked,dv,base_cache,wcache)
    vecfield_helmholtz!(v,masked_curlv,masked_divv,dv,vp,base_cache,velcache)

    return v
end

function velocity!(v::Edges{Primal},w::Nodes{Dual},sys::ILMSystem,t)
    @unpack phys_params, extra_cache, base_cache = sys
    @unpack dvb, velcache, divv_tmp, w_tmp = extra_cache

    prescribed_surface_jump!(dvb,t,sys)

    Vinf = evaluate_freestream(t,phys_params)

    fill!(divv_tmp,0.0)
    fill!(w_tmp,0.0)
    velocity!(v,w,divv_tmp,dvb,Vinf,base_cache,velcache,w_tmp)
end

velocity(w::Nodes{Dual},τ,sys::ILMSystem,t) = velocity!(zeros_grid(sys),w,sys,t)

@snapshotoutput velocity


function streamfunction!(ψ::Nodes{Dual},w::Nodes{Dual},vp::Tuple,base_cache::BasicILMCache,wcache::VectorPotentialCache)
    @unpack stemp = wcache

    fill!(ψ,0.0)
    vectorpotential_from_curlv!(stemp,w,base_cache)
    ψ .+= stemp

    vectorpotential_uniformvecfield!(stemp,vp...,base_cache)
    ψ .+= stemp

    return ψ
end

function streamfunction!(ψ::Nodes{Dual},w::Nodes{Dual},sys::ILMSystem,t)
    @unpack phys_params, extra_cache, base_cache = sys
    @unpack velcache = extra_cache
    @unpack wcache = velcache

    Vinf = evaluate_freestream(t,phys_params)

    streamfunction!(ψ,w,Vinf,base_cache,wcache)

end

streamfunction(w::Nodes{Dual},τ,sys::ILMSystem,t) = streamfunction!(zeros_gridcurl(sys),w,sys,t)

@snapshotoutput streamfunction

function pressure!(press::Nodes{Primal},w::Nodes{Dual},τ,sys::ILMSystem,t)
      @unpack extra_cache, base_cache, phys_params = sys
      @unpack velcache, v_tmp, dv_tmp, dv, divv_tmp, v_rot = extra_cache

      velocity!(v_tmp,w,sys,t)
      viscousflow_velocity_ode_rhs!(dv,v_tmp,sys,t)
      viscousflow_velocity_constraint_force!(dv_tmp,τ,sys)
      dv .-= dv_tmp

      divergence!(divv_tmp,dv,base_cache)
      inverse_laplacian!(press,divv_tmp,base_cache)

      # For rotating coordinate systems...
      if in_rotational_frame(phys_params)
        # Compute v̂ here
        velocity_rel_to_rotating_frame!(v_tmp,t,base_cache,phys_params)

        # Compute Ω × x̂ here
        fill!(v_rot,0.0)
        velocity_rel_to_rotating_frame!(v_rot,t,base_cache,phys_params)

        press .-= 0.5*magsq(v_tmp) - 0.5*magsq(v_rot)
      end

      return press
end

pressure(w::Nodes{Dual},τ,sys::ILMSystem,t) = pressure!(zeros_griddiv(sys),w,τ,sys,t)

@snapshotoutput pressure

function convective_acceleration!(vdv::Edges{Primal},w::Nodes{Dual},sys::ILMSystem,t)
    @unpack extra_cache, base_cache = sys
    @unpack v_tmp, cdcache = extra_cache
    velocity!(v_tmp,w,sys,t)
    convective_derivative!(vdv,v_tmp,base_cache,cdcache)
    return vdv
end

convective_acceleration(w::Nodes{Dual},τ,sys::ILMSystem,t) = convective_acceleration!(zeros_grid(sys),w,sys,t)

@snapshotoutput convective_acceleration


#= Surface fields =#

function traction!(tract::VectorData{N},τ::VectorData{N},sys::ILMSystem,t) where {N}
    @unpack bc, extra_cache, base_cache, phys_params = sys
    @unpack vb_tmp, dvb = extra_cache
    @unpack sscalar_cache = base_cache

    # (v̅ - Ẋ)⋅n -> sscalar_cache
    prescribed_surface_average!(vb_tmp,t,sys)
    surface_velocity!(dvb,sys,t)
    vb_tmp .-= dvb
    if in_rotational_frame(phys_params)
        velocity_rel_to_rotating_frame!(vb_tmp,t,base_cache,phys_params)
    end
    nrm = normals(sys)
    pointwise_dot!(sscalar_cache,nrm,vb_tmp)

    # [v] -> dvb
    prescribed_surface_jump!(dvb,t,sys)

    # [v](v̅ - Ẋ)⋅n
    product!(tract,dvb,sscalar_cache)

    tract .+= τ

    return tract

end
traction!(out::VectorData{0},τ::VectorData{0},sys::ILMSystem,t) = out
traction(w::Nodes{Dual},τ::VectorData,sys::ILMSystem,t) = traction!(zeros_surface(sys),τ,sys,t)
@snapshotoutput traction

function pressurejump!(dpb::ScalarData{N},τ::VectorData{N},sys::ILMSystem,t) where {N}
    @unpack base_cache = sys
    @unpack sdata_cache = base_cache

    nrm = normals(sys)
    traction!(sdata_cache,τ,sys,t)
    pointwise_dot!(dpb,nrm,sdata_cache)
    dpb .*= -1.0
    return dpb
end
pressurejump!(out::VectorData{0},τ::VectorData{0},sys::ILMSystem,t) = out
pressurejump(w::Nodes{Dual},τ::VectorData,sys::ILMSystem,t) = pressurejump!(zeros_surfacescalar(sys),τ,sys,t)
@snapshotoutput pressurejump


#= Integrated metrics =#

force(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing, nothing #Vector{Float64}(), Vector{Float64}()

"""
    force(sol,sys,bodyi[;inertial=true]) -> Tuple{Vector}

Calculated the force exerted by the fluid on body `bodyi` from the computational solution `sol` of system `sys`.
It returns the force history as a tuple of arrays: one array for each component.
If `inertial=true` (default), then the components are provided in the inertial
coordinate system. Otherwise, they are in the body coordinate system.
""" force(sol,sys,bodyi)

function force(w::Nodes{Dual},τ::VectorData{N},sys::ILMSystem{S,P,N},t,bodyi::Int;inertial=true) where {S,P,N}
    @unpack base_cache, phys_params = sys
    @unpack sdata_cache = base_cache
    traction!(sdata_cache,τ,sys,t)
    fx = integrate(sdata_cache.u,sys,bodyi)
    fy = integrate(sdata_cache.v,sys,bodyi)

    if (inertial && in_rotational_frame(phys_params))
        fx_r, fy_r = fx, fy
        fx, fy = transform_vector_to_inertial_coordinates(fx_r,fy_r,t,phys_params)
    end

    return fx, fy
end

@vectorsurfacemetric force

moment(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

"""
    moment(sol,sys,bodyi[;center=(0,0)]) -> Vector

Calculated the momemt exerted by the fluid on body `bodyi` from the computational solution `sol` of system `sys`.
It returns the moment history as an array. The moment is calculated about center
`center`, which defaults to `(0,0)` in whichever coordinate system the problem is solved.
""" moment(sol,sys,bodyi)


function moment(w::Nodes{Dual},τ::VectorData{N},sys::ILMSystem{S,P,N},t,bodyi::Int;center=(0.0,0.0)) where {S,P,N}
    @unpack base_cache = sys
    @unpack sdata_cache, sscalar_cache = base_cache
    xc, yc = center
    pts = points(sys)
    pts.u .-= xc
    pts.v .-= yc

    traction!(sdata_cache,τ,sys,t)
    pointwise_cross!(sscalar_cache,pts,sdata_cache)
    mom = integrate(sscalar_cache,sys,bodyi)

    return mom
end

@scalarsurfacemetric moment


power(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

"""
    power(sol,sys,bodyi)

Calculated the history of the total rate of work done by the flow on body `bodyi` (or on the flow by the body, if negative)
from the computational solution `sol` of system `sys`.
""" power(sol,sys,bodyi)


function power(w::Nodes{Dual},τ::VectorData{N},sys::ILMSystem{S,P,N},t,bodyi::Int;inertial=true) where {S,P,N}
    @unpack phys_params = sys
    mot = get_rotation_func(phys_params)

    Ω = angular_velocity(mot(t))
    U, V = translational_velocity(mot(t))

    mom = moment(w,τ,sys,t,bodyi)
    fx, fy = force(w,τ,sys,t,bodyi;inertial=inertial)

    pow = Ω*mom + fx*U + fy*V

    return pow

end

@scalarsurfacemetric power


extracted_power(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

"""
    extracted_power(sol,sys,bodyi)

Calculated the power extracted from the flow (or power input into the flow, if negative)
by body `bodyi` from the computational solution `sol` of system `sys`. This computes the
power from the force and velocity due to prescribed oscillations. It does not
include the power due to steady body motion.
""" extracted_power(sol,sys,bodyi)


function extracted_power(w::Nodes{Dual},τ::VectorData{N},sys::ILMSystem{S,P,N},t,bodyi::Int;inertial=true) where {S,P,N}
    @unpack phys_params = sys
    mot = get_rotation_func(phys_params)

    Ω = angular_velocity(mot(t))
    U, V = translational_velocity(mot(t))

    # Remove mean part, since only the unsteady motion is included
    U -= _mean_x_velocity(mot.kin)
    V -= _mean_y_velocity(mot.kin)

    mom = moment(w,τ,sys,t,bodyi)
    fx, fy = force(w,τ,sys,t,bodyi;inertial=inertial)

    pow = Ω*mom + fx*U + fy*V

    return pow

end

@scalarsurfacemetric extracted_power


# These routines calculate the mean velocity components
_mean_x_velocity(kin) = 0.0
_mean_y_velocity(kin) = 0.0
_mean_x_velocity(kin::Oscillation) = kin.Ux
_mean_x_velocity(kin::Pitchup) = kin.U₀
_mean_y_velocity(kin::Oscillation) = kin.Uy
_mean_y_velocity(kin::Pitchup) = 0.0


end
