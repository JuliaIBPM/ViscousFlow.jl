module ViscousFlow

#using DocStringExtensions
using Reexport
using UnPack
@reexport using ImmersedLayers
@reexport using GridUtilities


export ViscousIncompressibleFlowProblem
export setup_grid, viscousflow_system, setup_problem, surface_point_spacing

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
    Δx = gridRe/get_Reynolds_number(phys_params)
end


function get_Reynolds_number(phys_params)
  haskey(phys_params,"Re") || error("No Reynolds number set")
  return phys_params["Re"]
end

function default_timestep(sys)
    @unpack phys_params, motions, forcing = sys
    g = get_grid(sys)
    Fo = get(phys_params,"Fourier",DEFAULT_FOURIER_NUMBER)
    Co = get(phys_params,"CFL",DEFAULT_CFL_NUMBER)
    Re = get_Reynolds_number(phys_params)

    Uscale = 1.0
    if !isnothing(motions)
      Uscale, _ = maxlistvelocity(sys)
    end


    Δt = min(Fo*Re*cellsize(g)^2,Co*cellsize(g)/Uscale)
    return Δt
end

function default_freestream(t,phys_params)
    Vinfmag = get(phys_params,"freestream speed",0.0)
    Vinf_angle = get(phys_params,"freestream angle",0.0)
    Uinf = Vinfmag*cos(Vinf_angle)
    Vinf = Vinfmag*sin(Vinf_angle)
    return Uinf, Vinf
end

function default_rotation(t,phys_params)
    omega = get(phys_params,"angular velocity",0.0)
    return omega
end

function default_vsplus(t,base_cache,phys_params,motions)
  vsplus=get_surface_velocity(phys_params,t,base_cache)
  #vsplus = zeros_surface(base_cache)
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
const DEFAULT_ROTATION_FUNC = default_rotation
const DEFAULT_TIMESTEP_FUNC = default_timestep
const DEFAULT_VSPLUS_FUNC = default_vsplus
const DEFAULT_VSMINUS_FUNC = default_vsminus


#=
Process keywords
=#

function get_freestream_func(forcing::Dict)
    return get(forcing,"freestream",DEFAULT_FREESTREAM_FUNC)
end

get_freestream_func(::Nothing) = get_freestream_func(Dict())

function get_rotation_func(phys_params)
    return get(phys_params,"rotation",DEFAULT_ROTATION_FUNC)
end

get_rotation_func(::Nothing) = get_rotation_func(Dict())

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
    viscous_L = Laplacian(base_cache,gcurl_cache,1.0/Re)

    # Create cache for the convective derivative
    cdcache = ConvectiveDerivativeCache(base_cache)

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
    setup_grid(xlim::Tuple,ylim::Tuple,phys_params::Dict)

Construct a Cartesian grid with limits `xlim` and `ylim`
and spacing determined by the Reynolds number in the `phys_params`.
"""
function setup_grid(xlim::Tuple,ylim::Tuple,phys_params)
    PhysicalGrid(xlim,ylim,grid_spacing(phys_params))
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
- `phys_params = ` A dictionary to pass in physical parameters
- `bc = ` A dictionary to pass in boundary condition data or functions, using "external"
          and "internal" keys to pass in functions that provide the
          corresponding surface data outside and inside the surface(s).
- `forcing = ` A dictionary to pass in forcing models (via the "forcing models" key),
                or to pass in an alternative model for
               the freestream velocity (with the "freestream" key)
- `motions = ` to provide function(s) that specify the velocity of the immersed surface(s). Note: if this keyword is used, it is assumed that surfaces will move.
- `timestep_func =` to pass in a function for time-dependent problems that provides the time-step size.
                  It is expected that this function takes in two arguments,
                  the `grid::PhysicalGrid` and `phys_params`, and returns the time step.
                  It defaults to the basic Fourier/CFL type function `default_timestep`
"""
function viscousflow_system(args...;kwargs...)
  prob = setup_problem(args...;kwargs...)
  return construct_system(prob)
end



#= ODE functions =#

function convective_derivative_rot!(uw::Edges{Primal},u::Edges{Primal},w::Nodes{Dual},base_cache::BasicILMCache,extra_cache::ConvectiveDerivativeCache)
    fill!(uw,0.0)
    if !iszero(w) && !iszero(u)
      _unscaled_convective_derivative_rot!(uw,u,w,extra_cache)
      #ImmersedLayers._scale_derivative!(uw,base_cache)
    end
end

function _unscaled_convective_derivative_rot!(uw::Edges{Primal},u::Edges{Primal},w::Nodes{Dual},extra_cache::ConvectiveDerivativeCache)
    @unpack vt1_cache, vt2_cache, vt3_cache = extra_cache

    #ensure these caches are at primal edges
    u_dual = Nodes(Dual,u)
    ucrossw = Edges(Primal,u) #replace this with vt2_cache

    grid_interpolate!(ucrossw.u,grid_interpolate!(u_dual, u.v) ∘ w)
    grid_interpolate!(ucrossw.v,grid_interpolate!(u_dual,-u.u) ∘ w)
    uw .= -ucrossw
    return uw

end

function get_rotational_vel!(u_prime::Edges{Primal},t,sys)
    @unpack bc, forcing, phys_params, extra_cache, base_cache = sys

    rot_func=get_rotation_func(phys_params)
    omega=rot_func(t,phys_params)

    xc,yc=get(phys_params,"center of rotation",(0.0,0.0))
    #ensure center is a tuple containing (xc,yc)

    xg, yg = x_grid(base_cache), y_grid(base_cache)

    #calculate cross product of omega and (x-xc).
    #rot_vel= zeros_grid(base_cache)
    u_prime.u .-= omega.*(yg.u .- yc)
    u_prime.v .+= omega.*(xg.v .- xc)

    return u_prime
    #contains both components u and v stored at primal edges.

end

function get_rotational_vel_primalnode(t,sys)
    @unpack bc, forcing, phys_params, extra_cache, base_cache = sys

    rot_func=get_rotation_func(phys_params)
    omega=rot_func(t,phys_params)

    xc,yc=get(phys_params,"center of rotation",(0.0,0.0))
    #ensure center is a tuple containing (xc,yc)

    xg, yg = x_grid(base_cache), y_grid(base_cache)

    #Should i define surfacescalarcache and get xg yg of cell center and do w cross x and return that or should i do surfacevectorcache
    #and then interpolate?

    #grid_interpolate xg and yg to location of primal node? doesnt mag do the required interpolation for us?

    #compute the vector (x-xc)
    x=xg-xc;
    y=yg-yc;

    #calculate cross product of omega and (x-xc).
    rot_vel= zeros_grid(base_cache)
    rot_vel.u .= -omega*y.u;
    rot_vel.v .=  omega*x.v;

    return magsq(rot_vel)
end

function get_surface_velocity(phys_params,t,base_cache)

    rot_func=get_rotation_func(phys_params)
    omega=rot_func(t,phys_params)

    xc,yc=get(phys_params,"center of rotation",(0.0,0.0))

    pts=points(base_cache)
    vs_tmp = zeros_surface(base_cache)
    vs_tmp.u.= -omega.*(pts.v-yc)
    vs_tmp.v.= omega.*(pts.u-xc)
    vs_tmp
end

function viscousflow_vorticity_ode_rhs!(dw,w,sys::ILMSystem,t)
  @unpack extra_cache, base_cache = sys
  @unpack v_tmp, dv = extra_cache

  velocity!(v_tmp,w,sys,t)
  viscousflow_velocity_ode_rhs!(dv,v_tmp,w,sys,t)
  curl!(dw,dv,base_cache)

  return dw
end

function viscousflow_velocity_ode_rhs!(dv,v,w,sys::ILMSystem,t)
    @unpack bc, forcing, phys_params, extra_cache, base_cache = sys
    @unpack dvb, v_rot, dv_tmp, cdcache, fcache = extra_cache

    Re = get_Reynolds_number(phys_params)
    over_Re = 1.0/Re

    fill!(dv,0.0)

    # Calculate the convective derivative
    v_rot .= v
    get_rotational_vel!(v_rot,t,sys)

    convective_derivative_rot!(dv_tmp,v_rot,w,base_cache,cdcache)
    #convective_derivative!(dv_tmp,v_rot,base_cache,cdcache)
    #dv_tmp is at primal edges
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


function viscousflow_vorticity_bc_rhs!(vb,sys::ILMSystem,t)
    @unpack bc, forcing,extra_cache, base_cache, phys_params = sys
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
    freestream_func = get_freestream_func(forcing)
    Uinf, Vinf = freestream_func(t,phys_params)
    vb.u .-= Uinf
    vb.v .-= Vinf
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
    @unpack forcing, phys_params, extra_cache, base_cache = sys
    @unpack dvb, velcache, divv_tmp, w_tmp = extra_cache

    prescribed_surface_jump!(dvb,t,sys)

    freestream_func = get_freestream_func(forcing)
    Vinf = freestream_func(t,phys_params)

    fill!(divv_tmp,0.0)
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
    @unpack phys_params, forcing, extra_cache, base_cache = sys
    @unpack velcache = extra_cache
    @unpack wcache = velcache

    freestream_func = get_freestream_func(forcing)
    Vinf = freestream_func(t,phys_params)

    streamfunction!(ψ,w,Vinf,base_cache,wcache)

end

streamfunction(w::Nodes{Dual},τ,sys::ILMSystem,t) = streamfunction!(zeros_gridcurl(sys),w,sys,t)

@snapshotoutput streamfunction

function pressure!(press::Nodes{Primal},w::Nodes{Dual},τ,sys::ILMSystem,t)
      @unpack extra_cache, base_cache = sys
      @unpack velcache, v_tmp, dv_tmp, dv, divv_tmp = extra_cache

      velocity!(v_tmp,w,sys,t)
      viscousflow_velocity_ode_rhs!(dv,v_tmp,w,sys,t)
      viscousflow_velocity_constraint_force!(dv_tmp,τ,sys)
      dv .-= dv_tmp


      #change to v_rot here?
      divergence!(divv_tmp,dv,base_cache)
      inverse_laplacian!(press,divv_tmp,base_cache)
      #whats the velocity that goes into magsq?
      press.-=(0.5*magsq(get_rotational_vel!(dv,t,sys)) - (0.5*get_rotational_vel_primalnode(t,sys)))
      return press
end

pressure(w::Nodes{Dual},τ,sys::ILMSystem,t) = pressure!(zeros_griddiv(sys),w,τ,sys,t)

@snapshotoutput pressure


#= Surface fields =#

function traction!(tract::VectorData{N},τ::VectorData{N},sys::ILMSystem,t) where {N}
    @unpack bc, extra_cache, base_cache, phys_params = sys
    @unpack vb_tmp, dvb = extra_cache
    @unpack sscalar_cache = base_cache

    prescribed_surface_average!(vb_tmp,t,sys)
    surface_velocity!(dvb,sys,t)
    vb_tmp .-= dvb

    nrm = normals(sys)
    pointwise_dot!(sscalar_cache,nrm,vb_tmp)
    prescribed_surface_jump!(dvb,t,sys)
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

force(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int) where {S,P} = nothing, nothing #Vector{Float64}(), Vector{Float64}()

function force(w::Nodes{Dual},τ::VectorData{N},sys::ILMSystem{S,P,N},t,bodyi::Int) where {S,P,N}
    @unpack base_cache = sys
    @unpack sdata_cache = base_cache
    traction!(sdata_cache,τ,sys,t)
    fx = integrate(sdata_cache.u,sys,bodyi)
    fy = integrate(sdata_cache.v,sys,bodyi)

    return fx, fy
end

@vectorsurfacemetric force

moment(w::Nodes{Dual},τ::VectorData{0},sys::ILMSystem{S,P,0},t,bodyi::Int;kwargs...) where {S,P} = nothing # Vector{Float64}()

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


end
