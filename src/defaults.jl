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

function default_timestep(u,sys)
    @unpack phys_params = sys
    g = get_grid(sys)
    Fo = get(phys_params,"Fourier",DEFAULT_FOURIER_NUMBER)
    Co = get(phys_params,"CFL",DEFAULT_CFL_NUMBER)
    Re = get_Reynolds_number(phys_params)

    Umax = get_max_velocity(u,sys)

    Δt = min(Fo*Re*cellsize(g)^2,Co*cellsize(g)/Umax)
    return Δt
end

function get_max_velocity(u,sys)
    @unpack base_cache, phys_params, motions = sys
    @unpack m = motions

    x = aux_state(u)

    Umax = 0.0
    if !isnothing(m)
      Umax,i,tmax,bi = maxvelocity(u,sys)
    end

    for t in 0:0.05:10
      Uinf, Vinf = evaluate_freestream(t,x,sys)
      Umax = max(Umax,sqrt(Uinf^2+Vinf^2))
    end

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


function default_vsplus(t,x,base_cache,phys_params,motions)
  vsplus = zeros_surface(base_cache)
  return vsplus
end

function default_vsminus(t,x,base_cache,phys_params,motions)
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
