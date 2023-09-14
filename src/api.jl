#= API =#

setup_problem(g;kwargs...) =
    ViscousIncompressibleFlowProblem(g;timestep_func=DEFAULT_TIMESTEP_FUNC,
                                       kwargs...)

setup_problem(g,bl;bc=nothing,kwargs...) =
    ViscousIncompressibleFlowProblem(g,bl;timestep_func=DEFAULT_TIMESTEP_FUNC,
                                       bc=get_bc_func(bc),
                                       kwargs...)




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
- `bc = ` A dictionary to pass in boundary condition data or functions, using "external"
          and "internal" keys to pass in functions that provide the
          corresponding surface data outside and inside the surface(s).
- `forcing = ` A dictionary to pass in forcing models (via the "forcing models" key)
- `motions = ` to provide function(s) that specify the velocity of the immersed surface(s).
- `reference_body =` to solve the problem in a reference frame attached to a body rather than
              the inertial frame (with is the default, with body number 0).
- `timestep_func =` to pass in a function for time-dependent problems that provides the time-step size.
                  It is expected that this function takes in two arguments,
                  the `grid::PhysicalGrid` and `phys_params`, and returns the time step.
                  It defaults to the basic Fourier/CFL type function `default_timestep`
"""
function viscousflow_system(args...;phys_params=Dict(), kwargs...)
  prob = setup_problem(args...;phys_params=phys_params,kwargs...)
  return construct_system(prob)
end
