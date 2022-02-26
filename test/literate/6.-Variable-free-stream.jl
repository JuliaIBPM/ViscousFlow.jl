#=
# 6. Variable free stream
In this notebook we will simulate the flow with a time-varying free stream past a
stationary body. To demonstrate this, we will solve for oscillatory flow past a
rectangular object, in which the $x$ component of the free stream is

$$U_\infty(t) = A \sin(\Omega t)$$
=#
using ViscousFlow
#-
using Plots

#=
### Problem specification
We will set the Reynolds number to 200
=#
my_params = Dict()
my_params["Re"] = 200

#=
In order to set a time-varying free stream, we have to define a function
that provides the instantaneous free stream components and pass that
function into the system definition. In this function, we will
use the `Sinusoid` function (available via the `RigidBodyTools` module)
to create the modulated free stream. To demonstrate its possibilities,
we will pass in the parameters for the sinusoid via the `my_params` dictionary.
The "freestream average" specifies a mean free stream, if desired.
=#
my_params["freestream average"] = 0.0
my_params["freestream frequency"]  = 2.0
my_params["freestream amplitude"] = 1.0
my_params["freestream phase"] = 0.0

#=
Now we define the function. We can call it anything we want,
but it has to have the argument signature as shown. The
`Sinusoid` function is used, with the shift operator `>>`
to apply any phase lag.
=#
function my_freestream(t,phys_params)
    U = phys_params["freestream average"]
    Ω = phys_params["freestream frequency"]
    Ax = phys_params["freestream amplitude"]
    ϕx = phys_params["freestream phase"]
    Vinfmag = Ax*(RigidBodyTools.Sinusoid(Ω) >> (ϕx/Ω))
    Vinf_angle = get(phys_params,"freestream angle",0.0)

    Uinf = (U + Vinfmag(t))*cos(Vinf_angle)
    Vinf = (U + Vinfmag(t))*sin(Vinf_angle)
    return Uinf, Vinf
end

#=
The freestream function is generically considered a forcing function,
so we pass it in via the "freestream" key in the forcing dictionary.
=#
forcing_dict = Dict("freestream" => my_freestream)


# Now let us carry on with the other usual steps:

# ### Discretize
xlim = (-2.0,2.0)
ylim = (-1.5,1.5)
my_params["grid Re"] = 4.0
g = setup_grid(xlim,ylim,my_params)

#=
### Set up bodies
Here, we will set up a rectangle in the center of the domain
=#
Δs = surface_point_spacing(g,my_params)
body = Ellipse(0.25,0.5,Δs)

T = RigidTransform((0.0,0.0),0.0)
T(body)
#-
plot(body,xlim=xlim,ylim=ylim)

#=
### Construct the system structure
This step is like the previous notebook, but now we also provide the body and the freestream:
=#
sys = viscousflow_system(g,body,phys_params=my_params,forcing=forcing_dict);
#=
### Initialize
Now, we initialize with zero vorticity
=#
u0 = init_sol(sys)
# and create the integrator
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)

#=
### Solve
Now we are ready to solve the problem. Let's advance the solution to $t = 2.5$.
=#
@time step!(integrator,2.5)

#=
### Examine
Let's look at the flow field at the end of this interval
=#
sol = integrator.sol
plt = plot(layout = (2,2), size = (800, 600), legend=:false)
tsnap = 1.0:0.5:2.5
for (i, t) in enumerate(tsnap)
    plot!(plt[i],vorticity(sol,sys,t),sys,layers=false,title="t = $(round(t,digits=2))",clim=(-10,10),levels=range(-10,10,length=30),color = :RdBu)
    plot!(plt[i],surfaces(sol,sys,t))
end
plt

#=
#### Compute the force history
Just as we did for the stationary body in a constant free stream
=#
fx, fy = force(sol,sys,1);
# Plot them
plot(
plot(sol.t,2*fx,xlim=(0,Inf),ylim=(-6,6),xlabel="Convective time",ylabel="\$C_D\$",legend=:false),
plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-6,6),xlabel="Convective time",ylabel="\$C_L\$",legend=:false),
    size=(800,350)
)

# The mean drag and lift coefficients are
meanCD = GridUtilities.mean(2*fx[3:end])
#-
meanCL = GridUtilities.mean(2*fy[3:end])
