#=
# 7. Lid Driven Cavity Flow
In this notebook we will simulate the flow with a top moving wall. To demonstrate this, we will solve for internal flow in a square cavity by enforcing no-slip at the walls.
=#
using ViscousFlow
#-
using Plots

#=
##Problem specification
Take $Re=100$ for example. We will set the Reynolds number to 100
=#
Re = 100 # Reynolds number

#=
##Discretization
Note that the rectangle function used for making the cavity shape requires a specified half length. The immersed boundary projection method for internal flow requires the size of the domain to be at least a step size greater at the boundaries (i.e. half length + Δx).
=#
Δt,Δx = setstepsizes(Re,gridRe=1.6)
halflength=0.5 # to make rectangle with side length of 1
domain_lim=halflength+1.01*Δx; # 1.01 is just an abitrary factor chosen to be greater than 1
xlim, ylim = (-domain_lim,domain_lim),(-domain_lim,domain_lim)

#= 
##Cavity Geometry
A square cavity can be created using the $rectangle()$ function with the half length defined above.
=#
body = Rectangle(halflength,halflength,1.5*Δx) 
plot(body)

#=
##Boundary Condition at the moving wall
Assign velocity to the top boundary.

The $LidDrivenCavity()$ function can be used to specify the velocity value at the top wall.

Note : Non-dimensional velocity = 1
=#
m = ViscousFlow.LidDrivenCavity(1.0); # motion type

#=
##Construct the system structure
The last two input "flow_side" and "static_points" must specified so the default setting in the $NavierStokes()$ function can be overwritten.

"static_points" is true because the cavity is not moving.
=#
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,m,flow_side = InternalFlow,static_points = true)

#=
Initialize
=#
u0 = newstate(sys)

#=
Set up integrator
=#
tspan = (0.0,5.0)
integrator = init(u0,tspan,sys)

#=
## Solve
=#
step!(integrator,5)

#=
## Examine
### RE=100, ΔX=0.008
plot for vorticity and streamlines
=#

plot(
plot(vorticity(integrator),sys,title="Vorticity (Computed)",clim=(-5,5),color=:turbo,linewidth=1.5,ylim=ylim,fillrange=nothing),
plot(streamfunction(integrator),sys,title="Streamlines (Computed)",size=(700,300))
   )

sol = integrator.sol;
@gif for (u,t) in zip(sol.u,sol.t)
    plot(vorticity(u,sys,t),sys,clim=(-10,10),levels=range(-10,10,length=30),color=:turbo,fillrange=nothing)
end every 5
