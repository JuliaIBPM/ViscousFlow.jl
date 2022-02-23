#=
# 4. Multiple stationary bodies
Adding multiple bodies to a problem is easy, using the concepts of a `BodyList`
and `RigidTransformList`.
=#
using ViscousFlow
#-
using Plots

#=
In this example, we will set up a problem with three cylinders arranged in a
formation in a free stream.
=#
my_params = Dict()
my_params["Re"] = 200
my_params["freestream speed"] = 1.0
my_params["freestream angle"] = 0.0

#=
Set up the domain and surface point spacing
=#
xlim = (-2.0,4.0)
ylim = (-2.0,2.0)
my_params["grid Re"] = 4.0
g = setup_grid(xlim,ylim,my_params)
Δs = surface_point_spacing(g,my_params)

#=
### Set up bodies
We start by initializing a `BodyList` and an associated `RigidTransformList`.
Each member of the `RigidTransformList` will be used to place the respective body in
the correct position and orientation.
=#
bl = BodyList()
tl = RigidTransformList()

# Place the first body at (-1,0)
push!(bl,Circle(0.5,Δs))
push!(tl,RigidTransform((-1,0),0));

# Place the second body at (1,-1)
push!(bl,Circle(0.5,Δs))
push!(tl,RigidTransform((1,-1),0));

# and place the third body at (1,1)
push!(bl,Circle(0.5,Δs))
push!(tl,RigidTransform((1,1),0))

# Perform the actual transformation. Note that this operation works `in-place`:
tl(bl)

# #### Plot the initial configuration of the bodies
# Just to check they are in the right places
plot(bl,xlim=xlim,ylim=ylim)

#=
### Construct the system structure
We construct the system with the same syntax as for a single body:
=#
sys = viscousflow_system(g,bl,phys_params=my_params);

#-
u0 = init_sol(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)
#=
### Solve
Here, we run it for a little while, just to demonstrate:
=#
@time step!(integrator,3.0)
#=
### Examine
Let's inspect the results
=#
sol = integrator.sol
plt = plot(layout = (1,3), size = (800, 200), legend=:false)
tsnap = 1.0:1.0:3.0
for (i, t) in enumerate(tsnap)
    plot!(plt[i],vorticity(sol,sys,t),sys,clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu)
end
plt

# Now we will examine the force on each body
fx1, fy1 = force(sol,sys,1)
fx2, fy2 = force(sol,sys,2)
fx3, fy3 = force(sol,sys,3);
#-
plt = plot(layout = (2,1), size = (600, 400))
plot!(plt[1],sol.t,2*fx1,xlim=(0,3),ylim=(0,4),xlabel="Convective time",ylabel="\$C_D\$",label="Lead body",title="Drag force")
plot!(plt[2],sol.t,2*fy1,xlim=(0,3),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_L\$",label="Lead body",title="Side force")
plot!(plt[1],sol.t,2*fx2,xlim=(0,3),ylim=(0,4),xlabel="Convective time",ylabel="\$C_D\$",label="Trailing body",title="Drag force")
plot!(plt[2],sol.t,2*fy2,xlim=(0,3),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_L\$",label="Trailing body",title="Side force")
#-
println("Mean drag coefficient on lead body = ", GridUtilities.mean(2*fx1[3:end]))
#-
println("Mean drag coefficient on trailing body = ", GridUtilities.mean(2*fx2[3:end]))
