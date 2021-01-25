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
Re = 200 # Reynolds number
U = 1.0 # Free stream velocity
U∞ = (U,0.0);

#=
Set up the domain, grid spacing, and time step size
=#
xlim = (-2.0,4.0)
ylim = (-2.0,2.0)
Δx, Δt = setstepsizes(Re,gridRe=4)

#=
### Set up bodies
We start by initializing a `BodyList` and an associated `RigidTransformList`.
Each member of the `RigidTransformList` will be used to place the respective body in
the correct position and orientation.
=#
bl = BodyList()
tl = RigidTransformList()

# Place the first body at (-1,0)
push!(bl,Circle(0.5,1.5Δx))
push!(tl,RigidTransform((-1.,0.),0.));

# Place the second body at (1,-1)
push!(bl,Circle(0.5,1.5Δx))
push!(tl,RigidTransform((1.,-1.),0.));

# and place the third body at (1,1)
push!(bl,Circle(0.5,1.5Δx))
push!(tl,RigidTransform((1.,1.),0.))

# Perform the actual transformation. Note that this operation works `in-place`:
tl(bl)

# #### Plot the initial configuration of the bodies
# Just to check they are in the right places
plot(bl,xlim=xlim,ylim=ylim)

#=
### Construct the system structure
We construct the system with the same syntax as for a single body:
=#
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,bl,freestream = U∞)
#-
u0 = newstate(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)
#=
### Solve
Here, we run it for only a little while, just to demonstrate:
=#
@time step!(integrator,0.5)
#=
### Examine
Let's make an animation
=#
sol = integrator.sol;
@gif for (u,t) in zip(sol.u,sol.t)
    plot(vorticity(u,sys,t),sys,clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu)
end every 5

# Now we will examine the force on each body
fx1, fy1 = force(sol,sys,1)
fx2, fy2 = force(sol,sys,2)
fx3, fy3 = force(sol,sys,3);
#-
plt = plot(layout = (2,1), size = (600, 400))
plot!(plt[1],sol.t,2*fx1,xlim=(0,20),ylim=(0,4),xlabel="Convective time",ylabel="\$C_D\$",label="Lead body",title="Drag force")
plot!(plt[2],sol.t,2*fy1,xlim=(0,20),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_L\$",label="Lead body",title="Side force")
plot!(plt[1],sol.t,2*fx2,xlim=(0,20),ylim=(0,4),xlabel="Convective time",ylabel="\$C_D\$",label="Trailing body",title="Drag force")
plot!(plt[2],sol.t,2*fy2,xlim=(0,20),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_L\$",label="Trailing body",title="Side force")
#-
println("Mean drag coefficient on lead body = ", GridUtilities.mean(2*fx1))
#-
println("Mean drag coefficient on trailing body = ", GridUtilities.mean(2*fx2))
