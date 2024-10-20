#=
# Viscous flow about moving bodies
In this notebook we will demonstrate the simulation of a system of moving bodies.
As we will show, it is straightforward to set up a moving body, using the
tools in [RigidBodyTools.jl](https://github.com/JuliaIBPM/RigidBodyTools.jl).
The main caveat is that the simulation is slower,
because the integrator must update the operators continuously throughout the simulation.

We will demonstrate this on a system of three linked plates undergoing a flapping
motion, in which the middle plate heaves up and down, and the other two
bodies pitch back and forth on hinges connecting their edges to the middle plate.
=#

#md # ```@meta
#md # CurrentModule = ViscousFlow
#md # ```

using ViscousFlow
#-
#!jl using Plots

#=
### Problem specification and discretization
For simplicity, we will not create a free stream in this problem. Everything
here is the usual.
=#
my_params = Dict()
my_params["Re"] = 200
#-
xlim = (-1.0,1.0)
ylim = (-1.0,1.0)
my_params["grid Re"] = 4.0
g = setup_grid(xlim,ylim,my_params)

Δs = surface_point_spacing(g,my_params)

#=
### Set up body
Set up the plates.
=#
Lp = 0.5
body1 = Plate(Lp,Δs)
body2 = Plate(Lp,Δs)
body3 = Plate(Lp,Δs)
bl = BodyList([body1,body2,body3])

#=
### Set the body motions
Here, we make use of joints to prescribe the motion of every part of this system.
We will attach body 1 to the inertial system, making it oscillate up and down.
This is a special case of a joint with three degrees of freedom, called a `FreeJoint2d`.
Bodies 2 and 3 will each be connected by hinges (i.e., with a `RevoluteJoint`)
to body 1.

=#
parent_body, child_body = 0, 1
Xp = MotionTransform([0,0],0) # location of joint in inertial system
xpiv = [0,0] # place center of motion at center of the plate
Xc = MotionTransform(xpiv,0)

#=
Now the motion for joint 1, which we set up through the three degrees of freedom.
 The first and second
are meant to be fixed, so we give them zero velocity, and the third
we assign oscillatory kinematics
=#
adof = ConstantVelocityDOF(0.0)
xdof = ConstantVelocityDOF(0.0)

Ω = 1
A = 0.25  # amplitude/chord
ϕh = 0.0  # phase lag of heave
ydof = OscillatoryDOF(A,Ω,ϕh,0.0)

dofs = [adof,xdof,ydof]

#=
Now assemble the joint
=#
joint1 = Joint(FreeJoint2d,parent_body,Xp,child_body,Xc,dofs)

#=
Now the two hinges. Each of these is a `RevoluteJoint`. We place them
at either end of body 1. Joint 2 between bodies 1 and 2
=#
parent_body, child_body = 1, 2
Xp = MotionTransform([0.25,0],0) # right side of body 1
Xc = MotionTransform([-0.25,0],0) # left side of body 2
Δα = 20π/180 # amplitude of pitching
ϕp = π/2 # phase lead of pitch
θdof = OscillatoryDOF(Δα,Ω,ϕp,0.0)
joint2 = Joint(RevoluteJoint,parent_body,Xp,child_body,Xc,[θdof])


#=
and joint 3 between bodies 1 and 3
=#
parent_body, child_body = 1, 3
Xp = MotionTransform([-0.25,0],0) # left side of body 1
Xc = MotionTransform([0.25,0],0) # right side of body 3
ϕp = -π/2 # phase lead of pitch, but lagging rather than leading
θdof = OscillatoryDOF(Δα,Ω,ϕp,0.0)
joint3 = Joint(RevoluteJoint,parent_body,Xp,child_body,Xc,[θdof])

#=
Assemble everything together
=#
m = RigidBodyMotion([joint1,joint2,joint3],bl)

#=
We generate the initial joint state vector with `init_and update the body system and plot it
=#
x = init_motion_state(bl,m)
update_body!(bl,x,m)
#!jl plot(bl,xlim=xlim,ylim=ylim)

#=
Here is a useful macro to visualize the motion as a movie:

```
macro animate_motion(b,m,dt,tmax,xlim,ylim)
    return esc(quote
            bc = deepcopy($b)
            t0, x0 = 0.0, init_motion_state(bc,$m)
            dxdt = zero(x0)
            x = copy(x0)

            @gif for t in t0:$dt:t0+$tmax
                motion_rhs!(dxdt,x,($m,bc),t)
                global x += dxdt*$dt
                update_body!(bc,x,$m)
                plot(bc,xlim=$xlim,ylim=$ylim)
            end every 5
        end)
end
```
=#

#=
### Define the boundary condition functions
Instead of using the default boundary condition functions, we define
special ones here that provide the instantaneous surface velocity (i.e. the velocity
of every surface point) from the prescribed
motion. Every surface has an "exterior" and "interior" side. For
a flat plate, these two sides are the upper and lower sides, and both sides
are next to the fluid, so both of them are assigned the prescribed velocity
of the plate. (For closed bodies, we would assign this velocity to only
one of the sides, and zero to the other side. We will see an example of this in a later case.)
We pack these into a special dictionary and
pass these to the system construction.
=#
function my_vsplus(t,x,base_cache,phys_params,motions)
  vsplus = zeros_surface(base_cache)
  surface_velocity!(vsplus,x,base_cache,motions,t)
  return vsplus
end

function my_vsminus(t,x,base_cache,phys_params,motions)
  vsminus = zeros_surface(base_cache)
  surface_velocity!(vsminus,x,base_cache,motions,t)
  return vsminus
end

bcdict = Dict("exterior" => my_vsplus, "interior" => my_vsminus)

#=
### Construct the system structure
Here, we supply both the motion and boundary condition functions as additional arguments.
=#
sys = viscousflow_system(g,bl,phys_params=my_params,motions=m,bc=bcdict);

#=
and generate the initial condition
=#
u0 = init_sol(sys)

#=
Before we solve the problem, it is useful to note that the Reynolds number
we specified earlier may not be the most physically-meaningful Reynolds number.
More relevant in this problem is the Reynolds number based on the maximum
body speed and the total length of the plates
=#
Umax, imax, tmax, bmax = maxvelocity(u0,sys)
L = 3*Lp
Re_eff = my_params["Re"]*Umax*L

#=
In other problems, we have stored the inverse of the matrix system that we
have to solve at each time stage. But since that matrix system changes
as the body moves, it is faster here to use an iterative solver (congjugate gradient)
=#
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys,alg=LiskaIFHERK(saddlesolver=CG))

#=
### Solve
This takes a bit longer per time step than it does for stationary bodies. Here, we only
run it for 1.5 time units just to demonstrate it.
=#
@time step!(integrator,1.5)

#=
### Examine the solution
Let's look at a few snapshots of the vorticity field. Note that the
plotting here requires us to explicitly call the [`surfaces`](https://juliaibpm.github.io/ImmersedLayers.jl/stable/manual/problems/#ImmersedLayers.surfaces)
function to generate the instantaneous configuration of the plate.
=#
#!jl sol = integrator.sol
#!jl plt = plot(layout = (1,3), size = (800, 300), legend=:false)
#!jl tsnap = 0.5:0.5:1.5
#!jl for (i, t) in enumerate(tsnap)
#!jl     plot!(plt[i],vorticity(sol,sys,t),sys,layers=false,title="t = $(round(t,digits=2))",clim=(-5,5),levels=range(-5,5,length=30),color = :RdBu)
#!jl     plot!(plt[i],surfaces(sol,sys,t))
#!jl end
#!jl plt
# and the forces and moments
sol = integrator.sol
mom, fx, fy = force(sol,sys,1);
#-
#!jl plot(
#!jl plot(sol.t,2*fx,xlim=(0,Inf),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_D\$",legend=:false),
#!jl plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-2,2),xlabel="Convective time",ylabel="\$C_L\$",legend=:false),
#!jl     size=(800,350)
#!jl )
