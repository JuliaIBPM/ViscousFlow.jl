# # Exogenous degrees of freedom

#md # ```@meta
#md # CurrentModule = ViscousFlow
#md # ```

#=
Thus far, the motion of bodies or of incident flows has been prescribed
with time-varying functions known a priori. In other words, we knew
in advance how the body or free stream would behave for all time. However,
it is often the case that we do *not* know the motion in advance, because,
e.g., that motion depends on the evolving state of the fluid flow. That is
the case for fluid-body interactions, in which the body motion depends
on the forces imparted by the fluid on the body. It is also the
case for active feedback control or reinforcement learning, in which
a controller/agent might change some degree of freedom based on what
is observed.

We will refer to such a degree of freedom as *exogenous*. In this example,
we will demonstrate how the package can incorporate exogenous degrees
of freedom into the flow evolution. These exogenous states are specified by providing
the value of their acceleration.

The example will be a flat plate
in a steady free stream at nominally zero angle of incidence. However,
the y acceleration will vary randomly. Since there is only a single body
and it remains rigid, we will solve the problem in a reference frame moving
with the body. However, exogenous degrees of freedom can be used in problems
of arbitrary motion.

In this simple example, this exogenous
y acceleration will not depend on anything, but it will be obvious from the
example how it *could*.
=#
using ViscousFlow
#-
#!jl using Plots

#=
Set the Reynolds number
=#
my_params=Dict()
my_params["Re"] = 200

#=
Set up the grid
=#
xlim = (-2.0,4.0)
ylim = (-3.0,3.0)
my_params["grid Re"] = 3.0
g = setup_grid(xlim,ylim,my_params)

#=
Set up the body
=#
Δs = surface_point_spacing(g,my_params)
body = Plate(1.0,Δs)

#=
### Kinematics
There won't be any rotation, so we will put the
joint at the body's center.
=#
parent_body, child_body = 0, 1
Xp = MotionTransform([0,0],0) # transform from inertial system to joint
xpiv = [0.0,0.0]
Xc = MotionTransform(xpiv,0.0)

#=
Now set the kinematics. Rather than provide a free stream, we will set the x velocity of
the body to be -1. We will set the angular velocity to zero.
=#
adof = ConstantVelocityDOF(0)
xdof = ConstantVelocityDOF(-1)

#=
Now, to designate the y degree of freedom as exogenous, we simply use
=#
ydof = ExogenousDOF()

#=
Now assemble the joint and the motion
=#
dofs = [adof,xdof,ydof]
joint = Joint(FreeJoint2d,parent_body,Xp,child_body,Xc,dofs)
m = RigidBodyMotion(joint,body)

#=
Set the boundary condition functions
=#
function my_vsplus(t,x,base_cache,phys_params,motions)
  vsplus = zeros_surface(base_cache)
  surface_velocity_in_translating_frame!(vsplus,x,base_cache,motions,t)
  return vsplus
end

function my_vsminus(t,x,base_cache,phys_params,motions)
  vsminus = zeros_surface(base_cache)
  surface_velocity_in_translating_frame!(vsminus,x,base_cache,motions,t)
  return vsminus
end

bcdict = Dict("exterior" => my_vsplus, "interior" => my_vsminus)

#=
### Construct system and initialize
=#
sys = viscousflow_system(g,body,phys_params=my_params,bc=bcdict,motions=m,reference_body=1);

u0 = init_sol(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)
dt = timestep(u0,sys)

#=
### Solve
To solve the problem, we advance the solution inside a loop. The timestep
with which we advance is up to us, though it should be an integer multiple
of the timestep used by the underlying solver (`dt` above). This prevents
the need for expensive interpolations. We will simply set it to `dt` here.

Inside the loop, we update the value of the exogenous acceleration,
using the `update_exogenous!` function. This is done with a vector of
all such exogenous states, though there is only one state in this example.
We draw its value from a normal distribution.
=#
dt_advance = dt
tfinal = 2.5
for t in 0:dt_advance:tfinal
    a_y = randn()
    update_exogenous!(integrator,[a_y])
    step!(integrator,dt_advance)
end

#=
### Plot it
Plot the vorticity field
=#
sol = integrator.sol
#!jl plt = plot(layout = (4,3), size = (900, 900), legend=:false)
#!jl tsnap = tfinal/12:tfinal/12:tfinal
#!jl for (i, t) in enumerate(tsnap)
#!jl     plot!(plt[i],vorticity(sol,sys,t),sys,clim=(-20,20),levels=range(-20,20,length=16),ylim=(-1.5,1.5),title="$(round(t,digits=4))")
#!jl end
#!jl plt

#=
and compute and plot the force and moment
=#
mom, fx, fy = force(sol,sys,1)

#!jl plot(
#!jl plot(sol.t,2*fx,xlim=(0,Inf),ylim=(-5,5),xlabel="\$t\$",ylabel="\$C_x\$",legend=:false),
#!jl plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-5,5),xlabel="\$t\$",ylabel="\$C_y\$",legend=:false),
#!jl plot(sol.t,2*mom,xlim=(0,Inf),ylim=(-2,2),xlabel="\$t\$",ylabel="\$C_m\$",legend=:false),
#!jl     size=(800,350)
#!jl )

#=
We can also look at the y position and velocity. These are contained in
the joint state vector, `aux_state(u)`. We can use the functions `exogenous_position_vector`
and `exogenous_velocity_vector` to access them. (The index `[1]` below is required
because these functions return a vector of exogenous states for a particular joint.
There is only one component that is exogenous in this example.
=#
jointid = 1
y = map(u -> exogenous_position_vector(aux_state(u),m,jointid)[1],sol)
v = map(u -> exogenous_velocity_vector(aux_state(u),m,jointid)[1],sol)
#!jl plot(sol.t,y,label="\$y\$",xlabel="\$t\$")
#!jl plot!(sol.t,v,label="\$v\$")
