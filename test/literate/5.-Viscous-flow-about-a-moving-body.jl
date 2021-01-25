#=
# 5. Viscous flow about a moving body
In this notebook we will demonstrate the simulation of a moving body. It is straightforward
to set up a moving body. The main caveat is that the simulation is considerably slower,
because the integrator must update the operators continuously throughout the simulation.

We will demonstrate this on an oscillating flat plate.
=#

using ViscousFlow
#-
using Plots

#=
### Problem specification and discretization
=#
Re = 200; # Reynolds number
U = 1.0; # Free stream velocity
U∞ = (U,0.0);
#-
xlim = (-1.0,3.0)
ylim = (-1.0,1.0);
Δx, Δt = setstepsizes(Re,gridRe=4.0)

#=
### Set up body
Set up the plate and place it at the origin
=#
body = Plate(1.0,1.0Δx)
T = RigidTransform((0.,0.),0.)
T(body)

#=
### Set the body motion
Now we specify the body motion. We will use oscillatory pitch-heave kinematics for this:
=#
a = 0.25 # location of pitch axis, a = 0.5 is leading edge
ϕp = -π/2  # phase lag of pitch
ϕh = 0.0  # phase lag of heave
A = 0.25  # amplitude/chord
fstar = 1/π # fc/U
α₀ = 0 # mean angle of attack
Δα = 10π/180 # amplitude of pitching
U₀ = 0.0 # translational motion (set to zero in place of free stream)
K = π*fstar # reduced frequency, K = πfc/U

oscil1 = RigidBodyTools.PitchHeave(U₀,a,K,ϕp,α₀,Δα,A,ϕh)
motion = RigidBodyMotion(oscil1)

# We can inspect the kinematics in this `motion` by plotting them:
plot(motion)

#=
### Construct the system structure
Here, we supply the motion as an another argument.
=#
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,motion,freestream = U∞)
#-
u0 = newstate(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)

#=
### Solve
This takes longer than it does for stationary bodies. Here, we only run it
for a little while just to demonstrate it.
=#
step!(integrator,0.1)

# ### Examine the solution
plot(vorticity(integrator),sys,clim=(-10,10),levels=range(-10,10,length=30),color = :RdBu)

# and the forces
sol = integrator.sol
fx, fy = force(sol,sys,1);
#-
plot(
plot(sol.t,2*fx,xlim=(0,Inf),ylim=(-3,3),xlabel="Convective time",ylabel="\$C_D\$",legend=:false),
plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-6,6),xlabel="Convective time",ylabel="\$C_L\$",legend=:false),
    size=(800,350)
)
