#=
# 3. Applying pulse forcing to a flow
Here, we will show how to apply a transient force (a 'pulse') to a flow.
=#

using ViscousFlow
#-
using Plots

#=
For this case, we first seek to introduce a pulse to an otherwise quiescent fluid,
with no boundaries. The pulse shape will be a smooth dipole, directed upward.

We will start with the usual steps: specify the problem parameters and the discretization.
=#

Re = 200
xlim, ylim = (-2,2), (-2,4)
Δx,Δt = setstepsizes(Re,gridRe=4)

#=
### Construct the pulse
Now create the pulse. For the spatial shape of this pulse, we make use of a Gaussian
field, differentiated once in the `x` direction. We will center it at (0,0), and
give it a strength of 10.
=#
σx = 0.5
σy = 0.1
gauss = SpatialGaussian(σx,σy,0,0,10,deriv=1);

#=
Now we now 'shape' the pulse in time, also with a Gaussian. We center it at `t = 0.1`
with a half-width in time equal to 0.1.
=#
t0 = 0.1
σt = 0.1
pparams = PulseParams(gauss,t0,σt);

#=
### Construct the system structure
We supply these pulse characteristics with the keyword argument `pulses` as we
set up the usual system structure:
=#
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,pulses=pparams)

#=
The remaining steps go just as they did for the previous examples. We will simulate
the pulse for 5 time units and make an animation of the results:
=#
u0 = newstate(sys)
tspan = (0.0,10.0) # longer than we need, but just being safe
integrator = init(u0,tspan,sys)
#-
step!(integrator,5.0)

#=
### Examine
=#
sol = integrator.sol
@gif for (u,t) in zip(sol.u,sol.t)
    plot(vorticity(u,sys,t),sys)
end every 5
#-
plot(streamfunction(integrator),sys,title="Streamfunction at t = $(round(integrator.t,digits=2))")

#=
We can also supply more than one pulse. Let's give 3, each separated by 1 time unit.
=#
pparams = [PulseParams(gauss,0.1,0.1), PulseParams(gauss,1.1,0.1), PulseParams(gauss,2.1,0.1)];
#-
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,pulses=pparams)
#-
u0 = newstate(sys)
tspan = (0.0,10.0)
integrator = init(u0,tspan,sys)
#-
step!(integrator,4.0)
#=
### Examine
=#
sol = integrator.sol
@gif for (u,t) in zip(sol.u,sol.t)
    plot(vorticity(u,sys,t),sys)
end every 5
