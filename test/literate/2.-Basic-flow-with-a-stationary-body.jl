#=
# 2. Basic flow with a stationary body
In this notebook we will simulate the flow past a stationary body.
=#

using ViscousFlow
#-
using Plots

#=
### The basic steps
From the previous notebook, we add one additional step:
* **Specify the problem**: Set the Reynolds number and free stream
* **Discretize**: Set up a solution domain, grid cell size, time step size
* **Set up bodies**: *Create the body or bodies and specify their motions, if any*
* **Construct the system structure**: Create the operators that will be used to perform the simulation
* **Initialize**: Set the initial flow field and initialize the integrator
* **Solve**: Solve the flow field
* **Examine**: Examine the results
We will go through all of these here. For the examples we will carry out in this notebook,
the first three steps need only be carried out once.
=#

# ### Problem specification
# Set the Reynolds number and free stream
Re = 200 # Reynolds number
U = 1.0 # Free stream velocity
U∞ = (U,0.0);

#=
### Discretize
We set the grid Re to 4.0 here to get a quicker solution, though it is generally
better to make this smaller (it defaults to 2.0).
=#
xlim = (-1.0,3.0)
ylim = (-1.5,1.5)
Δx, Δt = setstepsizes(Re,gridRe=4.0)

#=
### Set up bodies
Here, we will set up a rectangle of half-height 0.5 and half-width 0.25
at 45 degrees angle of attack
=#
body = Rectangle(0.5,0.25,1.5Δx)

#=
We place the body at a desired location and orientation with the `RigidTransform`
function. This function creates an operator `T` that acts in-place on the body:
after the operation is applied, `body` is transformed to the correct location/orientation.
=#
cent = (0.0,0.0) # center of body
α = 45π/180 # angle
T = RigidTransform(cent,α)
T(body) # transform the body to the current configuration

# Let's plot it just to make sure
plot(body,xlim=xlim,ylim=ylim)

#=
### Construct the system structure
This step is like the previous notebook, but now we also provide the body and
the freestream:
=#
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,freestream = U∞)

#=
### Initialize
Now, we initialize with zero vorticity. Note that we do this by calling
`newstate` with no argument except for `sys` itself.
=#
u0 = newstate(sys)

#=
and now create the integrator, with a long enough time span to hold the whole
solution history:
=#
tspan = (0.0,20.0)
integrator = init(u0,tspan,sys)

#=
### Solve
Now we are ready to solve the problem. Let's advance the solution to $t = 1$.
=#
@time step!(integrator,1.0)

#=
### Examine
Let's look at the flow field at the end of this interval
=#
plot(
plot(vorticity(integrator),sys,title="Vorticity",clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu,ylim=ylim),
plot(streamfunction(integrator),sys,title="Streamlines",ylim=ylim,color = :Black),
    size=(700,300)
    )

#=
Now let's make a movie, like we did last time.
=#

sol = integrator.sol;
@gif for (u,t) in zip(sol.u,sol.t)
    plot(vorticity(u,sys,t),sys,clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu)
end every 5

#=
#### Compute the force history
To do this, we supply the solution history `sol`, the system `sys`, and the index
of the body (1).
=#
fx, fy = force(sol,sys,1);

#=
Plot the histories. Note that we are actually plotting the drag and lift
coefficient histories here:
$$ C_D = \dfrac{F_x}{\frac{1}{2}\rho U_\infty^2 L}, \quad C_L = \dfrac{F_y}{\frac{1}{2}\rho U_\infty^2 L} $$
Since the quantities in this simulation are already scaled by $\rho$, $U_\infty$, and $L$
(because $\rho$ has been scaled out of the equations, and the free stream speed is
set to 1 and the height of the shape to 1), then we obtain these coefficients by
simply dividing by 1/2, or equivalently, by multiplying by 2:
=#
plot(
plot(sol.t,2*fx,xlim=(0,Inf),ylim=(0,6),xlabel="Convective time",ylabel="\$C_D\$",legend=:false),
plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-6,6),xlabel="Convective time",ylabel="\$C_L\$",legend=:false),
    size=(800,350)
)

# The mean drag and lift coefficients are
meanCD = GridUtilities.mean(2*fx)
#-
meanCL = GridUtilities.mean(2*fy)
