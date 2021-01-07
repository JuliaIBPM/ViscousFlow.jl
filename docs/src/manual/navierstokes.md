# Navier-Stokes systems

```@meta
DocTestSetup = quote
using ViscousFlow
end
```

```math
\def\ddt#1{\frac{\mathrm{d}#1}{\mathrm{d}t}}

\renewcommand{\vec}{\boldsymbol}
\newcommand{\uvec}[1]{\vec{\hat{#1}}}
\newcommand{\utangent}{\uvec{\tau}}
\newcommand{\unormal}{\uvec{n}}

\renewcommand{\d}{\,\mathrm{d}}
```


```@setup create
using ViscousFlow
using Plots
```

Here, we will focus on putting tools together from the previous sections in order to set up and solve the Navier-Stokes system of equations. First, we will solve them in a completely unbounded domain (i.e., no bodies), and then we will solve them in the vicinity of a body.

To perform any simulation in `ViscousFlow`, we need to carry out a few basic steps:
* **Specify the problem**: Set the Reynolds number and free stream
* **Discretize**: Set up a solution domain, grid cell size, time step size
* **Create bodies**: If present, set up bodies and motions
* **Construct the system structure**: Create the operators that will be used to perform the simulation
* **Initialize**: Set the initial flow field and initialize the integrator
* **Solve**: Solve the flow field
* **Examine**: Examine the results

## Navier-Stokes without a body

The first thing we must do is specify the problem parameters: here, this simply means to set the Reynolds number. We will set the other parameters (the initial conditions) later

```@repl corotate
Re = 200
```

We will also set up two Gaussian vortices, centered at $(x_{01},y_{01} = (0.5,0)$
and $(x_{02},y_{02} = (-0.5,0)$, with strength $\Gamma = 1$ and radius $\sigma = 0.1$.
For this, we use the function `SpatialGaussian`:

```@repl corotate
σ = 0.1
x01, y01 = 0.5, 0.0
x02, y02 = -0.5, 0.0
Γ = 1
twogauss = SpatialGaussian(σ,x01,y01,Γ) + SpatialGaussian(σ,x02,y02,Γ)
```

Now we will discretize the problem: set the domain, the grid spacing, and time step size.
For these latter two, we use `setstepsizes`. There are some keyword arguments for this
that might come in handy. Below, we use of the them, `gridRe`, to set the grid
Reynolds number ($U_c\Delta x/\nu$), for characteristic velocity $U_c$ and
kinematic viscosity $\nu$, to 4. If we do not set this explicitly, it defaults to 2.
Other keyword arguments are `fourier` (to change the grid Fourier number $\nu\Delta t/\Delta x^2$
  from its default value of 0.5) and `cfl` (to change the grid CFL number $U_c\Delta t/\Delta x$
  from the default 0.5).

```@repl corotate
xlim = (-2.0,2.0)
ylim = (-2.0,2.0)
Δx, Δt = setstepsizes(Re,gridRe=4)
```

Now we set up the Navier-Stokes system. This sets the rest of the grid parameters,
(number of cells, immersed boundary operators), and creates some cache space on the grid.

```@repl corotate
sys = NavierStokes(Re,Δx,xlim,ylim,Δt)
```

For example, to check how many dual grid cells we have, we can use the `size` function, which has been extended to such systems:

```@repl corotate
size(sys)
```

Let's set up the initial state vector for this system. We initialize it with the
pair of gaussians `twogauss` that we set up earlier:

```@repl corotate
u0 = newstate(twogauss,sys);
```

Next we initialize the integrator. The integrator advances the solution and holds all
of the solution history, in the same manner as the
[`OrdinaryDiffEq`](https://github.com/SciML/OrdinaryDiffEq.jl) package. We will set a
time span of 8 time units:

```@repl corotate
tspan = (0.0,8.0);
integrator = init(u0,tspan,sys)
```

Now we are ready to solve the problem. We advance it by a desired amount of time
with the `step!` function. Let's do it a little at a time:

```@repl corotate
step!(integrator,2.0)
```

Let's look at the vorticity field at this point:

```@repl corotate
plot(vorticity(integrator),sys,levels=range(0.1,5,length=31),title="Vorticity at t = 2")
savefig("w0corotate.svg"); nothing # hide
```
![](w0corotate.svg)

They have orbited around each other a bit. Let's advance further. All we need
to do is run `step!` again with a desired new time interval. It automatically
builds on the previous results:

```@repl corotate
step!(integrator,6.0)
```

```@repl corotate
plot(vorticity(integrator),sys,levels=range(0.1,5,length=31),title="Vorticity at t = 8")
savefig("w1corotate.svg"); nothing # hide
```
![](w1corotate.svg)

## A problem with a stationary body

Now let's solve for flow past a body. The only additional step compared to the
previous example is creating a body and placing it in the desired location. Here, we will
solve for flow past an ellipse at 45 degrees at Re = 200.

```@setup cylflow
using ViscousFlow
using Plots
```

Set the Reynolds number and free stream velocity:
```@repl cylflow
Re = 200
U∞ = (1.0,0.0);
```

Now discretize as before: set the domain, and grid spacing and time step:
```@repl cylflow
xlim = (-1.0,3.0);
ylim = (-1.5,1.5);
Δx, Δt = setstepsizes(Re,gridRe=4.0);
```

We will place the ellipse at the origin. Placement and configuration is done with the `RigidTransform`
function:

```@repl cylflow
body = Ellipse(0.5,0.1,1.5Δx)
cent = (0.0,0.0) # center of body
α = -45π/180 # angle, in radians
T! = RigidTransform(cent,α)
T!(body) # transform the body to the desired configuration
```


Now we set up the system structure, supplying `body`. The freestream is supplied as
a keyword argument:

```@repl cylflow
sys = NavierStokes(Re,Δx,xlim,ylim,Δt,body,freestream = U∞)
```

The rest of the setup is done just as before, except that the initial condition
is set to zero:

```
@repl cylflow
u0 = newstate(sys);
tspan = (0.0,20.0);
integrator = init(u0,tspan,sys);
```

Advance by 1 time unit:
```
step!(integrator,1.0)
```

Let's plot this:

```@repl cylflow
plot(vorticity(integrator),sys,title="Vorticity at t = 1",clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu,ylim=ylim)
savefig("ellipse0.svg"); nothing # hide
```
![](ellipse0.svg)

Now let's advance further:

```@repl cylflow
step!(integrator,29.0)
```

and plot it now:

```@repl cylflow
plot(vorticity(integrator),sys,title="Vorticity at t = 11",clim=(-10,10),levels=range(-10,10,length=30), color = :RdBu,ylim=ylim)
savefig("ellipse1.svg"); nothing # hide
```
![](ellipse1.svg)


We can also compute the force history:

```@repl cylflow
fx, fy = force(sol,sys,1);
```

The last argument (1) specifies that we are obtaining the force components on body 1. Now
plot them:

```@repl cylflow
plot(
plot(sol.t,2*fx,xlim=(0,Inf),ylim=(0,6),xlabel="Convective time",ylabel="\$C_D\$",legend=:false),
plot(sol.t,2*fy,xlim=(0,Inf),ylim=(-6,6),xlabel="Convective time",ylabel="\$C_L\$",legend=:false),
    size=(800,350)
    savefig("force.svg"); nothing # hide
)
```
![](force.svg)


## Methods

```@autodocs
Modules = [ViscousFlow]
```


## Index

```@index
Pages = ["navierstokes.md"]
```
