# Navier-Stokes systems

```@meta
CurrentModule = ViscousFlow.Systems
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

## Navier-Stokes without a body

Here, we seek the solve the two-dimensional incompressible Navier-Stokes equations in their *discrete vorticity form*, in an unbounded domain:

$$\ddt w + N(v,w) = \frac{1}{Re} L w,$$

along with the initial condition

$$w(0) = w_0.$$

The field $w$ represents the discrete vorticity, which sits at the nodes of the dual cells. The velocity, $v$, lies on the edges of the primal cells. They are related
to each other by $v = Cs$, where $s = -L^{-1} w$ is the discrete streamfunction.

The second term on the left-hand side is the convective term, which we have
simply written as $N(v,w)$. There are several ways to write this term; here, we
will write it by using the discrete divergence,

$$N(v,w) = D(vw).$$

The `Systems` module has a function that is set up to compute this term; we will discuss it below. The right-hand side contains the viscous term, proportional to $1/Re$, where $Re$ is the Reynolds number. For this, we will use the integrating factor, described in [The integrating factor](@ref). For purposes of calculation, it is better to express the problem as

$$\ddt w - \frac{1}{Re} L w = r_1(w),$$

where $r_1(w) = -D(vw)$.

For demonstration, we will solve a problem consisting initially of two identical circular patches of vorticity.

```@setup corotate
using ViscousFlow
using Plots
pyplot()
```

The first thing we must do is set up a grid. We will make it square, with spacing equal to 0.02 in each cell.

```@repl corotate
xlim = (-2,2); ylim = (-2,2);
Δx = 0.02;
```

Now we will set the Reynolds number, and set the time step size so that it follows the so-called *CFL* condition (with CFL number set to 0.5). To be careful, we also make sure the time step size does not exceed a threshold in the grid Fourier number (also set to 0.5):

```@repl corotate
Re = 200
Δt = min(0.5*Δx,0.5*Δx^2*Re)
```

Now we set up the Navier-Stokes system. This sets the rest of the grid parameters,
(number of cells, etc), and creates some some buffer space on the grid.

```@repl corotate
sys = NavierStokes(Re,Δx,xlim,ylim,Δt)
```

For example, to check how many dual grid cells we have, we can use the `size` function, which has been extended to such systems:

```@repl corotate
size(sys)
```

Let's set up a set of dual nodes on this grid:

```@repl corotate
w₀ = Nodes(Dual,size(sys));
```

The physical grid coordinates of these dual nodes can be generated with the `coordinates` function:

```@repl corotate
xg, yg = coordinates(w₀,dx=sys.Δx,I0=Systems.origin(sys))
```

Now we are ready to set up the integrator for this problem. To account for the viscous diffusion, we need the integrating factor. There are no body constraints to enforce, so we will use the integrating factor Runge-Kutta method (`IFRK`). For this, we need to set up plans for the integrating factor and for the right-hand side ($r_1$). The `Systems` module has functions that do both for us, using the system data in `sys`. We just need to change their argument list so that they fit the template for the `IFRK` scheme:

```@repl corotate
plan_intfact(t,w) = Systems.plan_intfact(t,w,sys)
r₁(w,t) = Systems.r₁(w,t,sys)
```

Now we can construct the integrator. We will use 3rd-order Runge-Kutta:

```@repl corotate
ifrk = IFRK(w₀,sys.Δt,plan_intfact,r₁,rk=TimeMarching.RK31)
```

Note that we have only passed in `w₀` to this scheme to provide the form of data to be used for the state vector in the integrator. It does not matter that the data are still zeros.

Finally we are ready to solve the problem. We set up the initial condition. It is helpful to define a function first that specifies the vorticity distribution in each vortex patch. We will use a Gaussian:

```@repl corotate
using LinearAlgebra
gaussian(x,x0,σ) = exp(-LinearAlgebra.norm(x.-x0)^2/σ^2)/(π*σ^2)
```

Now the initial conditions. We will put one vortex at $(-0.5,0)$ and the other at $(0.5,0)$. They will each have a strength of $1$ and a radius of $0.2$:

```@repl corotate
t = 0.0
x01 = (-0.5,0); x02 = (0.5,0); σ = 0.2; Γ = 1
w₀ .= [Γ*gaussian((x,y),x01,σ) + Γ*gaussian((x,y),x02,σ) for x in xg, y in yg]*Δx;
w = deepcopy(w₀);
```

Note that we have multiplied the vorticity vector by the grid spacing. This is because the vector `w` is not actually the vorticity, but rather, a *grid* vorticity related to velocity through differencing. Let's plot it to see what we are starting with:

```@repl corotate
plot(xg,yg,w)
savefig("w0corotate.svg"); nothing # hide
```
![](w0corotate.svg)


We will integrate the problem for 1 time unit:

```@repl corotate
tf = 1
T = 0:Δt:tf
```

Now, do it. We will time it to see how long it takes:

```@repl corotate
@time for ti in T
    global t, w = ifrk(t,w)
end
```

and plot it again:

```@repl corotate
plot(xg,yg,w)
savefig("w1corotate.svg"); nothing # hide
```
![](w1corotate.svg)

Let's go further!

```@repl corotate
tf = 5
T = 0:Δt:tf
@time for ti in T
    global t, w = ifrk(t,w)
end
```

```@repl corotate
plot(xg,yg,w)
savefig("w2corotate.svg"); nothing # hide
```
![](w2corotate.svg)

## Methods

```@autodocs
Modules = [Systems]
```


## Index

```@index
Pages = ["navierstokes.md"]
```
