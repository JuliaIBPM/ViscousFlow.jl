# Time marching

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
`ViscousFlow` is equipped with a few classes of time marching schemes for advancing time-dependent
equations.

## Integrating factor systems

Integrating factor systems that we encounter in `ViscousFlow` are of the form

$$\ddt u = A u + r_1(u,t), \quad u(0) = u_0$$

The operator $A$ may be a matrix or a scalar, but is generally independent of time. (The
  method of integrating factors can deal with time-dependent $A$, but we don't encounter
  such systems in the `ViscousFlow` context so we won't discuss them.) For this purpose, we use the `IFRK` class of solver, which stands for Integrating Factor Runge-Kutta. This method solves
  the part associated with $A$ exactly, via the integrating factor, and advances a modified
  equation by Runge-Kutta method to account for the remaining part $r_1$.

  We discussed the construction
  of the integrating factor in the context of fields in [Fields](@ref). But first, let's
  give an example of how we can solve a simpler problem with just a single scalar-valued
  $u$. The example we will solve is

$$\ddt u = -\alpha u + \cos(\omega t),\quad u(0) = u_0$$

The exact solution is easily obtained:

$$u(t) = u_0 e^{-\alpha t} + \frac{1}{\alpha^2+\omega^2} \left[ \alpha(\cos(\omega t) - e^{-\alpha t}) + \omega \sin (\omega t)\right]$$

Let's solve it numerically, so we can evaluate the accuracy of the solver. We should note that the
integrating factor for this system is $e^{-\alpha t}$.

For demonstration, we will set $\alpha = 1$, $\omega = 4$, and $u_0 = 1$.

```@setup march
using ViscousFlow
using Plots
pyplot()
```

```@repl march
α = 1; ω = 4; u₀ = 1.0;
```

Here is the exact solution for later comparison
```@repl march
uex(t) = u₀*exp(-α*t) + (α*(cos(ω*t)-exp(-α*t))+ω*sin(ω*t))/(α^2+ω^2)
```

The first steps are to define operators that provide the integrating factor and the right-hand side
of the equations. For the integrating factor, we extend the definition of [`plan_intfact`](@ref)
from [Fields](@ref).

```@repl march
ViscousFlow.plan_intfact(t::Float64,u::Vector{Float64}) = exp(-α*t);
```

Note that we have defined this extended form of `plan_intfact` to adhere to the standard form,
accepting arguments for time `t` and the state vector `u`, even though the state vector isn't strictly needed here. The state 'vector' in this problem is actually only a scalar, of course. But
the time marching method does not accept scalar-type states currently, so we will
make `u` a 1-element vector to use the `ViscousFlow` tools.

Now let us define the right-hand side function. This function should also adhere to the standard
form, which requires the state vector `u` and the time `t` as arguments.

```@repl march
r₁(u::Vector{Float64},t::Float64) = cos(ω*t);
```

We also need to set the time-step size ($0.01$) and the initial condition. For the latter,
we set up the state vector as a 1-element vector, as discussed earlier:
```@repl march
Δt = 0.01;
u = [u₀];
```
We can now construct the integrator. We supply a form of the state vector (for use as a template
  for pre-allocating space for internal storage variables), the time-step size, and the
  definitions of the integrating factor and the right-hand side function:

```@repl march
ifrk = IFRK(u,Δt,plan_intfact,r₁,rk=TimeMarching.RK31)
```

We have set the time step size to $0.01$. We have also specified that the Runge-Kutta method to be used is a third-order method, `RK31`, specially designed for storing as few different versions of the integrating factor as necessary. This is actually the default method, so we could have omitted this keyword
argument. There are other choices, as well, such as `TimeMarching.Euler` for the
forward Euler method.

Now we can solve the system. The integrator has a simple form, accepting as arguments
the current time and state, and returning the updated versions of these at the end of the
step. We place this integrator inside of a loop and store the results. (Since `u` is set up
  as a 1-element vector, then we will store only the element of this vector.)

```@repl march
uhist = Float64[]; # for storing the solution
T = 0:Δt:10;
t = 0.0;
for ti in T
  push!(uhist,u[1]) # storage
  global t, u = ifrk(t,u) # advancement by one step by the integrator
end
```  

Now we can plot the result and compare it with the exact solution.

```@repl march
plot(T,uhist,label="numerical",xlabel="t",ylabel="u(t)")
plot!(T,uex.(T),label="exact soln")
savefig("ifrk.svg"); nothing # hide
```
![](ifrk.svg)

As we can see, the results are nearly indistinguishable.

## Constrained systems

## Constrained integrating factor systems

Constrained integrating factor systems that we encounter in `ViscousFlow` are of the form

$$\ddt u = A u - B_1^T f + r_1(u,t), \quad B_2 u = r_2(u,t), \quad u(0) = u_0$$

where $f$ is again the Lagrange multiplier for enforcing the constraints on $u$. Now, we combine the ideas of the last two sections into a single integrator.

Let's demonstrate this on the example of heat diffusion from a circular ring whose temperature
is held constant. In this case, $A$ is the discrete Laplace operator, $L$, times the heat diffusivity,
$r_1$ is zero (in the absence of volumetric heating sources), and $r_2$ is the temperature of
the ring. The operators $B_1^T$ and $B_2$ will be the regularization and interpolation
operators between discrete point-wise data on the ring and the field data.

The ring will have radius $1/2$ and fixed temperature $1$, and
the heat diffusivity is $1$. (In other words, the problem has been non-dimensionalized
by the diameter of the circle, the dimensional ring temperature, and the dimensional diffusivity.)

First, we will construct a field to accept the temperature on

```@repl march
nx = 129; ny = 129; Lx = 2.0; Δx = Lx/(nx-2);
u₀ = Nodes(Dual,(nx,ny)); # field initial condition
```

Now set up a ring of points on the circle at center $(1,1)$.

```@repl march
n = 128; θ = range(0,stop=2π,length=n+1);
R = 0.5; xb = 1.0 .+ R*cos.(θ); yb = 1.0 .+ R*sin.(θ);
X = VectorData(xb[1:n],yb[1:n]);
f = ScalarData(X); # to be used as the Lagrange multiplier
```

From this, construct the regularization and interpolation operators in their usual
symmetric form, and then set up a routine that will provide these operators inside the integrator:

```@repl march
reg = Regularize(X,Δx;issymmetric=true)
Hmat, Emat = RegularizationMatrix(reg,f,u₀);
plan_constraints(u::Nodes{Dual,nx,ny},t::Float64) = Hmat, Emat
```

Now set up the right-hand side operators. Both must take the standard form, with
arguments of the types of `u` and `t`. For $r_1$, we will simply set it to a field
of zeros in the same type as `u`. For $r_2$, we set the result uniformly to $1$.

```@repl march
r₁(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = Nodes(T,u); # sets to zeros
r₂(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = 1.0; # sets uniformly to 1.0
```

We will set the time-step size to a large value ($1.0$) for demonstration purposes.
The method remains stable for any choice. We also initialize time `t` and the state
`u`:

```@repl march
Δt = 1.0;
t = 0.0;
u = deepcopy(u₀);
```

Now we can construct the integrator. We supply examples for the state `u` and the
Lagrange multiplier data `f`, the time-step size, the constructor for the
integrating factor, a tuple of the operators for computing the actions of $B_1^T$ and $B_2$
on data of type `f` and `u`, respectively (which, in this case, are matrices `Hmat` and `Emat`),
and a tuple of the right-hand side functions.

```@repl march
ifherk = IFHERK(u,f,Δt,plan_intfact,plan_constraints,(r₁,r₂),rk=TimeMarching.Euler)
```

Here we've set the method to forward Euler. The resulting integrator accepts
as arguments the current time `t` and the current state `u`, and returns the
time, state, and Lagrange multiplier data at the end of the time step.

Now, let's advance the system. We'll also time it.

```@repl march
@time for i = 1:20
  global t, u, f = ifherk(t,u)
end
```

Now let's plot it

```@repl march
xg, yg = coordinates(u,dx=Δx);
plot(xg,yg,u)
plot!(xb,yb,linecolor=:black,linewidth=1.5)
savefig("ifherk.svg"); nothing # hide
```
![](ifherk.svg)

From a side view, we can see that it enforces the boundary condition:

```@repl march
plot(xg,u[65,:],xlabel="x",ylabel="u(x,1)")
savefig("ifherk-side.svg"); nothing # hide
```
![](ifherk-side.svg)

## Methods

```@autodocs
Modules = [TimeMarching]
Order   = [:type, :function]
```

## Index

```@index
Pages = ["timemarching.md"]
```
