# Navier-Stokes systems

```@meta
CurrentModule = Systems
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

The right-hand side contains the viscous term, proportional to $1/Re$, where $Re$ is the Reynolds number. Test.

## Methods

```@autodocs
Modules = [Systems]
```

```@docs
NavierStokes{NX,NY,N,isstatic}
```

## Index

```@index
Pages = ["navierstokes.md"]
```
