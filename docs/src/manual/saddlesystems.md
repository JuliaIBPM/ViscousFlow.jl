# Saddle point systems

```@meta
DocTestSetup = quote
using Whirl
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
using Whirl
using Plots
```
Saddle systems comprise an important part of solving mechanics problems with
constraints. In such problems, there is an underlying system to solve, and the
addition of constraints requires that the system is subjected to additional
forces (constraint forces, or Lagrange multipliers) that enforce these constraints
in the system. Examples of such constrained systems are the divergence-free
velocity constraint in incompressible flow (for which pressure is the associated
Lagrange multiplier field), the no-slip and/or no-flow-through condition in
general fluid systems adjacent to impenetrable bodies, and joint constraints in
rigid-body mechanics.

A general saddle-point system has the form

$$\left[ \begin{array}{cc} A & B_1^T \\ B_2 & 0\end{array}\right] \left(\begin{array}{c}u\\f \end{array}\right) = \left(\begin{array}{c}r_1\\r_2 \end{array}\right)$$

We are primarily interested in cases when the operator $A$ is symmetric and positive definite,
which is fairly typical. It is also fairly common for $B_1 = B_2$, so that the
whole system is symmetric.

`whirl` allows us to solve such systems for $u$ and $f$ in a fairly easy way.
We need only to provide rules for how to evaluate the actions of the various
operators in the system. Let us use an example to show how this can be done.


## A translating cylinder in potential flow

In irrotational, incompressible flow, the streamfunction $\psi$ satisfies Laplace's equation,

$$\nabla^2 \psi = 0$$

On the surface of an impenetrable body, the streamfunction must obey the constraint

$$\psi = \psi_b$$

where $\psi_b$ is the streamfunction associated with the body's motion. Let us
suppose the body is moving vertically with velocity 1. Then $\psi_b = -x$ for all
points inside or on the surface of the body. Thus, the streamfunction field outside
this body is governed by Laplace's equation subject to the constraint.

Let us solve this problem on a staggered grid, using the tools discussed in
the Fields section, including the regularization and interpolation methods to
immerse the body shape on the grid. Then our saddle-point system has the form

$$\left[ \begin{array}{cc} L & H \\ E & 0\end{array}\right] \left(\begin{array}{c}\psi\\f \end{array}\right) = \left(\begin{array}{c}0\\\psi_b \end{array}\right)$$

where $L$ is the discrete Laplacian, $H$ is the regularization operator, and
$E$ is the interpolation operator.

Physically, $f$ isn't really a force here, but
rather, represents the strengths of distributed singularities on the surface.
In fact, this strength represents the jump in normal derivative of $\psi$ across
the surface. Since this normal derivative is equivalent to the tangential velocity,
$f$ is the strength of the bound vortex sheet on the surface. This will be useful
to know when we check the value of $f$ obtained in our solution.

First, let us set up the body, centered at $(1,1)$ and of radius $1/2$. We will
also initialize a data structure for the force:

```@setup saddle
using Whirl
using Plots
pyplot()
```

```@repl saddle
n = 128; θ = linspace(0,2π,n+1);
xb = 1.0 + 0.5*cos.(θ[1:n]); yb = 1.0 + 0.5*sin.(θ[1:n]);
X = VectorData(xb,yb);
f = ScalarData(X);
```

Now let's set up a grid of size $102\times 102$ (including the usual layer
of ghost cells) and physical dimensions $2\times 2$.

```@repl saddle
nx = 102; ny = 102; Lx = 2.0; dx = Lx/(nx-2);
w = Nodes(Dual,(nx,ny));
```

We need to set up the operators now. First, the Laplacian:
```@repl saddle
L = plan_laplacian(size(w),with_inverse=true)
L⁻¹(w::T) where {T} = L\w
```
The last line just defines another operator for computing the inverse of $L$. We
have called it `L⁻¹` for useful shorthand.
This operator acts upon dual nodal data and returns data of the same type, e.g.
`ψ = L⁻¹(w)`. The saddle point system structure requires operators that have this
sort of form.

Now we need to set up the regularization `H` and interpolation `E` operators.
```@repl saddle
regop = Regularize(X,dx;issymmetric=true)
Hmat, Emat = RegularizationMatrix(regop,f,w);
H(f::ScalarData) = Hmat*f;
E(w::Nodes) = Emat*w;
```
The last two lines do as we did for the Laplacian: generate operators of the correct
form for the saddle point system structure.

Now we are ready to set up the system.
```@repl saddle
S = SaddleSystem((w,f),(L⁻¹,H,E),issymmetric=true,isposdef=true)
```
Note that we have provided a tuple of the types of data, `w` and `f`, that we want the solver
to work with, along with a tuple of the definitions of the three operators. A more
physically appealing (and easier to remember) way to supply the operators is to
supply them as a $2\times 2$ array:
```@repl saddle
Ã = [L⁻¹ H
     E   0];
S = SaddleSystem((w,f),Ã,issymmetric=true,isposdef=true)
```

We have also
set two optional flags, to specify that the system is symmetric and positive definite.
This instructs on which solver to use. (This is actually not quite true for the Laplacian: it is only positive semi-definite, since this operator has a null space. It is
  adequate criteria for using the conjugate gradient method, but we will have to
  be careful of some aspects of the solution, we we will see below.)

Let's solve the system. We need to supply the right-hand side.
```@repl saddle
w = Nodes(Dual,(nx,ny));
ψb = ScalarData(X);
ψb .= -(xb-1);
```
The right-hand side of the Laplace equation is zero. The right-hand side of the
constraint is the specified streamfunction on the body. Note that we have
subtracted the circle center from the $x$ positions on the body. The reason for
this will be discussed in a moment.

Finally, we solve it:
```@repl saddle
@time (ψ,f) = S\(w,ψb)
plot(ψ)
savefig("sfunc.svg"); nothing # hide
```
![](sfunc.svg)

The solution shows the streamlines for a circle in vertical motion, as expected.
All of the streamlines inside the circle are vertical.

## Methods

```@autodocs
Modules = [SaddlePointSystems]
Order   = [:type, :function]
```

## Index

```@index
Pages = ["saddlesystems.md"]
```
