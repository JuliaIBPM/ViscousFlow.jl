# Fields

```@meta
DocTestSetup = quote
using Whirl
srand(1)
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
pyplot()
```
In `Whirl`, field data, such as velocity, vorticity and pressure, are stored on
a staggered uniform grid. Such a grid is divided into *cells*, with *edges* (which,
on a two-dimensional grid, are the same as *faces*) and *nodes* (cell centers).
Nodes hold scalar-valued data. Edges, on the other hand, hold the components of
vector-valued data that are normal to the respective edges; one component lies
on the vertical edges, while the other is on the horizontal edges.

Furthermore, there are two different cell types: *primal* and *dual*. On
the physical grid, these cell types are offset with respect to each other by half
a cell spacing in each direction. In other words, the four corners of the primal
(resp. dual) cell are the nodes of four dual (resp. primary) cells.

Thus, on a two-dimensional staggered grid, there are four distinct vector spaces,
associated with where the data are held on the grid:
- dual nodes,
- dual edges,
- primal nodes, and
- primal edges.
In `Whirl`, these are each distinct data types. Furthermore, the relationships between these types are
defined by an underlying grid shared by all. By convention, this grid is defined by
the number of dual cells `NX` and `NY` in each direction; we will often refer to it
as the *dual grid*. For example, `Nodes{Dual,NX,NY}` is the type for dual node data
on this grid; `Edges{Primal,NX,NY}` is the type for edge data on the primal cells
within this same `NX` by `NY` dual grid. Note that, even though this latter type is
parameterized by `NX` and `NY`, these values do *not* correspond to the number of primal
edges in each direction on this dual grid. These values always correspond to the
number of dual cells on the grid, for any data type. This makes it clear the
grid is shared by all data.

## Setting up field data

Let's see an example of creating a blank set of dual node data and filling it with
something:

```@repl create
w = Nodes(Dual,(5,4))
w .= reshape(1:20,5,4)
```

Other data types on the same grid can be set up in similar fashion. To ensure
that they have a size that is consistent with the dual node data `w`, we can use
this in place of the size:
```@repl create
q = Edges(Primal,w);
q.u[2,3] = 1;
q
```

## Field differencing operations

Field operations transform one data type to another. Some of these are differencing
operations, analogous to differential counterparts in continuum calculus: `curl`,
`divergence`, and `gradient`. For example, a `curl` operation can act upon dual nodal data
(like streamfunction) and return primal edge data (i.e. velocity); a `divergence`
operation acts on edge data (primal or dual) and returns nodal data of the same cell
type. Note that these operations are *mimetic*: they maintain some of the same properties as the
continuous counterparts. For example, the divergence of the curl of any dual nodal
data is exactly zero. The curl of the gradient of primal nodal data is also zero.

Let's take the curl of the dual nodal data we constructed:
```@repl create
curl(w)
```

We could also make this a little more cute by giving the curl operator a symbol
and then acting upon the data as though it were a matrix-vector operation:
```@repl create
C = Curl()
C*w
```
Note that `C` is not actually a matrix. Rather, it is simply another name for the
`curl` operator, and `*` is defined in this context to apply `curl` to whatever is
to the right of it. The other operators have similar constructs.

Suppose we wish to apply the `curl` operation over and over. The `curl()` function
allocates memory for the result whenever it is used; this would become expensive
if it is done often. So it makes sense to preallocate space for this result and
use the `curl!()` function, which simply fills in the elements:
```@repl create
q = Edges(Primal,w)
curl!(q,w)
```
Note that we used a convenience function for setting up primal edge data `q` of a
size that corresponds with `w`.

Let's check that divergence of the curl is indeed zero:
```@repl create
D = Divergence()
D*(C*w)
```

## The Laplacian and its inverse

`Whirl` also makes heavy use of the discrete Laplacian operator, $L$. This mimics the
continuous operator, $\nabla^2$, and acts upon data of any type. Let's apply
this to the original data:
```@repl create
laplacian(w)
```

As with the other operators, we can also construct a shorthand of the discrete
Laplacian operator,
```@repl create
L = Laplacian(size(w))
L*w
```

An important part of `whirl` is the *inverse* of this operator. That is, we need
the ability to solve the discrete Poisson system

$$Ls = w$$

for $s$, for given data $w$. We achieve this in `whirl` with the *lattice Green's
function*. To outfit the operator with its inverse, we simply set the optional
flag:
```@repl create
L = Laplacian(size(w),with_inverse=true)
```

Then, the Poisson system is solved with the backslash (`\`),
```@repl create
s = L\w
L*s
```

It should be observed that the cells on the perimeter have not recovered the original values
of `w`. These are the ghost cells, and the Laplacian operation does not apply
to these.

It is also important to note that, although it looks as though we've constructed a
matrix `L` and performed various matrix-vector operations with it, this is not
actually the case. In fact, the `\`
operation associated with `L` is significantly faster than a matrix inversion.
Internally, it carries out a fast convolution between the data in `w` and the
lattice Green's function, via fast Fourier transform. The lattice Green's function
(LGF) table is pre-computed and pre-transformed in the original construction of `L`.
(In fact, because this table is not dependent on the size of the grid, it is
actually computed once for all time and stored in a file; subsequent applications
of it just load it in and use the portion of it necessary for a certain grid.)

The lattice Green's function has the advantage that it is independent of the grid
size. Let's solve the Poisson system when $w$ is a unit field, i.e. a field
of zeros, except for a single $1$ entry at one node. The solution $s$ represents
the influence of this point on all nodes. To see that the LGF does
not depend on the grid size, let's use a grid that is long and skinny and plot
the solution on it
```@repl create
w = Nodes(Dual,(50,10));
w[20,5] = 1.0
L = Laplacian(w,with_inverse=true)
plot(L\w)
savefig("Linvw.svg"); nothing # hide
```
![](Linvw.svg)

The influence is not affected by the narrow grid dimensions.

## The integrating factor

An operator related to the lattice Green's function is the *integrating factor*.
Suppose we have the system of ODEs

$$\ddt u = L u + f(u,t), \quad u(0) = u_0,$$

where $L$ is the discrete Laplacian (on an infinite uniform grid), and $u$ are
nodal data (and $f$ is a nodal-valued function acting on this nodal data). The
exact solution of this problem is

$$u(t) = E(t)u_0 + \int_0^t E(t-\tau) f(u(\tau),\tau)\,\mathrm{d}\tau,$$

where $E(t)$ is the integrating factor (or matrix exponential) for the system. The
easiest way to understand the role of $E(t)$ is to consider its behavior when $f$
is zero and $u_0$ contains a field of zeros except for a single $1$ entry at one
cell. Let's set up this initial data:
```@repl create
w = Nodes(Dual,(100,100));
w[40,50] = 1.0
plot(w)
savefig("w1.svg"); nothing # hide
```
![](w1.svg)

Then, $E(t)u_0$ diffuses this initial unit perturbation in each direction. Here, we apply it
with $t = 5$:

```@repl create
E = IntFact(5,w)
plot(E*w)
savefig("Ew1.svg"); nothing # hide
```
![](Ew1.svg)

Note that $E(0) = I$, where $I$ is the identity. Also, the integrating factor has the useful property that $E(t+\tau) = E(t)E(\tau)$. From these properties, it
follows that $E^{-1}(t) = E(-t)$. Let's suppose we wish to advance $u$ from time
$t = \tau-h$ to time $t = \tau$. Then we can define an auxiliary quantity,
$v(t;\tau) = E(\tau-t)u(t)$, and this new quantity satisfies the modified set
of ODEs

$$\ddt v = E(\tau-t) f\left[ E(t-\tau) v(t;\tau),t\right],\quad v(\tau-h;\tau) = E(h)u(\tau-h)$$    

The result of integrating this set of ODEs to $t = \tau$ is $v(\tau;\tau) = u(\tau)$. In
other words, the integrating factor allows us to solve a somewhat reduced set
of ODEs.

## Other field operations

Other field operations shift the data, by local averaging, from one data type to
another. These operations are all called `shift!`, and they require that the
target data be preallocated. For example, to shift dual node data to the dual edges,

```@repl create
w = Nodes(Dual,(5,4));
w .= reshape(1:20,5,4)
Ww = Edges(Dual,w);
shift!(Ww,w)
```
Note that the edges in the ghost cells are 0; these edges are not assigned any
values in the shift operation.

We can then shift this to primal edges:
```@repl create
q = Edges(Primal,w);
shift!(q,Ww)
```

We can also compute the Hadamard (i.e. element by element) product of any data
of the same type, e.g.,
```@repl create
q∘q
```


## Methods

```@autodocs
Modules = [Fields]
Order   = [:type, :function]
```

## Index

```@index
Pages = ["fields.md"]
```
