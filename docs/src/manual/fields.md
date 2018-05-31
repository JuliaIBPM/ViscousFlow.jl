# Fields

```@meta
DocTestSetup = quote
using Whirl
srand(1)
end
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
parameterized by `NX` and `NY`, these values do not correspond to the number of primal
edges in each direction on this dual grid.

Field operations transform one data type to another. Some of these are differencing
operations, analogous to differential counterparts in continuum calculus: `curl`,
`divergence`, and `gradient`. For example, a `curl` operation can act upon dual nodal data
(like streamfunction) and return primal edge data (i.e. velocity); a `divergence`
operation acts on edge data (primal or dual) and returns nodal data of the same cell
type. Note that these operations are *mimetic*: they maintain some of the same properties as the
continuous counterparts. For example, the divergence of the curl of any dual nodal
data is exactly zero. The curl of the gradient of primal nodal data is also zero.

There is also a Laplacian operator, `laplacian`, which acts upon data of any type
and returns the 5-point (in 2-d) discrete Laplacian. Importantly, this operator
also comes outfitted with an *inverse* operation, based on the lattice Green's function.
The inverse does not depend on the size of the grid.

Other field operations shift the data, by local averaging, from one data type to
another. These operations are all called `shift`.

## Methods

```@autodocs
Modules = [Fields]
Order   = [:type, :function]
```

## Index

```@index
Pages = ["fields.md"]
```
