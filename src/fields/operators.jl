include("convolution.jl")
include("lgf.jl")
include("intfact.jl")
include("ddf.jl")

import Base: *, \, A_mul_B!, A_ldiv_B!

# laplacian

"""
    laplacian!(v,w)

Evaluate the discrete Laplacian of `w` and return it as `v`. The data `w` can be
of type dual/primary nodes or edges; `v` must be of the same type.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> v = deepcopy(w);

julia> w[4,3] = 1.0;

julia> laplacian!(v,w)
Whirl.Fields.Nodes{Whirl.Fields.Dual,8,6} data
Printing in grid orientation (lower left is (1,1)):
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   1.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  -4.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0   1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function laplacian!(out::Nodes{Dual,NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 2:NX-1
        out[x,y] = w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1]
    end
    out
end

function laplacian!(out::Nodes{Primal,NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-2, x in 2:NX-2
        out[x,y] = w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1]
    end
    out
end

function laplacian!(out::Edges{Dual,NX, NY}, w::Edges{Dual,NX, NY}) where {NX, NY}
  @inbounds for y in 2:NY-1, x in 2:NX-2
      out.u[x,y] = w.u[x,y-1] + w.u[x-1,y] - 4w.u[x,y] + w.u[x+1,y] + w.u[x,y+1]
  end
  @inbounds for y in 2:NY-2, x in 2:NX-1
      out.v[x,y] = w.v[x,y-1] + w.v[x-1,y] - 4w.v[x,y] + w.v[x+1,y] + w.v[x,y+1]
  end
  out
end

function laplacian!(out::Edges{Primal,NX, NY}, w::Edges{Primal,NX, NY}) where {NX, NY}
  @inbounds for y in 2:NY-2, x in 2:NX-1
      out.u[x,y] = w.u[x,y-1] + w.u[x-1,y] - 4w.u[x,y] + w.u[x+1,y] + w.u[x,y+1]
  end
  @inbounds for y in 2:NY-1, x in 2:NX-2
      out.v[x,y] = w.v[x,y-1] + w.v[x-1,y] - 4w.v[x,y] + w.v[x+1,y] + w.v[x,y+1]
  end
  out
end

"""
    laplacian(w)

Evaluate the discrete Laplacian of `w`. The data `w` can be of type dual/primary
nodes or edges. The returned result is of the same type as `w`.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[2,2] = 1.0;

julia> laplacian(q)
Whirl.Fields.Edges{Whirl.Fields.Primal,8,6} data
u (in grid orientation):
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0   1.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  -4.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function laplacian(w::Nodes{T,NX,NY}) where {T<:CellType,NX,NY}
  laplacian!(Nodes(T,(NX,NY)), w)
end

function laplacian(w::Edges{T,NX,NY}) where {T<:CellType,NX,NY}
  laplacian!(Edges(T,(NX,NY)), w)
end


struct Laplacian{NX, NY, R}
    conv::Nullable{CircularConvolution{NX, NY}}
end

"""
    Laplacian(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE])

Constructor to set up an operator for evaluating the discrete Laplacian on
dual or primal nodal data of dimension `dims`. If the optional keyword
`with_inverse` is set to `true`, then it also sets up the inverse Laplacian
(the lattice Green's function). These can then be applied, respectively, with
`*` and `\\` operations on data of the appropriate size.

# Example

```jldoctest
julia> w = Nodes(Dual,(7,7));

julia> w[4,4] = 1.0;

julia> L = Laplacian(7,7;with_inverse=true)
Discrete Laplacian (and inverse) on a (nx = 7, ny = 7) grid

julia> s = L\\w
Whirl.Fields.Nodes{Whirl.Fields.Dual,7,7} data
Printing in grid orientation (lower left is (1,1)):
 0.488075  0.462207  0.440376   0.430281     0.440376  0.462207  0.488075
 0.462207  0.424413  0.38662    0.36338      0.38662   0.424413  0.462207
 0.440376  0.38662   0.31831    0.25         0.31831   0.38662   0.440376
 0.430281  0.36338   0.25      -1.26132e-16  0.25      0.36338   0.430281
 0.440376  0.38662   0.31831    0.25         0.31831   0.38662   0.440376
 0.462207  0.424413  0.38662    0.36338      0.38662   0.424413  0.462207
 0.488075  0.462207  0.440376   0.430281     0.440376  0.462207  0.488075

julia> L*s ≈ w
true
```
"""
function Laplacian(dims::Tuple{Int,Int};
                   with_inverse = false, fftw_flags = FFTW.ESTIMATE)
    NX, NY = dims
    if !with_inverse
        return Laplacian{NX, NY, false}(Nullable())
    end

    G = view(LGF_TABLE, 1:NX, 1:NY)
    Laplacian{NX, NY, true}(Nullable(CircularConvolution(G, fftw_flags)))
end

function Laplacian(nx::Int, ny::Int; with_inverse = false, fftw_flags = FFTW.ESTIMATE)
    Laplacian((nx, ny), with_inverse = with_inverse, fftw_flags = fftw_flags)
end

"""
    Laplacian(w::Nodes,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE])

Construct the Laplacian operator with size appropriate for nodal data `w`.
"""
function Laplacian(nodes::Nodes{T,NX,NY}; with_inverse = false, fftw_flags = FFTW.ESTIMATE) where {T<:CellType,NX,NY}
    Laplacian(node_inds(T,(NX,NY)), with_inverse = with_inverse, fftw_flags = fftw_flags)
end


function Base.show(io::IO, L::Laplacian{NX, NY, R}) where {NX, NY, R}
    nodedims = "(nx = $NX, ny = $NY)"
    inverse = R ? " (and inverse)" : ""
    print(io, "Discrete Laplacian$inverse on a $nodedims grid")
end

A_mul_B!(out::Nodes{T,NX,NY}, L::Laplacian, s::Nodes{T,NX,NY}) where {T<:CellType,NX,NY} = laplacian!(out, s)
L::Laplacian * s::Nodes{T,NX,NY} where {T <: CellType, NX,NY} = laplacian(s)

function A_ldiv_B!(out::Nodes{T,NX, NY},
                   L::Laplacian{MX, MY, true},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY}

    A_mul_B!(out.data, get(L.conv), s.data)
    out
end

function (\)(L::Laplacian,s::Nodes{T,NX,NY}) where {T <: CellType, NX,NY}
  A_ldiv_B!(Nodes(T,s), L, s)
end

# Integrating factor

struct IntFact{NX, NY, a}
    conv::Nullable{CircularConvolution{NX, NY}}
end

"""
    IntFact(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])

Constructor to set up an operator for evaluating the integrating factor with
real-valued parameter `a`. This can then be applied with the `*` operation on
data of the appropriate size.

# Example

```jldoctest
julia> w = Nodes(Dual,(6,6));

julia> w[4,4] = 1.0;

julia> E = IntFact(1.0,(6,6))
Integrating factor with parameter 1.0 on a (nx = 6, ny = 6) grid

julia> E*w
Whirl.Fields.Nodes{Whirl.Fields.Dual,6,6} data
Printing in grid orientation (lower left is (1,1)):
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00888233   0.028765    0.0664124   0.0951774   0.0664124   0.028765
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.000828935  0.00268447  0.00619787  0.00888233  0.00619787  0.00268447
```
"""
function IntFact(a::Real,dims::Tuple{Int,Int};fftw_flags = FFTW.ESTIMATE)
    NX, NY = dims

    if a == 0
      return Identity()
    end

    #qtab = [intfact(x, y, a) for x in 0:NX-1, y in 0:NY-1]
    Nmax = 0
    while intfact(Nmax,0,a) > eps(Float64)
      Nmax += 1
    end
    qtab = [max(x,y) <= Nmax ? intfact(x, y, a) : 0.0 for x in 0:NX-1, y in 0:NY-1]
    IntFact{NX, NY, a}(Nullable(CircularConvolution(qtab, fftw_flags)))
end

"""
    IntFact(a::Real,w::Nodes,[fftw_flags=FFTW.ESTIMATE])

Construct the integrating factor operator with size appropriate for nodal data `w`.
"""
function IntFact(a::Real,nodes::Nodes{T,NX,NY}; fftw_flags = FFTW.ESTIMATE) where {T<:CellType,NX,NY}
    IntFact(a,node_inds(T,(NX,NY)), fftw_flags = fftw_flags)
end

function Base.show(io::IO, E::IntFact{NX, NY, a}) where {NX, NY, a}
    nodedims = "(nx = $NX, ny = $NY)"
    print(io, "Integrating factor with parameter $a on a $nodedims grid")
end

function A_mul_B!(out::Nodes{T,NX, NY},
                   E::IntFact{MX, MY, a},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, a}

    A_mul_B!(out.data, get(E.conv), s.data)
    out
end

function (*)(E::IntFact,s::Nodes{T,NX,NY}) where {T <: CellType, NX,NY}
  A_mul_B!(Nodes(T,s), E, s)
end

# Identity

struct Identity end


(*)(::Identity,s::Union{Nodes,Edges}) = s


"""
    curl!(q::Edges{Primal},w::Nodes{Dual})

Evaluate the discrete curl of `w` and return it as `q`.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Primal,w);

julia> curl!(q,w)
Whirl.Fields.Edges{Whirl.Fields.Primal,8,6} data
u (in grid orientation):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  -1.0  1.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function curl!(edges::Edges{Primal, NX, NY},
               s::Nodes{Dual,NX, NY}) where {NX, NY}

    @inbounds for y in 1:NY-1, x in 1:NX
        edges.u[x,y] = s[x,y+1] - s[x,y]
    end

    @inbounds for y in 1:NY, x in 1:NX-1
        edges.v[x,y] = s[x,y] - s[x+1,y]
    end
    edges
end

"""
    curl(w::Nodes{Dual}) --> Edges{Primal}

Evaluate the discrete curl of `w`. Another way to perform this operation is
to construct a `Curl` object and apply it with `*`.

# Example

```jldoctest
julia> C = Curl();

julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> C*w
Whirl.Fields.Edges{Whirl.Fields.Primal,8,6} data
u (in grid orientation):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  -1.0  1.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
```
"""
curl(nodes::Nodes{Dual,NX,NY}) where {NX,NY} = curl!(Edges(Primal, nodes), nodes)

function curl!(nodes::Nodes{Dual,NX, NY},
               edges::Edges{Primal, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v
    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = u[x,y-1] - u[x,y] - v[x-1,y] + v[x,y]
    end
    nodes
end

function curl(edges::Edges{Primal, NX, NY}) where {NX, NY}
    curl!(Nodes(Dual,(NX, NY)), edges)
end

struct Curl end

(*)(::Curl,w::Union{Nodes{Dual,NX,NY},Edges{Primal,NX,NY}}) where {NX,NY} = curl(w)

"""
    divergence!(w::Nodes,q::Edges)

Evaluate the discrete divergence of edge data `q` and return it as nodal data `w`.
Note that `q` can be either primal or dual edge data, but `w` must be of the same
cell type.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> w = Nodes(Primal,(8,6));

julia> divergence!(w,q)
Whirl.Fields.Nodes{Whirl.Fields.Primal,8,6} data
Printing in grid orientation (lower left is (1,1)):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function divergence!(nodes::Nodes{Primal, NX, NY},
                     edges::Edges{Primal, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 1:NY-1, x in 1:NX-1
        nodes[x,y] = - u[x,y] + u[x+1,y] - v[x,y] + v[x,y+1]
    end
    nodes
end

function divergence!(nodes::Nodes{Dual, NX, NY},
                     edges::Edges{Dual, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = - u[x-1,y] + u[x,y] - v[x,y-1] + v[x,y]
    end
    nodes
end

"""
    divergence(q::Edges) --> Nodes

Evaluate the discrete divergence of edge data `q`. Can also perform this operation
by creating an object of Divergence type and applying it with `*`.

# Example

```jldoctest
julia> D = Divergence();

julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> D*q
Whirl.Fields.Nodes{Whirl.Fields.Primal,8,6} data
Printing in grid orientation (lower left is (1,1)):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function divergence(edges::Edges{T, NX, NY}) where {T <: CellType, NX, NY}
    divergence!(Nodes(T, NX, NY), edges)
end

struct Divergence end

(*)(::Divergence,w::Edges{T,NX,NY}) where {T<:CellType,NX,NY} = divergence(w)


"""
    grad!(q::Edges{Primal},w::Nodes{Primal})

Evaluate the discrete gradient of primal nodal data `w` and return it as primal
edge data `q`.

# Example

```jldoctest
julia> w = Nodes(Primal,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Primal,(8,6));

julia> grad!(q,w)
Whirl.Fields.Edges{Whirl.Fields.Primal,8,6} data
u (in grid orientation):
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function grad!(edges::Edges{Primal, NX, NY},
                     p::Nodes{Primal, NX, NY}) where {NX, NY}

    @inbounds for y in 1:NY-1, x in 2:NX-1
        edges.u[x,y] = - p[x-1,y] + p[x,y]
    end
    @inbounds for y in 2:NY-1, x in 1:NX-1
        edges.v[x,y] = - p[x,y-1] + p[x,y]
    end
    edges
end

"""
    grad(w::Nodes{Primal}) --> Edges{Primal}

Evaluate the discrete gradient of primal nodal data `w`. Can also perform this
operation by creating an object of Grad type and applying it with `*`.

# Example

```jldoctest
julia> w = Nodes(Primal,(8,6));

julia> w[3,4] = 1.0;

julia> G = Grad();

julia> G*w
Whirl.Fields.Edges{Whirl.Fields.Primal,8,6} data
u (in grid orientation):
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function grad(p::Nodes{Primal, NX, NY}) where {NX, NY}
  grad!(Edges(Primal,(NX,NY)),p)
end


function grad!(d::EdgeGradient{Primal, Dual, NX, NY},
                     edges::Edges{Primal, NX, NY}) where {NX, NY}

    @inbounds for y in 1:NY-1, x in 1:NX-1
        d.dudx[x,y] = - edges.u[x,y] + edges.u[x+1,y]
        d.dvdy[x,y] = - edges.v[x,y] + edges.v[x,y+1]
    end
    @inbounds for y in 2:NY-1, x in 2:NX-1
        d.dudy[x,y] = - edges.u[x,y-1] + edges.u[x,y]
        d.dvdx[x,y] = - edges.v[x-1,y] + edges.v[x,y]
    end
    d
end

function grad!(d::EdgeGradient{Dual, Primal, NX, NY},
                     edges::Edges{Dual, NX, NY}) where {NX, NY}

    @inbounds for y in 2:NY-1, x in 2:NX-1
        d.dudx[x,y] = - edges.u[x-1,y] + edges.u[x,y]
        d.dvdy[x,y] = - edges.v[x,y-1] + edges.v[x,y]
    end
    @inbounds for y in 1:NY-1, x in 1:NX-1
        d.dudy[x,y] = - edges.u[x,y] + edges.u[x,y+1]
        d.dvdx[x,y] = - edges.v[x,y] + edges.v[x+1,y]
    end
    d
end

function grad(edges::Edges{C, NX, NY}) where {C<:CellType,NX,NY}
  grad!(EdgeGradient(C,(NX,NY)),edges)
end

struct Grad end

(*)(::Grad,w::Union{Nodes{Primal, NX, NY},Edges{Dual,NX,NY},Edges{Primal,NX,NY}}) where {NX,NY} = grad(w)


struct Regularize{N,DV}

  "x values of points, normalized to grid index space"
  x :: Vector{Float64}

  "y values of points, normalized to grid index space"
  y :: Vector{Float64}

  "Discrete Delta function"
  ddf :: DDF
end

"""
    Regularize(x,y,dx,[ddftype=Roma],[I0=(1,1)])

Constructor to set up an operator for regularizing and interpolating data from/to
points immersed in the grid to/from fields on the grid itself. The supplied
`x` and `y` represent physical coordinates of the immersed points, and `dx`
denotes a uniform physical cell size of the grid. The separate arguments `x` and
`y` can be replaced by a single argument `X` of type `VectorData` holding the
coordinates.

The operations of regularization and interpolation are carried out with a discrete
delta function (ddf), which defaults to the type `Roma`. Others are also possible,
such as `Goza`. The optional tuple
`I0` represents the indices of the primary node that coincides with `(x,y) = (0,0)`.
This defaults to `(1,1)`, which leaves one layer of ghost (dual) cells and sets
the physical origin in the lower left corner of the grid of interior dual cells.

The resulting operator can be used in either direction, with the first argument
representing the *target* (the entity to regularize/interpolate to), and the second argument
the *source* (the entity to regularize/interpolate from).

# Example

In the example below, we set up a 12 x 12 grid. Using the default value for `I0`
and setting `dx = 0.1`, the physical dimensions of the non-ghost part of the grid
are 1.0 x 1.0. Three points are set up in the interior, and a vector field is assigned
to them, with the x component of each of them set to 1.0. These data are regularized
to a field of primal edges on the grid.

```jldoctest
julia> x = [0.25,0.75,0.25]; y = [0.75,0.25,0.25];

julia> X = VectorData(x,y);

julia> q = Edges(Primal,(12,12));

julia> dx = 0.1;

julia> H = Regularize(x,y,dx)
Regularization operator with 3 points

julia> f = VectorData(X);

julia> fill!(f.u,1.0);

julia> H(q,f)
Whirl.Fields.Edges{Whirl.Fields.Primal,12,12} data
u (in grid orientation):
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0
 0.0  0.0  8.33333  33.3333  8.33333     8.33333  33.3333  8.33333  0.0  0.0
 0.0  0.0  0.0       0.0     0.0         0.0       0.0     0.0      0.0  0.0
 0.0  0.0  0.0       0.0     0.0      …  0.0       0.0     0.0      0.0  0.0
v (in grid orientation):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function Regularize(x::Vector{T},y::Vector{T},dx::T;
                    ddftype=Roma,
                    I0::Tuple{Int,Int}=(1,1)) where {T<:Real}
  @assert length(x)==length(y)
  Regularize{length(x),dx*dx}(x/dx+I0[1],y/dx+I0[2],DDF(ddftype=ddftype,dx=1.0))
end

Regularize(x::T,y::T,a...;b...) where {T<:Real} = Regularize([x],[y],a...;b...)

Regularize(x::VectorData,a...;b...) = Regularize(x.u,x.v,a...;b...)

function Base.show(io::IO, H::Regularize{N}) where {N}
    print(io, "Regularization operator with $N points")
end

# These regularization operations should be easy to macro-generate

function (H::Regularize{N,DV})(out::Edges{Primal,NX,NY},f::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX
    out.u[x,y] = dot(H.ddf.(x-0.5-H.x,y-H.y),f.u)/DV
    out.v[x,y] = dot(H.ddf.(x-H.x,y-0.5-H.y),f.v)/DV
  end
  out
end

function (H::Regularize{N,DV})(out::Edges{Dual,NX,NY},f::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX-1
    out.u[x,y] = dot(H.ddf.(x-H.x,y-0.5-H.y),f.u)/DV
    out.v[x,y] = dot(H.ddf.(x-0.5-H.x,y-H.y),f.v)/DV
  end
  out
end

function (H::Regularize{N,DV})(out::Tuple{Nodes{Dual,NX,NY},Nodes{Primal,NX,NY}},f::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    out[1][x,y] = dot(H.ddf.(x-0.5-H.x,y-0.5-H.y),f.u)/DV
  end
  @inbounds for y in 1:NY-1, x in 1:NX-1
    out[2][x,y] = dot(H.ddf.(x-H.x,y-H.y),f.v)/DV
  end
  out
end

function (H::Regularize{N,DV})(out::Tuple{Nodes{Primal,NX,NY},Nodes{Dual,NX,NY}},f::VectorData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    out[1][x,y] = dot(H.ddf.(x-H.x,y-H.y),f.u)/DV
  end
  @inbounds for y in 1:NY, x in 1:NX
    out[2][x,y] = dot(H.ddf.(x-0.5-H.x,y-0.5-H.y),f.v)/DV
  end
  out
end

function (H::Regularize{N,DV})(out::Nodes{Primal,NX,NY},f::ScalarData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    out[x,y] = dot(H.ddf.(x-H.x,y-H.y),f.data)/DV
  end
  out
end

function (H::Regularize{N,DV})(out::Nodes{Dual,NX,NY},f::ScalarData{N}) where {N,DV,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    out[x,y] = dot(H.ddf.(x-0.5-H.x,y-0.5-H.y),f.data)/DV
  end
  out
end

function (H::Regularize)(out,f::VectorData,ds::Union{Float64,Vector{Float64}})
  fds = deepcopy(f)
  fds.u .*= ds
  fds.v .*= ds
  H(out,fds)
end

function (H::Regularize)(out,f::ScalarData,ds::Union{Float64,Vector{Float64}})
  fds = deepcopy(f)
  fds .*= ds
  H(out,fds)
end


function (H::Regularize{N})(f::VectorData{N},out::Edges{Primal,NX,NY}) where {N,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX
    f.u .+= H.ddf.(x-0.5-H.x,y-H.y)*out.u[x,y]
    f.v .+= H.ddf.(x-H.x,y-0.5-H.y)*out.v[x,y]
  end
  f
end

function (H::Regularize{N})(f::VectorData{N},out::Edges{Dual,NX,NY}) where {N,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX-1
    f.u .+= H.ddf.(x-H.x,y-0.5-H.y)*out.u[x,y]
    f.v .+= H.ddf.(x-0.5-H.x,y-H.y)*out.v[x,y]
  end
  f
end

function (H::Regularize{N})(f::VectorData{N},out::Tuple{Nodes{Dual,NX,NY},Nodes{Primal,NX,NY}}) where {N,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    f.u .+= H.ddf.(x-0.5-H.x,y-0.5-H.y)*out[1][x,y]
  end
  @inbounds for y in 1:NY-1, x in 1:NX-1
    f.v .+= H.ddf.(x-H.x,y-H.y)*out[2][x,y]
  end
  f
end

function (H::Regularize{N})(f::VectorData{N},out::Tuple{Nodes{Primal,NX,NY},Nodes{Dual,NX,NY}}) where {N,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    f.u .+= H.ddf.(x-H.x,y-H.y)*out[1][x,y]
  end
  @inbounds for y in 1:NY, x in 1:NX
    f.v .+= H.ddf.(x-0.5-H.x,y-0.5-H.y)*out[2][x,y]
  end
  f
end


function (H::Regularize{N})(f::ScalarData{N},out::Nodes{Primal,NX,NY}) where {N,NX,NY}
  @inbounds for y in 1:NY-1, x in 1:NX-1
    f.data .+= H.ddf.(x-H.x,y-H.y)*out[x,y]
  end
  f
end

function (H::Regularize{N})(f::ScalarData{N},out::Nodes{Dual,NX,NY}) where {N,NX,NY}
  @inbounds for y in 1:NY, x in 1:NX
    f.data .+= H.ddf.(x-0.5-H.x,y-0.5-H.y)*out[x,y]
  end
  f
end
