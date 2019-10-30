# Differential calculus mimetic operators

"""
    curl!(q::Edges{Primal},w::Nodes{Dual})

Evaluate the discrete curl of `w` and return it as `q`.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Primal,w);

julia> curl!(q,w)
Edges{Primal,8,6,Float64} data
u (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
6×7 Array{Float64,2}:
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
Edges{Primal,8,6,Float64} data
u (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
6×7 Array{Float64,2}:
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0  -1.0  1.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0
```
"""
curl(nodes::Nodes{Dual,NX,NY,T}) where {NX,NY,T} = curl!(Edges(Primal, nodes), nodes)

function curl!(nodes::Nodes{Dual,NX, NY},
               edges::Edges{Primal, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v
    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = u[x,y-1] - u[x,y] - v[x-1,y] + v[x,y]
    end
    nodes
end

function curl(edges::Edges{Primal, NX, NY}) where {NX, NY}
    curl!(Nodes(Dual,edges), edges)
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
Nodes{Primal,8,6,Float64} data
Printing in grid orientation (lower left is (1,1))
5×7 Array{Float64,2}:
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
Nodes{Primal,8,6,Float64} data
Printing in grid orientation (lower left is (1,1))
5×7 Array{Float64,2}:
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function divergence(edges::Edges{T, NX, NY}) where {T <: CellType, NX, NY}
    divergence!(Nodes(T, edges), edges)
end

"""
    divergence!(w::XEdges/YEdges,q::NodePair)

Evaluate the discrete divergence of node pair data `q` and return it as data `w`.
Note that `q` can be either primal/dual or dual/primal edge data, and `w`
must be, respectively, primal x-edges/dual y-edges or primal y-edges/dual x-edges type.
"""
function divergence!(nodes::Union{XEdges{Primal, NX, NY},YEdges{Dual, NX, NY}},
                     edges::NodePair{Primal, Dual,NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 1:NY-1, x in 2:NX-1
        nodes[x,y] = - u[x-1,y] + u[x,y] - v[x,y] + v[x,y+1]
    end
    nodes
end

function divergence!(nodes::Union{YEdges{Primal, NX, NY},XEdges{Dual, NX, NY}},
                     edges::NodePair{Dual, Primal,NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 2:NY-1, x in 1:NX-1
        nodes[x,y] = - u[x,y] + u[x+1,y] - v[x,y-1] + v[x,y]
    end
    nodes
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
Edges{Primal,8,6,Float64} data
u (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
6×7 Array{Float64,2}:
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
Edges{Primal,8,6,Float64} data
u (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0   0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0  -1.0  0.0  0.0  0.0  0.0
 0.0  0.0   1.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
 0.0  0.0   0.0  0.0  0.0  0.0  0.0
```
"""
function grad(p::Nodes{Primal, NX, NY,T}) where {NX, NY,T}
  grad!(Edges(Primal,p),p)
end

"""
    grad!(q::Edges{Dual},w::Nodes{Dual})

Evaluate the discrete gradient of dual nodal data `w` and return it as dual
edge data `q`.
"""
function grad!(edges::Edges{Dual, NX, NY},
                   p::Nodes{Dual, NX, NY}) where {NX, NY}

    @inbounds for y in 1:NY, x in 1:NX-1
        edges.u[x,y] = - p[x,y] + p[x+1,y]
    end
    @inbounds for y in 1:NY-1, x in 1:NX
        edges.v[x,y] = - p[x,y] + p[x,y+1]
    end
    edges
end

"""
    grad(w::Nodes{Dual}) --> Edges{Dual}

Evaluate the discrete gradient of dual nodal data `w`. Can also perform this
operation by creating an object of Grad type and applying it with `*`.
"""
function grad(p::Nodes{Dual, NX, NY}) where {NX, NY}
  grad!(Edges(Dual,p),p)
end

"""
    grad!(d::EdgeGradient{Primal,Dual},q::Edges{Primal})

Evaluate the discrete gradient of primal edge data `q` and return it as edge
gradient data `d`, where the diagonal entries of the gradient lie on primal
nodes and the off-diagonal entries lie at dual nodes.
"""
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

"""
    grad!(d::EdgeGradient{Dual,Primal},q::Edges{Dual})

Evaluate the discrete gradient of dual edge data `q` and return it as edge
gradient data `d`, where the diagonal entries of the gradient lie on dual
nodes and the off-diagonal entries lie at primal nodes.
"""
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

"""
    grad(q::Edges{Primal/Dual}) --> EdgeGradient{Dual/Primal,Primal/Dual}

Evaluate the discrete gradient of primal or dual edge data `q`. Can also perform this
operation by creating an object of Grad type and applying it with `*`.
"""
function grad(edges::Edges{C, NX, NY}) where {C<:CellType,NX,NY}
  grad!(EdgeGradient(C,edges),edges)
end

struct Grad end

(*)(::Grad,w::Union{Nodes{Primal, NX, NY},Nodes{Dual, NX, NY},Edges{Dual,NX,NY},Edges{Primal,NX,NY}}) where {NX,NY} = grad(w)



include("differencing1d.jl")
