include("lgf.jl")
include("convolution.jl")

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

function Base.show(io::IO, L::Laplacian{NX, NY, R}) where {NX, NY, R}
    nodedims = "(nx = $NX, ny = $NY)"
    inverse = R ? " (and inverse)" : ""
    print(io, "Discrete Laplacian$inverse on a $nodedims grid")
end

A_mul_B!(out::Nodes{Dual,NX,NY}, L::Laplacian, s::Nodes{Dual,NX,NY}) where {NX,NY} = laplacian!(out, s)
L::Laplacian * s::Nodes{Dual,NX,NY} where {NX,NY} = laplacian(s)

function A_ldiv_B!(out::Nodes{Dual,NX, NY},
                   L::Laplacian{NX, NY, true},
                   s::Nodes{Dual, NX, NY}) where {NX, NY}

    A_mul_B!(out.data, get(L.conv), s.data)
    out
end
L::Laplacian \ s::Nodes{Dual,NX,NY} where {NX,NY} = A_ldiv_B!(Nodes(Dual,size(s)), L, s)

# curl
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

# divergence

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

function divergence(edges::Edges{T, NX, NY}) where {T <: CellType, NX, NY}
    divergence!(Nodes(T, NX, NY), edges)
end

# grad
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
