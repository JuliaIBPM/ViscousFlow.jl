include("lgf.jl")
include("convolution.jl")

import Base: *, \, A_mul_B!, A_ldiv_B!

function laplacian!(out::DualNodes{NX, NY}, w::DualNodes{NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 2:NX-1
        out[x,y] = w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1]
    end
    out
end

function laplacian(w::DualNodes{NX, NY}) where {NX, NY}
    laplacian!(DualNodes(NX, NY), w)
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

A_mul_B!(out::DualNodes, L::Laplacian, s::DualNodes) = laplacian!(out, s)
L::Laplacian * s::DualNodes = laplacian(s)

function A_ldiv_B!(out::DualNodes{NX, NY},
                   L::Laplacian{NX, NY, true},
                   s::DualNodes{NX, NY}) where {NX, NY}

    A_mul_B!(out.data, get(L.conv), s.data)
end
L::Laplacian \ s::DualNodes = A_ldiv_B!(DualNodes(size(s)), L, s)

function curl!(edges::Edges{Primal, NX, NY},
               s::DualNodes{NX, NY}) where {NX, NY}

    @inbounds for y in 1:NY-1, x in 1:NX
        edges.u[x,y] = s[x,y+1] - s[x,y]
    end

    @inbounds for y in 1:NY, x in 1:NX-1
        edges.v[x,y] = s[x,y] - s[x+1,y]
    end
    edges
end

curl(nodes::DualNodes) = curl!(Edges(Primal, nodes), nodes)

function curl!(nodes::DualNodes{NX, NY},
               edges::Edges{Primal, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v
    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = u[x,y-1] - u[x,y] - v[x-1,y] + v[x,y]
    end
    nodes
end

function curl(edges::Edges{Primal, NX, NY}) where {NX, NY}
    curl!(DualNodes(NX, NY), edges)
end

function divergence!(nodes::DualNodes{NX, NY},
                     edges::Edges{Dual, NX, NY}) where {NX, NY}

    u, v = edges.u, edges.v

    @inbounds for y in 2:NY-1, x in 2:NX-1
        nodes[x,y] = - u[x-1,y-1] + u[x,y-1] - v[x-1,y-1] + v[x-1,y]
    end
    nodes
end

function divergence(edges::Edges{Dual, NX, NY}) where {NX, NY}
    divergence!(DualNodes(NX, NY), edges)
end
