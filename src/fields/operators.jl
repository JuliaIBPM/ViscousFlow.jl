include("convolution.jl")
include("lgf.jl")
include("intfact.jl")
include("ddf.jl")

import Base: *, \

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
Nodes{Dual,8,6} data
Printing in grid orientation (lower left is (1,1))
6×8 Array{Float64,2}:
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
Edges{Primal,8,6} data
u (in grid orientation)
5×8 Array{Float64,2}:
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0   1.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  -4.0  1.0  0.0  0.0  0.0  0.0  0.0
 0.0   0.0  0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
6×7 Array{Float64,2}:
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


"""
    plan_laplacian(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],
                          [dx=1.0])

Constructor to set up an operator for evaluating the discrete Laplacian on
dual or primal nodal data of dimension `dims`. If the optional keyword
`with_inverse` is set to `true`, then it also sets up the inverse Laplacian
(the lattice Green's function, LGF). These can then be applied, respectively, with
`*` and `\\` operations on data of the appropriate size. The optional parameter
`dx` is used in adjusting the uniform value of the LGF to match the behavior
of the continuous analog at large distances; this is set to 1.0 by default.

Instead of the first argument, one can also supply `w::Nodes` to specify the
size of the domain.

# Example

```jldoctest
julia> w = Nodes(Dual,(5,5));

julia> w[3,3] = 1.0;

julia> L = plan_laplacian(5,5;with_inverse=true)
Discrete Laplacian (and inverse) on a (nx = 5, ny = 5) grid with spacing 1.0

julia> s = L\\w
Nodes{Dual,5,5} data
Printing in grid orientation (lower left is (1,1))
5×5 Array{Float64,2}:
 0.16707    0.129276     0.106037     0.129276    0.16707
 0.129276   0.0609665   -0.00734343   0.0609665   0.129276
 0.106037  -0.00734343  -0.257343    -0.00734343  0.106037
 0.129276   0.0609665   -0.00734343   0.0609665   0.129276
 0.16707    0.129276     0.106037     0.129276    0.16707

julia> L*s ≈ w
true
```
"""
function plan_laplacian end

"""
    plan_laplacian!(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],
                          [dx=1.0])

Same as [`plan_laplacian`](@ref), but operates in-place on data.
"""
function plan_laplacian! end

struct Laplacian{NX, NY, R, DX, inplace}
    conv::Union{CircularConvolution{NX, NY},Nothing}
end



for (lf,inplace) in ((:plan_laplacian,false),
                     (:plan_laplacian!,true))
    @eval function $lf(dims::Tuple{Int,Int};
                   with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0)
        NX, NY = dims
        if !with_inverse
            #return Laplacian{NX, NY, false, dx, $inplace}(Nullable())
            return Laplacian{NX, NY, false, dx, $inplace}(nothing)
        end

        G = view(LGF_TABLE, 1:NX, 1:NY)
        #Laplacian{NX, NY, true, dx, $inplace}(Nullable(CircularConvolution(G, fftw_flags)))
        Laplacian{NX, NY, true, dx, $inplace}(CircularConvolution(G, fftw_flags))
    end

    @eval function $lf(nx::Int, ny::Int;
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0)
        $lf((nx, ny), with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx)
    end

    @eval function $lf(nodes::Nodes{T,NX,NY};
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0) where {T<:CellType,NX,NY}
        $lf(node_inds(T,(NX,NY)), with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx)
    end
end



function Base.show(io::IO, L::Laplacian{NX, NY, R, DX, inplace}) where {NX, NY, R, DX, inplace}
    nodedims = "(nx = $NX, ny = $NY)"
    inverse = R ? " (and inverse)" : ""
    isinplace = inplace ? " in-place" : ""
    print(io, "Discrete$isinplace Laplacian$inverse on a $nodedims grid with spacing $DX")
end

mul!(out::Nodes{T,NX,NY}, L::Laplacian, s::Nodes{T,NX,NY}) where {T<:CellType,NX,NY} = laplacian!(out, s)
*(L::Laplacian{MX,MY,R,DX,false}, s::Nodes{T,NX,NY}) where {MX,MY,R,DX,T <: CellType,NX,NY} =
      laplacian(s)
function (*)(L::Laplacian{MX,MY,R,DX,true}, s::Nodes{T,NX,NY}) where {MX,MY,R,DX,T <: CellType,NX,NY}
    laplacian!(s,deepcopy(s))
end
#L::Laplacian * s::Nodes{T,NX,NY} where {T <: CellType, NX,NY} = laplacian(s)


function ldiv!(out::Nodes{T,NX, NY},
                   L::Laplacian{MX, MY, true, DX, inplace},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, DX, inplace}

    mul!(out.data, L.conv, s.data)

    # Adjust the behavior at large distance to match continuous kernel
    out.data .-= (sum(s.data)/2π)*(GAMMA+log(8)/2-log(DX))
    out
end

\(L::Laplacian{MX,MY,R,DX,false},s::Nodes{T,NX,NY}) where {MX,MY,R,DX,T <: CellType,NX,NY} =
  ldiv!(Nodes(T,s), L, s)

\(L::Laplacian{MX,MY,R,DX,true},s::Nodes{T,NX,NY}) where {MX,MY,R,DX,T <: CellType,NX,NY} =
  ldiv!(s, L, deepcopy(s))

# Integrating factor

"""
    plan_intfact(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])

Constructor to set up an operator for evaluating the integrating factor with
real-valued parameter `a`. This can then be applied with the `*` operation on
data of the appropriate size.

The `dims` argument can be replaced with `w::Nodes` to specify the size of the
domain.

# Example

```jldoctest
julia> w = Nodes(Dual,(6,6));

julia> w[4,4] = 1.0;

julia> E = plan_intfact(1.0,(6,6))
Integrating factor with parameter 1.0 on a (nx = 6, ny = 6) grid

julia> E*w
Nodes{Dual,6,6} data
Printing in grid orientation (lower left is (1,1))
6×6 Array{Float64,2}:
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00888233   0.028765    0.0664124   0.0951774   0.0664124   0.028765
 0.00619787   0.0200715   0.0463409   0.0664124   0.0463409   0.0200715
 0.00268447   0.00869352  0.0200715   0.028765    0.0200715   0.00869352
 0.000828935  0.00268447  0.00619787  0.00888233  0.00619787  0.00268447
```
"""
function plan_intfact end

"""
    plan_intfact!(a::Real,dims::Tuple,[fftw_flags=FFTW.ESTIMATE])

Same as [`plan_intfact`](@ref), but the resulting operator performs an in-place
operation on data.
"""
function plan_intfact! end


struct IntFact{NX, NY, a, inplace}
    conv::Union{CircularConvolution{NX, NY},Nothing}
end

for (lf,inplace) in ((:plan_intfact,false),
                     (:plan_intfact!,true))

    @eval function $lf(a::Real,dims::Tuple{Int,Int};fftw_flags = FFTW.ESTIMATE)
        NX, NY = dims

        if a == 0
          return IntFact{NX, NY, 0.0, $inplace}(nothing)
        end

        #qtab = [intfact(x, y, a) for x in 0:NX-1, y in 0:NY-1]
        Nmax = 0
        while abs(intfact(Nmax,0,a)) > eps(Float64)
          Nmax += 1
        end
        qtab = [max(x,y) <= Nmax ? intfact(x, y, a) : 0.0 for x in 0:NX-1, y in 0:NY-1]
        #IntFact{NX, NY, a, $inplace}(Nullable(CircularConvolution(qtab, fftw_flags)))
        IntFact{NX, NY, a, $inplace}(CircularConvolution(qtab, fftw_flags))
      end

      @eval $lf(a::Real,nodes::Nodes{T,NX,NY}; fftw_flags = FFTW.ESTIMATE) where {T<:CellType,NX,NY} =
          $lf(a,node_inds(T,(NX,NY)), fftw_flags = fftw_flags)


end


function Base.show(io::IO, E::IntFact{NX, NY, a, inplace}) where {NX, NY, a, inplace}
    nodedims = "(nx = $NX, ny = $NY)"
    isinplace = inplace ? "In-place integrating factor" : "Integrating factor"
    print(io, "$isinplace with parameter $a on a $nodedims grid")
end

function mul!(out::Nodes{T,NX, NY},
                   E::IntFact{MX, MY, a, inplace},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, a, inplace}

    mul!(out.data, E.conv, s.data)
    out
end

function mul!(out::Nodes{T,NX, NY},
                   E::IntFact{MX, MY, 0.0, inplace},
                   s::Nodes{T, NX, NY}) where {T <: CellType, NX, NY, MX, MY, inplace}
    out .= deepcopy(s)
end

*(E::IntFact{MX,MY,a,false},s::Nodes{T,NX,NY}) where {MX,MY,a,T <: CellType, NX,NY} =
  mul!(Nodes(T,s), E, s)

*(E::IntFact{MX,MY,a,true},s::Nodes{T,NX,NY}) where {MX,MY,a,T <: CellType, NX,NY} =
    mul!(s, E, deepcopy(s))


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
Edges{Primal,8,6} data
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
Edges{Primal,8,6} data
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
Nodes{Primal,8,6} data
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
Nodes{Primal,8,6} data
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
Edges{Primal,8,6} data
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
Edges{Primal,8,6} data
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
function grad(p::Nodes{Primal, NX, NY}) where {NX, NY}
  grad!(Edges(Primal,(NX,NY)),p)
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
  grad!(Edges(Dual,(NX,NY)),p)
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
  grad!(EdgeGradient(C,(NX,NY)),edges)
end

struct Grad end

(*)(::Grad,w::Union{Nodes{Primal, NX, NY},Nodes{Dual, NX, NY},Edges{Dual,NX,NY},Edges{Primal,NX,NY}}) where {NX,NY} = grad(w)



include("regularization.jl")
