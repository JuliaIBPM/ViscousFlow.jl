# Laplacian

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
Nodes{Dual,8,6,Float64} data
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
Edges{Primal,8,6,Float64} data
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
function laplacian(w::Nodes{C,NX,NY}) where {C<:CellType,NX,NY}
  laplacian!(Nodes(C,w), w)
end

function laplacian(w::Edges{C,NX,NY}) where {C<:CellType,NX,NY}
  laplacian!(Edges(C,w), w)
end

"""
    laplacian_symm!(v,w)

Evaluate the symmetric 5-point discrete Laplacian of `w` and return it as `v`. The data `w` can be
of type dual nodes only for now. This symmetric Laplacian also evaluates the
partial Laplacians (using only available stencil data) on the ghost nodes.
"""
function laplacian_symm!(out::Nodes{Dual,NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 2:NX-1
        out[x,y] = w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1]
    end
    @inbounds for y in 2:NY-1
        out[1,y]  = w[1,y-1]            - 4w[1,y] + w[2,y] + w[1,y+1]
        out[NX,y] = w[NX,y-1] + w[NX-1,y]- 4w[NX,y]        + w[NX,y+1]
    end
    @inbounds for x in 2:NX-1
        out[x,1]  = w[x-1,1] + w[x+1,1] - 4w[x,1] + w[x,2]
        out[x,NY] = w[x-1,NY]+ w[x+1,NY]- 4w[x,NY] + w[x,NY-1]
    end
    out[1,1] = -4w[1,1] + w[1,2] + w[2,1]
    out[NX,1] = -4w[NX,1] + w[NX-1,1] + w[NX,2]
    out[1,NY] = -4w[1,NY] + w[1,NY-1] + w[2,NY]
    out[NX,NY] = -4w[NX,NY] + w[NX,NY-1] + w[NX-1,NY]
    out
end


"""
    plan_laplacian(dims::Tuple,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],
                          [dx=1.0],[dtype=Float64])

Constructor to set up an operator for evaluating the discrete Laplacian on
dual or primal nodal data of dimension `dims`. If the optional keyword
`with_inverse` is set to `true`, then it also sets up the inverse Laplacian
(the lattice Green's function, LGF). These can then be applied, respectively, with
`*` and `\\` operations on data of the appropriate size. The optional parameter
`dx` is used in adjusting the uniform value of the LGF to match the behavior
of the continuous analog at large distances; this is set to 1.0 by default. The
type of data on which to act is floating point by default, but can also be ComplexF64.
This is specified with the optional parameter `dtype`

Instead of the first argument, one can also supply `w::Nodes` to specify the
size of the domain.

# Example

```jldoctest
julia> w = Nodes(Dual,(5,5));

julia> w[3,3] = 1.0;

julia> L = plan_laplacian(5,5;with_inverse=true)
Discrete Laplacian (and inverse) on a (nx = 5, ny = 5) grid acting on Float64 data with spacing 1.0

julia> s = L\\w
Nodes{Dual,5,5,Float64} data
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

struct Laplacian{NX, NY, T, R, DX, inplace}
    conv::Union{CircularConvolution{NX, NY, T},Nothing}
end



for (lf,inplace) in ((:plan_laplacian,false),
                     (:plan_laplacian!,true))
    @eval function $lf(dims::Tuple{Int,Int};
                   with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0, dtype = Float64)
        NX, NY = dims
        if !with_inverse
            return Laplacian{NX, NY, dtype, false, dx, $inplace}(nothing)
        end

        G = view(LGF_TABLE, 1:NX, 1:NY)
        Laplacian{NX, NY, dtype, true, dx, $inplace}(CircularConvolution(G, fftw_flags,dtype=dtype))
    end

    @eval function $lf(nx::Int, ny::Int;
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0, dtype = Float64)
        $lf((nx, ny), with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx, dtype = dtype)
    end

    @eval function $lf(nodes::Nodes{T,NX,NY};
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0, dtype = Float64) where {T<:CellType,NX,NY}
        $lf(node_inds(T,(NX,NY)), with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx, dtype = dtype)
    end
end



function Base.show(io::IO, L::Laplacian{NX, NY, T, R, DX, inplace}) where {NX, NY, T, R, DX, inplace}
    nodedims = "(nx = $NX, ny = $NY)"
    inverse = R ? " (and inverse)" : ""
    isinplace = inplace ? " in-place" : ""
    print(io, "Discrete$isinplace Laplacian$inverse on a $nodedims grid acting on $T data with spacing $DX")
end

mul!(out::Nodes{C,NX,NY}, L::Laplacian, s::Nodes{C,NX,NY}) where {C<:CellType,NX,NY} = laplacian!(out, s)
*(L::Laplacian{MX,MY,T,R,DX,false}, s::Nodes{C,NX,NY}) where {MX,MY,T,R,DX,C <: CellType,NX,NY} =
      laplacian(s)
function (*)(L::Laplacian{MX,MY,T,R,DX,true}, s::Nodes{C,NX,NY}) where {MX,MY,T,R,DX,C <: CellType,NX,NY}
    laplacian!(s,deepcopy(s))
end


function ldiv!(out::Nodes{C,NX, NY,T},
                   L::Laplacian{MX, MY, T, true, DX, inplace},
                   s::Nodes{C, NX, NY,T}) where {C <: CellType, NX, NY, MX, MY, T, DX, inplace}

    mul!(out.data, L.conv, s.data)

    # Adjust the behavior at large distance to match continuous kernel
    out.data .-= (sum(s.data)/2π)*(GAMMA+log(8)/2-log(DX))
    out
end

\(L::Laplacian{MX,MY,T,R,DX,false},s::Nodes{C,NX,NY}) where {MX,MY,T,R,DX,C <: CellType,NX,NY} =
  ldiv!(Nodes(C,s), L, s)

\(L::Laplacian{MX,MY,T,R,DX,true},s::Nodes{C,NX,NY}) where {MX,MY,T,R,DX,C <: CellType,NX,NY} =
  ldiv!(s, L, deepcopy(s))
