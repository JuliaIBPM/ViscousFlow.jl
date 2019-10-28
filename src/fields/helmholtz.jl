# Helmholtz operators

"""
    helmholtz!(v,w,α)

Evaluate the discrete Helmholtz operator (iα - L) on `w` and return it as `v`.
The data `w` can be of type dual/primary nodes or edges; `v` must be of the same type.
However, both have to be of complex data type.
"""
function helmholtz!(out::Nodes{Dual,NX, NY,T}, w::Nodes{Dual,NX, NY,T}, α::Number) where {NX, NY, T<:ComplexF64}
    @inbounds for y in 2:NY-1, x in 2:NX-1
        out[x,y] = im*α*w[x,y] - (w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1])
    end
    out
end

function helmholtz!(out::Nodes{Primal,NX, NY,T}, w::Nodes{Primal,NX, NY,T}, α::Number) where {NX, NY, T<:ComplexF64}
    @inbounds for y in 2:NY-2, x in 2:NX-2
        out[x,y] = im*α*w[x,y] - (w[x,y-1] + w[x-1,y] - 4w[x,y] + w[x+1,y] + w[x,y+1])
    end
    out
end

function helmholtz!(out::Edges{Dual,NX, NY,T}, w::Edges{Dual,NX, NY,T}, α::Number) where {NX, NY, T<:ComplexF64}
  @inbounds for y in 2:NY-1, x in 2:NX-2
      out.u[x,y] = im*α*w.u[x,y] - (w.u[x,y-1] + w.u[x-1,y] - 4w.u[x,y] + w.u[x+1,y] + w.u[x,y+1])
  end
  @inbounds for y in 2:NY-2, x in 2:NX-1
      out.v[x,y] = im*α*w.v[x,y] - (w.v[x,y-1] + w.v[x-1,y] - 4w.v[x,y] + w.v[x+1,y] + w.v[x,y+1])
  end
  out
end

function helmholtz!(out::Edges{Primal,NX, NY, T}, w::Edges{Primal,NX, NY, T}, α::Number) where {NX, NY, T<:ComplexF64}
  @inbounds for y in 2:NY-2, x in 2:NX-1
      out.u[x,y] = im*α*w.u[x,y] - (w.u[x,y-1] + w.u[x-1,y] - 4w.u[x,y] + w.u[x+1,y] + w.u[x,y+1])
  end
  @inbounds for y in 2:NY-1, x in 2:NX-2
      out.v[x,y] = im*α*w.v[x,y] - (w.v[x,y-1] + w.v[x-1,y] - 4w.v[x,y] + w.v[x+1,y] + w.v[x,y+1])
  end
  out
end

"""
    helmholtz(w,α)

Evaluate the discrete Helmholtz operator (iα - L) on `w`. The data `w` can be of complex type dual/primary
nodes or edges. The returned result is of the same type as `w`.

# Example
"""
function helmholtz(w::Nodes{C,NX,NY,T},α::Number) where {C<:CellType,NX,NY,T<:ComplexF64}
  helmholtz!(Nodes(C,w), w, α)
end

function helmholtz(w::Edges{C,NX,NY,T},α::Number) where {C<:CellType,NX,NY,T<:ComplexF64}
  helmholtz!(Edges(C,w), w, α)
end

"""
    plan_helmholtz(dims::Tuple,α::Number,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],
                          [dx=1.0])

Constructor to set up an operator for evaluating the discrete Helmholtz operator on
complex dual or primal nodal data of dimension `dims`. If the optional keyword
`with_inverse` is set to `true`, then it also sets up the inverse Helmholtz operator
(the lattice Green's function, LGF). These can then be applied, respectively, with
`*` and `\\` operations on data of the appropriate size. The optional parameter
`dx` is used in adjusting the uniform value of the LGF to match the behavior
of the continuous analog at large distances; this is set to 1.0 by default.

Instead of the first argument, one can also supply `w::Nodes` to specify the
size of the domain.

# Example

"""
function plan_helmholtz end

"""
    plan_helmholtz!(dims::Tuple,α::Number,[with_inverse=false],[fftw_flags=FFTW.ESTIMATE],
                          [dx=1.0])

Same as [`plan_helmholtz`](@ref), but forms an operator that works in-place on the
data it operates on.
"""
function plan_helmholtz! end

struct Helmholtz{NX, NY, R, DX, inplace}
    α::Number
    conv::Union{CircularConvolution{NX, NY},Nothing}
end



for (lf,inplace) in ((:plan_helmholtz,false),
                     (:plan_helmholtz!,true))
    @eval function $lf(dims::Tuple{Int,Int},α::Number;
                   with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0)
        NX, NY = dims
        if !with_inverse
            return Helmholtz{NX, NY, false, dx, $inplace}(nothing)
        end
        lgfh_table = load_lgf_helmholtz(NX+1,α)
        G = view(lgfh_table, 1:NX, 1:NY)
        Helmholtz{NX, NY, true, dx, $inplace}(α,CircularConvolution(G, fftw_flags,dtype=ComplexF64))
    end

    @eval function $lf(nx::Int, ny::Int,α::Number;
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0)
        $lf((nx, ny), α, with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx)
    end

    @eval function $lf(nodes::Nodes{C,NX,NY,T},α::Number;
        with_inverse = false, fftw_flags = FFTW.ESTIMATE, dx = 1.0) where {C<:CellType,NX,NY,T<:ComplexF64}
        $lf(node_inds(C,(NX,NY)), α, with_inverse = with_inverse, fftw_flags = fftw_flags, dx = dx)
    end
end



function Base.show(io::IO, L::Helmholtz{NX, NY, R, DX, inplace}) where {NX, NY, R, DX, inplace}
    nodedims = "(nx = $NX, ny = $NY)"
    inverse = R ? " (and inverse)" : ""
    isinplace = inplace ? " in-place" : ""
    print(io, "Discrete$isinplace Helmholtz operator$inverse on a $nodedims grid with spacing $DX")
end

mul!(out::Nodes{C,NX,NY,T}, L::Helmholtz, s::Nodes{C,NX,NY,T}) where {C<:CellType,NX,NY,T<:ComplexF64} = helmholtz!(out, s, L.α)
*(L::Helmholtz{MX,MY,R,DX,false}, s::Nodes{C,NX,NY,T}) where {MX,MY,R,DX,C <: CellType,NX,NY,T<:ComplexF64} =
      helmholtz(s,L.α)
function (*)(L::Helmholtz{MX,MY,R,DX,true}, s::Nodes{C,NX,NY,T}) where {MX,MY,R,DX,C <: CellType,NX,NY,T<:ComplexF64}
    helmholtz!(s,deepcopy(s),L.α)
end


function ldiv!(out::Nodes{C,NX, NY,T},
                   L::Helmholtz{MX, MY, true, DX, inplace},
                   s::Nodes{C, NX, NY,T}) where {C <: CellType, NX, NY, T<:ComplexF64, MX, MY, DX, inplace}

    mul!(out.data, L.conv, s.data)

    # Adjust the behavior at large distance to match continuous kernel
    #out.data .-= (sum(s.data)/2π)*(GAMMA+log(8)/2-log(DX))
    out
end

\(L::Helmholtz{MX,MY,R,DX,false},s::Nodes{C,NX,NY,T}) where {MX,MY,R,DX,C <: CellType,NX,NY,T<:ComplexF64} =
  ldiv!(Nodes(C,s), L, s)

\(L::Helmholtz{MX,MY,R,DX,true},s::Nodes{C,NX,NY,T}) where {MX,MY,R,DX,C <: CellType,NX,NY,T<:ComplexF64} =
  ldiv!(s, L, deepcopy(s))
