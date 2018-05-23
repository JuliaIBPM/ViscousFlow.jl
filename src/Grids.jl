"""
The `Grids` module provides a data structure:

[`Grids.DualPatch`](@ref)

and functions that act on this structure.

"""
module Grids

export DualPatch

import Whirl
import Whirl:@get, MappedVector
# import Whirl.NavierStokes
using FastGaussQuadrature

abstract type Grid end

const nodes, weights = gausslegendre(100)

mutable struct DualPatch <: Grid
    "array of number of dual cells in each direction in patch interior"
    N::Vector{Int}

    "uniform grid spacing (in all directions)"
    Δx::Float64

    "coordinates of lower left-hand corner"
    xmin::Vector{Float64}
    xmax::Vector{Float64}

    "first interior cell indices in each direction"
    ifirst::Vector{Int}

    "range of interior indices"
    cellint::Array{UnitRange{Int},1}
    nodeint::Array{UnitRange{Int},1}
    facexint::Array{UnitRange{Int},1}
    faceyint::Array{UnitRange{Int},1}

#=
    "mapping from grid index to matrix row/column"
    cellmap
    nodemap
    facexmap
    faceymap
=#

    "cell center variable"
    cell::Array{Float64,Whirl.ndim}

    "face variable components"
    facex::Array{Float64,Whirl.ndim}
    facey::Array{Float64,Whirl.ndim}

    "dual velocity field components"
    dualfacex::Array{Float64,Whirl.ndim}
    dualfacey::Array{Float64,Whirl.ndim}

    "cell node variable"
    node::Array{Float64,Whirl.ndim}

    "LGF table"
    lgftab::Array{Float64,Whirl.ndim}
    lgfhat::Array{Complex{Float64},Whirl.ndim}

    "integrating factor table"
    qtab::Array{Float64,Whirl.ndim}
    qhat::Array{Complex{Float64},Whirl.ndim}

    "LGF*Q"
    gqhat::Array{Complex{Float64},Whirl.ndim}
    α::Float64

    "preplanned fft"
    fftop::FFTW.rFFTWPlan

end


function DualPatch(N,Δx,xmin)

    xmax = xmin + N*Δx

    # set up grid arrays with ghosts
    cell = zeros(Float64,N[1]+2,N[2]+2)
    node = zeros(Float64,N[1]+1,N[2]+1)
    facex = zeros(Float64,N[1]+2,N[2]+1)
    facey = zeros(Float64,N[1]+1,N[2]+2)
    dualfacex = zeros(Float64,N[1]+1,N[2]+2)
    dualfacey = zeros(Float64,N[1]+2,N[2]+1)

    # set the first interior indices in each direction
    # (this implicitly sets the number of ghost cell layers)
    ifirst = [2,2]

    # set the index ranges of interior data for each type of grid data array
    cellint =  ifirst-1+[1:N[1],1:N[2]];
    nodeint =  ifirst-1+[1:N[1]-1,1:N[2]-1];
    facexint = ifirst-1+[1:N[1],1:N[2]-1];
    faceyint = ifirst-1+[1:N[1]-1,1:N[2]];

#=
    """
        construct maps of grid interior indices to matrix operator
        row/column indices.

	e.g., cellmap accepts entries from cellint[1].start to
	cellint[1].stop in the first argument, and cellint[2].start to
	cellint[2].stop in the second argument, and returns a single
	integer providing the row index for a data vector that would
	hold cell data, or the column index for a matrix that would
	act upon such a data vector

    """
    cellmap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*N[1]
    nodemap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*(N[1]-1)
    facexmap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*N[1]
    faceymap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*(N[1]-1)
    =#

    lgftab = Array{Float64}(0,0)
    lgfhat = Array{Complex{Float64}}(0,0)
    qtab = Array{Float64}(0,0)
    qhat = Array{Complex{Float64}}(0,0)
    gqhat = Array{Complex{Float64}}(0,0)
    α = 0.0

    fftop = setFFTPlan(N)

    DualPatch(N,Δx,xmin,xmax,ifirst,cellint,nodeint,facexint,faceyint,
    cell,facex,facey,dualfacex,dualfacey,node,
    lgftab,lgfhat,qtab,qhat,gqhat,α,fftop)
	      #cellmap,nodemap,facexmap,faceymap,


end

setFFTPlan(N) = FFTW.plan_rfft(zeros(2*N[1]+3,2*N[2]+3))

function setFFTPlan!(g::DualPatch)
  g.fftop = setFFTPlan(g.N)
  nothing
end

######## Grid data types




########

xcell(g::DualPatch) = [g.xmin[1]+g.Δx*(i-g.ifirst[1]+1/2) for i=g.cellint[1]]
ycell(g::DualPatch) = [g.xmin[2]+g.Δx*(j-g.ifirst[2]+1/2) for j=g.cellint[2]]
xnode(g::DualPatch) = [g.xmin[1]+g.Δx*(i-g.ifirst[1]+1) for i=g.nodeint[1]]
ynode(g::DualPatch) = [g.xmin[2]+g.Δx*(j-g.ifirst[2]+1) for j=g.nodeint[2]]
xfacex(g::DualPatch) = [g.xmin[1]+g.Δx*(i-g.ifirst[1]+1/2) for i=g.facexint[1]]
yfacex(g::DualPatch) = [g.xmin[2]+g.Δx*(j-g.ifirst[2]+1) for j=g.facexint[2]]
xfacey(g::DualPatch) = [g.xmin[1]+g.Δx*(i-g.ifirst[1]+1) for i=g.faceyint[1]]
yfacey(g::DualPatch) = [g.xmin[2]+g.Δx*(j-g.ifirst[2]+1/2) for j=g.faceyint[2]]

struct GridOperators{TL}

  L⁻¹ :: TL
  curl :: Function

end


# Differential operations
function curl!(cell,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    @. cell[ir,jr] = -facex[ir,jr]+facex[ir,jr-1]+facey[ir,jr]-facey[ir-1,jr]
    nothing
end

function curl!(facex,facey,ir::UnitRange{Int},jr::UnitRange{Int},cell)
    irx = ir.start:ir.stop+1
    jry = jr.start:jr.stop+1
    @. facex[irx,jr] = cell[irx,jr+1]-cell[irx,jr]
    @. facey[ir,jry] = cell[ir,jry]-cell[ir+1,jry]
    nothing
end

function diverg!(node,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    @. node[ir,jr] = facex[ir+1,jr]-facex[ir,jr]+facey[ir,jr+1]-facey[ir,jr]
    nothing
end

function dualdiverg!(cell,ir::UnitRange{Int},jr::UnitRange{Int},dualfacex,dualfacey)
    @. cell[ir,jr] = dualfacex[ir,jr]-dualfacex[ir-1,jr]+dualfacey[ir,jr]-dualfacey[ir,jr-1]
    nothing
end

function grad!(facex,facey,ir::UnitRange{Int},jr::UnitRange{Int},node)
    iry = ir.start-1:ir.stop
    jrx = jr.start-1:jr.stop
    @. facex[ir,jrx] = node[ir,jrx]-node[ir-1,jrx]
    @. facey[iry,jr] = node[iry,jr]-node[iry,jr-1]
    nothing
end

function grad!(gradface,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    irnode = ir.start-1:ir.stop
    jrnode = jr.start-1:jr.stop
    @. gradface[1,1][irnode,jrnode] =
                           facex[irnode+1,jrnode]-facex[irnode,jrnode] # node
    @. gradface[1,2][ir,jr] = facex[ir,jr]-facex[ir,jr-1] # cell
    @. gradface[2,1][ir,jr] = facey[ir,jr]-facey[ir-1,jr] # cell
    @. gradface[2,2][irnode,jrnode] =
                           facey[irnode,jrnode+1]-facey[irnode,jrnode] # node
    nothing
end

function lap!(lapf,ir::UnitRange{Int},jr::UnitRange{Int},f)
    @. lapf[ir,jr] = f[ir+1,jr]+f[ir-1,jr]+f[ir,jr+1]+f[ir,jr-1]-4f[ir,jr]
    nothing
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    irx = ir.start-1:ir.stop
    jry = jr.start-1:jr.stop
    @. vx[irx,jr] = 0.25(facex[irx,jr]+facex[irx+1,jr]+facex[irx,jr-1]+facex[irx+1,jr-1])
    @. vy[ir,jry] = 0.25(facey[ir-1,jry]+facey[ir-1,jry+1]+facey[ir,jry]+facey[ir,jry+1])
    nothing
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},v)
    irx = ir.start-1:ir.stop
    jry = jr.start-1:jr.stop
    @. vx[irx,jr] = 0.5(v[irx,jr]+v[irx+1,jr])
    @. vy[ir,jry] = 0.5(v[ir,jry]+v[ir,jry+1])
    nothing
end

# Differential operations with grid interface
function curl(g::DualPatch,facex,facey)
    cell = zeros(g.cell)
    curl!(cell,g.cellint[1],g.cellint[2],facex,facey)
    cell
end

curl(g::DualPatch,q::Tuple{Array{Float64,2},Array{Float64,2}}) = curl(g,q[1],q[2])

function curl(g::DualPatch,cell)
    facex = zeros(g.facex)
    facey = zeros(g.facey)
    # also calculate curl in ghost cells
    curl!(facex,facey,g.cellint[1].start-1:g.cellint[1].stop,
	        g.cellint[2].start-1:g.cellint[2].stop,cell)
    facex, facey
end

function diverg(g::DualPatch,facex,facey)
    node = zeros(g.node)
    diverg!(node,g.cellint[1].start-1:g.cellint[1].stop,
	      g.cellint[2].start-1:g.cellint[2].stop,facex,facey)
    node
end

diverg(g::DualPatch,q::Tuple{Array{Float64,2},Array{Float64,2}}) = diverg(g,q[1],q[2])

function dualdiverg(g::DualPatch,dualfacex,dualfacey)
    cell = zeros(g.cell)
    dualdiverg!(cell,g.cellint[1],g.cellint[2],dualfacex,dualfacey)
    cell
end

dualdiverg(g::DualPatch,q::Tuple{Array{Float64,2},Array{Float64,2}}) = dualdiverg(g,q[1],q[2])



function grad(g::DualPatch,node)
    facex = zeros(g.facex)
    facey = zeros(g.facey)
    grad!(facex,facey,g.cellint[1],g.cellint[2],node)
    facex,facey
end

function grad(g::DualPatch,facex,facey)
    gradface = Array{Array{Float64,2}}(2,2)
    gradface[1,1] = gradface[2,2] = zeros(g.node)
    gradface[1,2] = gradface[2,1] = zeros(g.cell)
    grad!(gradface,g.cellint[1],g.cellint[2],facex,facey)
    gradface
end

function lap(g::DualPatch,cell)
    lapcell = zeros(g.cell)
    lap!(lapcell,g.cellint[1],g.cellint[2],cell)
    lapcell
end

function lap(g::DualPatch,facex,facey)
    lapfacex = zeros(g.facex)
    lapfacey = zeros(g.facey)
    lap!(lapfacex,g.facexint[1],g.facexint[2],facex)
    lap!(lapfacey,g.faceyint[1],g.faceyint[2],facey)
    lapfacex,lapfacey
end

lap(g::DualPatch,q::Tuple{Array{Float64,2},Array{Float64,2}}) = lap(g,q[1],q[2])

function shift(g::DualPatch,facex,facey)
    vx = zeros(g.dualfacex)
    vy = zeros(g.dualfacey)
    shift!(vx,vy,g.cellint[1],g.cellint[2],facex,facey)
    vx, vy
end

shift(g::DualPatch,q::Tuple{Array{Float64,2},Array{Float64,2}}) = shift(g,q[1],q[2])

function shift(g::DualPatch,cell)
    cellx = zeros(g.dualfacex)
    celly = zeros(g.dualfacey)
    shift!(cellx,celly,g.cellint[1],g.cellint[2],cell)
    cellx, celly
end

function cross(g::DualPatch,cell,facex,facey)
    vx,vy = shift(g,facex,facey)
    ux,uy = shift(g,cell)
    -uy.*vy, ux.*vx
end

function dualshiftx!(vx,ir::UnitRange{Int},jr::UnitRange{Int},facex)
     irx = ir.start:ir.stop
     jry = jr.start:jr.stop
     @. vx[irx,jry] = 0.5(facex[irx,jry-1]+facex[irx,jry])
     nothing
 end
 function dualshifty!(vy,ir::UnitRange{Int},jr::UnitRange{Int},facey)
     irx = ir.start:ir.stop
     jry = jr.start:jr.stop
     @. vy[irx,jry] = 0.5(facey[irx-1,jry]+facey[irx,jry])
     nothing
 end

 function dualshiftx(g::DualPatch,facex)
     cellx = zeros(g.cell)
     dualshiftx!(cellx,g.cellint[1],g.cellint[2],facex)
     cellx
 end

 function dualshifty(g::DualPatch,facey)
     celly = zeros(g.cell)
     dualshifty!(celly,g.cellint[1],g.cellint[2],facey)
     celly
 end
function addgh!(u)
utp=zeros(size(u,1)+2,size(u,2)+2)
utp[2:1+size(u,1),2:1+size(u,2)]=u
u=utp
end

function rmgh!(u)
utp=zeros(size(u,1)-2,size(u,2)-2)
utp=u[2:size(u,1)-1,2:size(u,2)-1]
u=utp
end
# Lattice Green's function
function lgf(n)

    if n[1]==n[2]==0
       return 0.0
    elseif n[1]>=n[2]
       v = quadgauss() do x
         if x == -1
	          return sqrt(2)abs(n[1])
         else
            t = (x+1)/2
	          return 0.5real((1 -
              ( (t-sqrt(1im))./(t+sqrt(1im)) ).^(n[2]+abs(n[1])).*
     	        ( (t+sqrt(-1im))./(t-sqrt(-1im)) ).^(n[2]-abs(n[1])) ))./t
	        end

        end
        return 0.5v/pi
    else
       return lgf([n[2],n[1]])
    end

end

# Integrating factor
intfact(n,a) = exp(-4a)besseli(n[1],2a)besseli(n[2],2a)

"""
    G = lgf(g::Grid)

Set up a table of lgf values in the upper right quadrant, centered at
the lower left ghost cell in grid 'g'.
"""
lgf(g::Grid) = reshape([lgf([i,j]) for i=0:g.N[1]+1,j=0:g.N[2]+1],
       g.N[1]+2,g.N[2]+2)

"""
    GE = intfact(g::Grid,a::Float64)

Set up a table of integrating factor values in the upper right
quadrant, centered at the lower left ghost cell in grid 'g'.
This zeros out the values lower than machine epsilon.
"""
function intfact(g::Grid,a::Float64)

  Nmax = 0
  while intfact([Nmax,0],a) > eps(Float64)
    Nmax += 1
  end

  return reshape([max(i,j) <= Nmax ?
        intfact([i,j],a) : 0.0
        for i=0:g.N[1]+1,j=0:g.N[2]+1], g.N[1]+2,g.N[2]+2)
end

# intfact(g::Grid,a::Float64) = reshape([intfact([i,j],a) > eps(Float64) ?
#         intfact([i,j],a) : 0.0
#         for i=0:g.N[1]+1,j=0:g.N[2]+1], g.N[1]+2,g.N[2]+2)



"""
    s = convolve(G,cell)

Perform a discrete convolution of grid data `cell` with one of the Green's
function tables (the LGF or the integrating factor). This exploits the
symmetries in these functions.
"""
function convolve(G::Array{T,2},cell::AbstractArray{T,2}) where T
    N = size(cell)

    Gw = zeros(cell)
    for isrc=1:N[1],jsrc=1:N[2]
      if abs(cell[isrc,jsrc])<eps(Float64)
	         continue
	    end
    	for itarg=1:N[1], jtarg=1:N[2]
	    m = max(abs(isrc-itarg),abs(jsrc-jtarg))
	    n = min(abs(isrc-itarg),abs(jsrc-jtarg))
	    Gw[itarg,jtarg] += G[m+1,n+1]cell[isrc,jsrc]
	end
    end
    Gw
end

struct Identity{T, M, N}
end

function Base.show(io::IO, c::Identity{T,M,N}) where {T,M,N}
    print(io, "Identity operator for field on a $M × $N grid")
end

function Identity(cell::Array{T,2}) where {T}
  m, n = size(cell)
  res = zeros(T,m,n)
  Identity{T,m,n}()
end

(c::Identity{T,M,N})(cell) where {T,M,N} = cell


"""
    Convolution{T, M, N}

A preplanned, zero-padded convolution operator for an M × N matrix.

# Fields
- `F̂`: preallocated space to store the DFT of the zero-padded input matrix
- `Ĝ`: DFT coefficients of the convolution kernel
- `fftop`: preplanned rFFT on a (2M-1) × (2N-1) matrix
- `paddedSpace`: preallocated space to zero-pad the input matrix

# Constructors:

- `Convolution(fftop::FFTW.rFFTWPlan{T}, Ĝ) where T`

# Example:
```jldoctest
julia> G = mirror(repmat(1:3,1,4))
5×7 Array{Int64,2}:
 3  3  3  3  3  3  3
 2  2  2  2  2  2  2
 1  1  1  1  1  1  1
 2  2  2  2  2  2  2
 3  3  3  3  3  3  3

julia> C = Convolution(plan_rfft(G), rfft(G))
Zero-padded convolution for a 3 × 4 matrix

julia> C(reshape(1:12, 3, 4))
3×4 Array{Float64,2}:
 164.0  164.0  164.0  164.0
 130.0  130.0  130.0  130.0
 148.0  148.0  148.0  148.0
```
"""
struct Convolution{T, M, N}
    F̂::Array{Complex{T}, 2}
    Ĝ::Array{Complex{T}, 2}
    fftop::FFTW.rFFTWPlan{T}
    paddedSpace::Array{T, 2}
end

function Base.show(io::IO, c::Convolution{T, M, N}) where {T, M, N}
    print(io, "Zero-padded convolution for a $M × $N matrix")
end

function Convolution(fftop::FFTW.rFFTWPlan{T}, Ĝ) where T
    m, n = size(fftop)
    #m̂, n̂ = @. ((m,n) + 1) ÷ 2
    m̂, n̂ = @. ((m,n) ÷ 2) + 1

    paddedSpace = zeros(T, m, n)

    F̂ = zeros(Complex{T}, m̂, n)

    Convolution{T, m̂, n̂}(F̂, Ĝ, fftop, paddedSpace)
end

function (c::Convolution{T,M,N})(cell) where {T,M,N}
    inds = CartesianRange((M,N))
    fill!(c.paddedSpace, zero(T))
    copy!(c.paddedSpace, inds, cell, inds)
    A_mul_B!(c.F̂, c.fftop, c.paddedSpace)

    c.F̂ .*= c.Ĝ

    A_ldiv_B!(c.paddedSpace, c.fftop, c.F̂)
    c.paddedSpace[M:end, N:end]
end

L⁻¹(g::Grid) = Convolution(g.fftop, g.lgfhat)
Q(g::Grid) = Convolution(g.fftop,g.qhat)
Id(g::DualPatch) = Identity(g.cell)

L⁻¹_slow(g::Grid,cell) = convolve(g.lgftab,cell)
Q_slow(g::Grid,cell) = convolve(g.qtab,cell)

"""
    mirror(A::AbstractMatrix{T}) where {T}

Mirror a matrix about its first column and row, i.e.
```
                    N … 1 … N
    1 … N        M |   |     |
 1 |     |       ⋮ |___|_____|
 ⋮ |  A  |  ⟶    1 |   |     |
 M |     |       ⋮ |   |  A  |
                 M |   |     |

```

Example

```jldoctest
julia> x = reshape(1:12, 3,4)
3×4 Base.ReshapedArray{Int64,2,UnitRange{Int64},Tuple{}}:
 1  4  7  10
 2  5  8  11
 3  6  9  12

julia> mirror(x)
5×7 Array{Int64,2}:
 12  9  6  3  6  9  12
 11  8  5  2  5  8  11
 10  7  4  1  4  7  10
 11  8  5  2  5  8  11
 12  9  6  3  6  9  12
```
"""
function mirror(a::AbstractArray{T,2}) where {T}
    Nr, Nc = size(a)
    A = zeros(T, 2Nr-1, 2Nc-1)
    A[1:Nr-1, 1:Nc-1] .= a[Nr:-1:2, Nc:-1:2]
    A[1:Nr-1, Nc:end] .= a[Nr:-1:2, 1:Nc]
    A[Nr:end, 1:Nc-1] .= a[1:Nr, Nc:-1:2]
    A[Nr:end, Nc:end] .= a
    A
end

function lgf_table!(g::Grid)
  g.lgftab = lgf(g)

  # Construct DFT of the table for use in fast application
  g.lgfhat = g.fftop * mirror(g.lgftab)

end

function q_table!(g::Grid,α::Float64)
  g.qtab = intfact(g,0.5α)
  g.α = α

  g.qhat = g.fftop * mirror(g.qtab)
end

function q_table!(g::Grid,α::Float64,c::Float64)
  g.qtab = intfact(g,c*α)
  g.α = α

  g.qhat = g.fftop * mirror(g.qtab)
end

quadgauss(f::Function) = dot(weights,f(nodes))


function Base.show(io::IO, g::Grid)
    println(io, "Grid: number of cells = ($(g.N[1]),$(g.N[2])), Δx = "*
    "$(g.Δx), xmin = $(g.xmin), xmax = $(g.xmax)")
end


end
