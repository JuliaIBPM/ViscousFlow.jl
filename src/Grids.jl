"""
The `Grids` module provides a data structure:

[`Grids.DualPatch`](@ref)

and functions that act on this structure.

"""
module Grids

export DualPatch

import Whirl2d
import Whirl2d:@get, MappedVector

using FastGaussQuadrature

abstract type Grid end

nodes, weights = gausslegendre(100)

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
    xfaceint::Array{UnitRange{Int},1}
    yfaceint::Array{UnitRange{Int},1}

    "mapping from grid index to matrix row/column"
    cellmap
    nodemap
    xfacemap
    yfacemap

    "vorticity field"
    w::Array{Float64,Whirl2d.ndim}

    "velocity field components"
    qx::Array{Float64,Whirl2d.ndim}
    qy::Array{Float64,Whirl2d.ndim}

    "dual velocity field components"
    dualqx::Array{Float64,Whirl2d.ndim}
    dualqy::Array{Float64,Whirl2d.ndim}

    "streamfunction"
    psi::Array{Float64,Whirl2d.ndim}

    "pressure"
    p::Array{Float64,Whirl2d.ndim}

    "LGF table"
    lgftab::Array{Float64,Whirl2d.ndim}
    lgfhat::Array{Complex{Float64},Whirl2d.ndim}

    "integrating factor table"
    qtab::Array{Float64,Whirl2d.ndim}
    qhat::Array{Complex{Float64},Whirl2d.ndim}

    "LGF*Q"
    gqhat::Array{Complex{Float64},Whirl2d.ndim}
    α::Float64

    "preplanned fft"
    fftop::FFTW.rFFTWPlan

end

function DualPatch(N,Δx,xmin)

    xmax = xmin + N*Δx

    # set up grid arrays with ghosts
    w = zeros(Float64,N[1]+2,N[2]+2)
    psi = zeros(Float64,N[1]+2,N[2]+2)
    p = zeros(Float64,N[1]+1,N[1]+1)
    qx = zeros(Float64,N[1]+2,N[2]+1)
    qy = zeros(Float64,N[1]+1,N[2]+2)
    dualqx = zeros(Float64,N[1]+1,N[2]+2)
    dualqy = zeros(Float64,N[1]+2,N[2]+1)

    # set the first interior indices in each direction
    # (this implicitly sets the number of ghost cell layers)
    ifirst = [2,2]

    # set the index ranges of interior data for each type of grid data array
    cellint =  ifirst-1+[1:N[1],1:N[2]];
    nodeint =  ifirst-1+[1:N[1]-1,1:N[2]-1];
    xfaceint = ifirst-1+[1:N[1],1:N[2]-1];
    yfaceint = ifirst-1+[1:N[1]-1,1:N[2]];

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
    xfacemap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*N[1]
    yfacemap(i,j) = i-ifirst[1]+1+(j-ifirst[2])*(N[1]-1)

    lgftab = Array{Float64}(0,0)
    lgfhat = Array{Complex{Float64}}(0,0)
    qtab = Array{Float64}(0,0)
    qhat = Array{Complex{Float64}}(0,0)
    gqhat = Array{Complex{Float64}}(0,0)
    α = 0.0

    fftop = FFTW.plan_rfft(zeros(2*N[1]+3,2*N[2]+3))

    DualPatch(N,Δx,xmin,xmax,ifirst,cellint,nodeint,xfaceint,yfaceint,
	      cellmap,nodemap,xfacemap,yfacemap,
	      w,qx,qy,dualqx,dualqy,psi,p,
        lgftab,lgfhat,qtab,qhat,gqhat,α,fftop)

end

# Differential operations
function curl!(w,ir::UnitRange{Int},jr::UnitRange{Int},qx,qy)
    for i=ir, j=jr
    	w[i,j] = -qx[i,j]+qx[i,j-1]+qy[i,j]-qy[i-1,j]
    end
end

function curl!(qx,qy,ir::UnitRange{Int},jr::UnitRange{Int},psi)
    for j=jr
    	qx[ir.stop+1,j] = psi[ir.stop+1,j+1]-psi[ir.stop+1,j]
    end
    for i=ir
    	qy[i,jr.stop+1] = psi[i,jr.stop+1]-psi[i+1,jr.stop+1]
    end
    for i=ir, j=jr
    	qx[i,j] = psi[i,j+1]-psi[i,j]
	qy[i,j] = psi[i,j]-psi[i+1,j]
    end
end

function diverg!(p,ir::UnitRange{Int},jr::UnitRange{Int},qx,qy)
    for i=ir, j=jr
    	p[i,j] = qx[i+1,j]-qx[i,j]+qy[i,j+1]-qy[i,j]
    end
end

function grad!(qx,qy,ir::UnitRange{Int},jr::UnitRange{Int},p)
    for i=ir, j=jr
    	qx[i,j] = p[i,j]-p[i-1,j]
	qy[i,j] = p[i,j]-p[i,j-1]
    end
    for j=jr
    	qy[ir.start-1,j] = p[ir.start-1,j]-p[ir.start-1,j-1]
    end
    for i=ir
    	qx[i,jr.start-1] = p[i,jr.start-1]-p[i-1,jr.start-1]
    end
end

function lap!(lapf,ir::UnitRange{Int},jr::UnitRange{Int},f)
    for i=ir, j=jr
    	lapf[i,j] = f[i+1,j]+f[i-1,j]+f[i,j+1]+f[i,j-1]-4f[i,j]
    end
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},qx,qy)
    for j=jr
    	vx[ir.start-1,j] = 0.25(qx[ir.start-1,j]+qx[ir.start,j]+
			        qx[ir.start-1,j-1]+qx[ir.start,j-1])
    end
    for i=ir
	vy[i,jr.start-1] = 0.25(qy[i-1,jr.start-1]+qy[i-1,jr.start]+
			        qy[i,jr.start-1]+qy[i,jr.start])
    end
    for i=ir, j=jr
    	vx[i,j] = 0.25(qx[i,j]+qx[i+1,j]+qx[i,j-1]+qx[i+1,j-1])
	vy[i,j] = 0.25(qy[i-1,j]+qy[i-1,j+1]+qy[i,j]+qy[i,j+1])
    end
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},w)
    for j=jr
    	vx[ir.start-1,j]=0.5(w[ir.start-1,j]+w[ir.start,j])
    end
    for i=ir
    	vy[i,jr.start-1]=0.5(w[i,jr.start-1]+w[i,jr.start])
    end
    for i=ir,j=jr
    	vx[i,j] = 0.5(w[i,j]+w[i+1,j])
    	vy[i,j] = 0.5(w[i,j]+w[i,j+1])
    end
end

# Differential operations with grid interface
function curl(g::DualPatch,qx,qy)
    w = zeros(g.w)
    curl!(w,g.cellint[1],g.cellint[2],qx,qy)
    w
end

function curl(g::DualPatch,psi)
    qx = zeros(g.qx)
    qy = zeros(g.qy)
    # also calculate curl in ghost cells
    curl!(qx,qy,g.cellint[1].start-1:g.cellint[1].stop,
	        g.cellint[2].start-1:g.cellint[2].stop,psi)
    qx, qy
end

function diverg(g::DualPatch,qx,qy)
    p = zeros(g.p)
    diverg!(p,g.cellint[1].start-1:g.cellint[1].stop,
	      g.cellint[2].start-1:g.cellint[2].stop,qx,qy)
    p
end

function grad(g::DualPatch,p)
    qx = zeros(g.qx)
    qy = zeros(g.qy)
    grad!(qx,qy,g.cellint[1],g.cellint[2],p)
    qx,qy
end

function lap(g::DualPatch,psi)
    lappsi = zeros(g.psi)
    lap!(lappsi,g.cellint[1],g.cellint[2],psi)
    lappsi
end

function lap(g::DualPatch,qx,qy)
    lapqx = zeros(g.qx)
    lapqy = zeros(g.qy)
    lap!(lapqx,g.xfaceint[1],g.xfaceint[2],qx)
    lap!(lapqy,g.yfaceint[1],g.yfaceint[2],qy)
    lapqx,lapqy
end

function shift(g::DualPatch,qx,qy)
    vx = zeros(g.dualqx)
    vy = zeros(g.dualqy)
    shift!(vx,vy,g.cellint[1],g.cellint[2],qx,qy)
    vx, vy
end

function shift(g::DualPatch,w)
    vx = zeros(g.dualqx)
    vy = zeros(g.dualqy)
    shift!(vx,vy,g.cellint[1],g.cellint[2],w)
    vx, vy
end

function cross(g::DualPatch,w,qx,qy)
    vx,vy = shift(g,qx,qy)
    ux,uy = shift(g,w)
    -uy.*vy, ux.*vx
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
lgf(g::Grid) = reshape([lgf([i,j]) for i=0:g.N[1]+1 for j=0:g.N[2]+1],
       g.N[1]+2,g.N[2]+2)

"""
    GE = intfact(g::Grid,a::Float64)

Set up a table of integrating factor values in the upper right
quadrant, centered at the lower left ghost cell in grid 'g'.
This zeros out the values lower than machine epsilon.
"""
intfact(g::Grid,a::Float64) = reshape([intfact([i,j],a) > eps(Float64) ?
        intfact([i,j],a) : 0.0
        for i=0:g.N[1]+1 for j=0:g.N[2]+1], g.N[1]+2,g.N[2]+2)


"""
    s = convolve(G,w)

Perform a discrete convolution of grid data w with one of the Green's
function tables (the LGF or the integrating factor). This exploits the
symmetries in these functions.
"""
function convolve(G::Array{T,2},w::AbstractArray{T,2}) where T
    N = size(w)

    Gw = zeros(w)
    for isrc=1:N[1],jsrc=1:N[2]
      if abs(w[isrc,jsrc])<eps(Float64)
	         continue
	    end
    	for itarg=1:N[1], jtarg=1:N[2]
	    m = max(abs(isrc-itarg),abs(jsrc-jtarg))
	    n = min(abs(isrc-itarg),abs(jsrc-jtarg))
	    Gw[itarg,jtarg] += G[m+1,n+1]w[isrc,jsrc]
	end
    end
    Gw
end

function convolve_fft(fftop::FFTW.FFTWPlan{T,K,false},
        Ghat::Array{Complex{T},2},w::AbstractArray{T,2}) where {T<:AbstractFloat,K}

  wbig = mirror(zeros(w))
  n = size(w)
  wbig[1:n[1],1:n[2]] = w
  what = fftop * wbig

  Gwbig = fftop \ (Ghat .* what)
  Gwbig[n[1]:end,n[2]:end]

end

convolve_fft(g::Grid,Ghat::Array{Complex{T},2},w) where T =
        convolve_fft(g.fftop,Ghat,w)


L⁻¹(g::Grid,w) = convolve_fft(g,g.lgfhat,w)
Q(g::Grid,w) = convolve_fft(g,g.qhat,w)

L⁻¹_slow(g::Grid,w) = convolve(g.lgftab,w)
Q_slow(g::Grid,w) = convolve(g.qtab,w)

mirror(a::AbstractArray{T,2} where T) = hcat(
  flipdim(vcat(flipdim(a[:,2:end],1),a[2:end,2:end]),2),
  vcat(flipdim(a,1),a[2:end,:])
  )

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

function quadgauss(f::Function)
    dot(weights,f(nodes))
end


function Base.show(io::IO, g::Grid)
    println(io, "Grid: number of cells = ($(g.N[1]),$(g.N[2])), Δx = "*
    "$(g.Δx), xmin = $(g.xmin), xmax = $(g.xmax)")
end


end
