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
    facexint::Array{UnitRange{Int},1}
    faceyint::Array{UnitRange{Int},1}

    "mapping from grid index to matrix row/column"
    cellmap
    nodemap
    facexmap
    faceymap

    "cell center variable"
    cell::Array{Float64,Whirl2d.ndim}

    "face variable components"
    facex::Array{Float64,Whirl2d.ndim}
    facey::Array{Float64,Whirl2d.ndim}

    "dual velocity field components"
    dualfacex::Array{Float64,Whirl2d.ndim}
    dualfacey::Array{Float64,Whirl2d.ndim}

    "cell node variable"
    node::Array{Float64,Whirl2d.ndim}

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
    cell = zeros(Float64,N[1]+2,N[2]+2)
    node = zeros(Float64,N[1]+1,N[1]+1)
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

    lgftab = Array{Float64}(0,0)
    lgfhat = Array{Complex{Float64}}(0,0)
    qtab = Array{Float64}(0,0)
    qhat = Array{Complex{Float64}}(0,0)
    gqhat = Array{Complex{Float64}}(0,0)
    α = 0.0

    fftop = FFTW.plan_rfft(zeros(2*N[1]+3,2*N[2]+3))

    DualPatch(N,Δx,xmin,xmax,ifirst,cellint,nodeint,facexint,faceyint,
	      cellmap,nodemap,facexmap,faceymap,
	      cell,facex,facey,dualfacex,dualfacey,node,
        lgftab,lgfhat,qtab,qhat,gqhat,α,fftop)

end

# Differential operations
function curl!(cell,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    for i=ir, j=jr
    	cell[i,j] = -facex[i,j]+facex[i,j-1]+facey[i,j]-facey[i-1,j]
    end
end

function curl!(facex,facey,ir::UnitRange{Int},jr::UnitRange{Int},cell)
    for j=jr
    	facex[ir.stop+1,j] = cell[ir.stop+1,j+1]-cell[ir.stop+1,j]
    end
    for i=ir
    	facey[i,jr.stop+1] = cell[i,jr.stop+1]-cell[i+1,jr.stop+1]
    end
    for i=ir, j=jr
    	facex[i,j] = cell[i,j+1]-cell[i,j]
	    facey[i,j] = cell[i,j]-cell[i+1,j]
    end
end

function diverg!(node,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    for i=ir, j=jr
    	node[i,j] = facex[i+1,j]-facex[i,j]+facey[i,j+1]-facey[i,j]
    end
end

function grad!(facex,facey,ir::UnitRange{Int},jr::UnitRange{Int},node)
    for i=ir, j=jr
    	facex[i,j] = node[i,j]-node[i-1,j]
	    facey[i,j] = node[i,j]-node[i,j-1]
    end
    for j=jr
    	facey[ir.start-1,j] = node[ir.start-1,j]-node[ir.start-1,j-1]
    end
    for i=ir
    	facex[i,jr.start-1] = node[i,jr.start-1]-node[i-1,jr.start-1]
    end
end

function lap!(lapf,ir::UnitRange{Int},jr::UnitRange{Int},f)
    for i=ir, j=jr
    	lapf[i,j] = f[i+1,j]+f[i-1,j]+f[i,j+1]+f[i,j-1]-4f[i,j]
    end
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},facex,facey)
    for j=jr
    	vx[ir.start-1,j] = 0.25(facex[ir.start-1,j]+facex[ir.start,j]+
			        facex[ir.start-1,j-1]+facex[ir.start,j-1])
    end
    for i=ir
	    vy[i,jr.start-1] = 0.25(facey[i-1,jr.start-1]+facey[i-1,jr.start]+
			        facey[i,jr.start-1]+facey[i,jr.start])
    end
    for i=ir, j=jr
    	vx[i,j] = 0.25(facex[i,j]+facex[i+1,j]+facex[i,j-1]+facex[i+1,j-1])
	    vy[i,j] = 0.25(facey[i-1,j]+facey[i-1,j+1]+facey[i,j]+facey[i,j+1])
    end
end

function shift!(vx,vy,ir::UnitRange{Int},jr::UnitRange{Int},cell)
    for j=jr
    	vx[ir.start-1,j]=0.5(cell[ir.start-1,j]+cell[ir.start,j])
    end
    for i=ir
    	vy[i,jr.start-1]=0.5(cell[i,jr.start-1]+cell[i,jr.start])
    end
    for i=ir,j=jr
    	vx[i,j] = 0.5(cell[i,j]+cell[i+1,j])
    	vy[i,j] = 0.5(cell[i,j]+cell[i,j+1])
    end
end

# Differential operations with grid interface
function curl(g::DualPatch,facex,facey)
    cell = zeros(g.cell)
    curl!(cell,g.cellint[1],g.cellint[2],facex,facey)
    cell
end

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

function grad(g::DualPatch,node)
    facex = zeros(g.facex)
    facey = zeros(g.facey)
    grad!(facex,facey,g.cellint[1],g.cellint[2],node)
    facex,facey
end

function lap(g::DualPatch,cell)
    lapcell = zeros(g.cell)
    lap!(lappsi,g.cellint[1],g.cellint[2],cell)
    lapcell
end

function lap(g::DualPatch,facex,facey)
    lapfacex = zeros(g.facex)
    lapfacey = zeros(g.facey)
    lap!(lapfacex,g.facexint[1],g.facexint[2],facex)
    lap!(lapfacey,g.faceyint[1],g.faceyint[2],facey)
    lapfacex,lapfacey
end

function shift(g::DualPatch,facex,facey)
    vx = zeros(g.dualfacex)
    vy = zeros(g.dualfacey)
    shift!(vx,vy,g.cellint[1],g.cellint[2],facex,facey)
    vx, vy
end

function shift(g::DualPatch,cell)
    vx = zeros(g.dualfacex)
    vy = zeros(g.dualfacey)
    shift!(vx,vy,g.cellint[1],g.cellint[2],cell)
    vx, vy
end

function cross(g::DualPatch,cell,facex,facey)
    vx,vy = shift(g,facex,facey)
    ux,uy = shift(g,cell)
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

function convolve_fft(fftop::FFTW.FFTWPlan{T,K,false},
        Ghat::Array{Complex{T},2},cell::AbstractArray{T,2}) where {T<:AbstractFloat,K}

  cellbig = mirror(zeros(cell))
  n = size(cell)
  cellbig[1:n[1],1:n[2]] = cell
  cellhat = fftop * cellbig

  Gwbig = fftop \ (Ghat .* cellhat)
  Gwbig[n[1]:end,n[2]:end]

end

convolve_fft(g::Grid,Ghat::Array{Complex{T},2},cell) where T =
        convolve_fft(g.fftop,Ghat,cell)


L⁻¹(g::Grid,cell) = convolve_fft(g,g.lgfhat,cell)
Q(g::Grid,cell) = convolve_fft(g,g.qhat,cell)

L⁻¹_slow(g::Grid,cell) = convolve(g.lgftab,cell)
Q_slow(g::Grid,cell) = convolve(g.qtab,cell)

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
