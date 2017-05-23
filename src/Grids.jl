module Grids

export DualPatch

import Whirl2d
import Whirl2d:@get, MappedVector

using FastGaussQuadrature

abstract type Grid end
const ndim = 2

nodes, weights = gausslegendre(100)

struct DualPatch <: Grid
    "array of number of dual cells in each direction in patch interior"
    N::Vector{Int}

    "uniform grid spacing (in all directions)"
    Δx::Float64

    "coordinates of lower left-hand corner"
    xmin::Vector{Float64}

    "first interior cell indices in each direction"
    ifirst::Vector{Int}

    "vorticity field"
    w::Array{Float64,ndim}

    "velocity field components"
    qx::Array{Float64,ndim}
    qy::Array{Float64,ndim}

    "streamfunction"
    psi::Array{Float64,ndim}    
    
    "preplanned DST used to solve Poisson equation"
    dst!::FFTW.r2rFFTWPlan

end

function DualPatch(N,Δx,xmin)

    # set up grid arrays with ghosts
    w = zeros(Float64,N[1]+2,N[2]+2)
    psi = zeros(Float64,N[1]+2,N[2]+2)
    qx = zeros(Float64,N[1]+2,N[2]+1)
    qy = zeros(Float64,N[1]+1,N[2]+2)

    # set the first interior indices in each direction
    ifirst = [2,2]

    dst! = FFTW.plan_r2r!(w, FFTW.RODFT00);

    DualPatch(N,Δx,xmin,ifirst,w,qx,qy,psi,dst!)

end

# Differential operations
function curl!(w,ix::NTuple{2,Int},iy::NTuple{2,Int},qx,qy)
    for i=ix[1]:ix[2], j=iy[1]:iy[2]
    	w[i,j] = -qx[i,j]+qx[i,j-1]+qy[i,j]-qy[i-1,j]
    end
end

function curl!(qx,qy,ix::NTuple{2,Int},iy::NTuple{2,Int},psi)
    for j=iy[1]:iy[2]
    	qx[ix[2]+1,j] = psi[ix[2]+1,j+1]-psi[ix[2]+1,j]
    end
    for i=ix[1]:ix[2]
    	qy[i,iy[2]+1] = psi[i,iy[2]+1]-psi[i+1,iy[2]+1]
    end
    for i=ix[1]:ix[2], j=iy[1]:iy[2]
    	qx[i,j] = psi[i,j+1]-psi[i,j]
	qy[i,j] = psi[i,j]-psi[i+1,j]
    end
end

function diverg!(p,ix::NTuple{2,Int},iy::NTuple{2,Int},qx,qy)
    for i=ix[1]:ix[2], j=iy[1]:iy[2]
    	p[i,j] = qx[i+1,j]-qx[i,j]+qy[i,j+1]-qy[i,j]
    end      
end

function grad!(qx,qy,ix::NTuple{2,Int},iy::NTuple{2,Int},p)
    for i=ix[1]:ix[2], j=iy[1]:iy[2]
    	qx[i,j] = p[i,j]-p[i-1,j]
	qy[i,j] = p[i,j]-p[i,j-1]
    end
    for j=iy[1]:iy[2]
    	qy[ix[1]-1,j] = p[ix[1]-1,j]-p[ix[1]-1,j-1]
    end
    for i=ix[1]:ix[2]
    	qx[i,iy[1]-1] = p[i,iy[1]-1]-p[i-1,iy[1]-1]
    end 
end

function lap!(lapf,ix::NTuple{2,Int},iy::NTuple{2,Int},f)
    for i=ix[1]:ix[2], j=iy[1]:iy[2]
    	lapf[i,j] = f[i+1,j]+f[i-1,j]+f[i,j+1]+f[i,j-1]-4f[i,j]
    end
end



# Differential operations with grid interface
function curl(g::DualPatch,qx,qy)
    w = zeros(g.w)
    curl!(w,(g.ifirst[1],g.ifirst[1]+g.N[1]-1),
	    (g.ifirst[2],g.ifirst[2]+g.N[2]-1),qx,qy)
    w
end

function curl(g::DualPatch,psi)
    qx = zeros(g.qx)
    qy = zeros(g.qy)
    curl!(qx,qy,(g.ifirst[1]-1,g.ifirst[1]+g.N[1]-1),
	        (g.ifirst[2]-1,g.ifirst[2]+g.N[2]-1),psi)
    qx, qy
end

function diverg(g::DualPatch,qx,qy)
    p = zeros(Float64,g.N[1]+1,g.N[2]+1)
    diverg!(p,(g.ifirst[1]-1,g.ifirst[1]+g.N[1]-1),
	      (g.ifirst[2]-1,g.ifirst[2]+g.N[2]-1),qx,qy)
    p
end

function grad(g::DualPatch,p)
    qx = zeros(g.qx)
    qy = zeros(g.qy)
    grad!(qx,qy,(g.ifirst[1],g.ifirst[1]+g.N[1]-1),
	        (g.ifirst[2],g.ifirst[2]+g.N[2]-1),p)
    qx,qy
end

function lap(g::DualPatch,psi)
    lappsi = zeros(g.psi)
    lap!(lappsi,(g.ifirst[1],g.ifirst[1]+g.N[1]-1),
	        (g.ifirst[2],g.ifirst[2]+g.N[2]-1),psi)
    lappsi
end

function lap(g::DualPatch,qx,qy)
    lapqx = zeros(g.qx)
    lapqy = zeros(g.qy)
    lap!(lapqx,(g.ifirst[1],g.ifirst[1]+g.N[1]-1),
	       (g.ifirst[2],g.ifirst[2]+g.N[2]-2),qx)
    lap!(lapqy,(g.ifirst[1],g.ifirst[1]+g.N[1]-2),
	       (g.ifirst[2],g.ifirst[2]+g.N[2]-1),qy)
    lapqx,lapqy
end


# Lattice Green's function
function lgf(n)

    if n[1]==n[2]==0
       return 0.0
    else
       v = quadgauss() do x
        if x == -1
	   return sqrt(2)abs(n[1])
        else
           t = (x+1)/2
	   return 0.5real((1 - ( (t-sqrt(1im))./(t+sqrt(1im)) ).^(n[2]+abs(n[1])) .* 
     	     ( (t+sqrt(-1im))./(t-sqrt(-1im)) ).^(n[2]-abs(n[1])) )) ./t
	end
	
       end
       return 0.5v/pi
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
"""
intfact(g::Grid,a::Float64) = reshape([intfact([i,j],a) 
        for i=0:g.N[1]+1 for j=0:g.N[2]+1], g.N[1]+2,g.N[2]+2)



function lap_inv(w,G)
    N = size(w)

    s = zeros(w)
    for itarg=1:N[1], jtarg=1:N[2]
    	for isrc=1:N[1],jsrc=1:N[2]
	    s[itarg,jtarg] += G[abs(isrc-itarg)+1,abs(jsrc-jtarg)+1]w[isrc,jsrc]
	end
    end
    s
end

function quadgauss(f::Function)
    dot(weights,f(nodes))
end


function Base.show(io::IO, g::Grid)
    println(io, "Grid: number of cells = ($(g.N[1]),$(g.N[2])), Δx = "*
    "$(g.Δx), xmin = ($(g.xmin[1]), $(g.xmin[2])), xmax = "*
    "($(g.N[1]*g.Δx+g.xmin[1]),$(g.N[2]*g.Δx+g.xmin[2]))") 
end


end
