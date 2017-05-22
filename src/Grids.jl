module Grids

export Grid
export grad,div,curl

import Whirl2d
import Whirl2d:@get, MappedVector

abstract type Mesh end
abstract type GridData2d <: AbstractArray{Float64,2} end


struct Grid <: Mesh
    "number of cells in each direction"
    N::Vector{Int}

    "uniform grid spacing (in all directions)"
    Δx::Float64

    "coordinates of lower left-hand corner"
    xmin::Vector{Float64}

end

XFaceData(g::Grid) = zeros(Float64,g.N[1]+1,g.N[2]);
YFaceData(g::Grid) = zeros(Float64,g.N[1],g.N[2]+1);
CellData(g::Grid) = zeros(Float64,g.N[1],g.N[2]);
NodeData(g::Grid) = zeros(Float64,g.N[1]+1,g.N[2]+1);

function curl(qx,qy)
    Nx = size(qy,1)
    Ny = size(qx,2)
    w = zeros(eltype(qx),Nx+1,Ny+1)
    for i=2:Nx, j=2:Ny
    	w[i,j] = -qx[i,j]+qx[i,j-1]+qy[i,j]-qy[i-1,j]
    end
    w
end

function curl(psi)
    Nx = size(psi,1)-1
    Ny = size(psi,2)-1
    qx = zeros(eltype(psi),Nx+1,Ny)
    qy = zeros(eltype(psi),Nx,Ny+1)
    
    for i=1:Nx, j=1:Ny
    	qx[i,j] = psi[i,j+1]-psi[i,j]
	qy[i,j] = psi[i,j]-psi[i+1,j]
    end
    for j=1:Ny
    	qx[Nx+1,j] = psi[Nx+1,j+1]-psi[Nx+1,j]
    end
    for i=1:Nx
    	qy[i,Ny+1] = psi[i,Ny+1]-psi[i+1,Ny+1]
    end 
    qx, qy
end

function div(qx,qy)
    Nx = size(qy,1)
    Ny = size(qx,2)
    w = zeros(eltype(qx),Nx,Ny)
    for i=1:Nx, j=1:Ny
    	w[i,j] = qx[i+1,j]-qx[i,j]+qy[i,j+1]-qy[i,j]
    end
    w       
end

function grad(psi)
    Nx = size(psi,1)
    Ny = size(psi,2)
    qx = zeros(eltype(psi),Nx+1,Ny)
    qy = zeros(eltype(psi),Nx,Ny+1)
    
    for i=2:Nx, j=2:Ny
    	qx[i,j] = psi[i,j]-psi[i-1,j]
	qy[i,j] = psi[i,j]-psi[i,j-1]
    end
    for j=2:Ny
    	qy[1,j] = psi[1,j]-psi[1,j-1]
    end
    for i=2:Nx
    	qx[i,1] = psi[i,1]-psi[i-1,1]
    end 
    qx, qy
end

function Base.show(io::IO, g::Grid)
    println(io, "Grid: number of cells = ($(g.N[1]),$(g.N[2])), Δx = "*
    "$(g.Δx), xmin = ($(g.xmin[1]), $(g.xmin[2])), xmax = "*
    "($(g.N[1]*g.Δx+g.xmin[1]),$(g.N[2]*g.Δx+g.xmin[2]))") 
end


end
