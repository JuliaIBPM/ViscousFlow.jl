module Grids

export Grid

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

type XFaceData <: GridData2d end
type YFaceData <: GridData2d end
type CellData <: GridData2d end
type NodeData <: GridData2d end

function GridData(g::Grid,n::NTuple{2,Int})
    zeros(Float64,g.N[1]+n[1],g.N[2]+n[2]);
end

XFaceData(g::Grid) = GridData(g,(1,0));
YFaceData(g::Grid) = GridData(g,(0,1));
CellData(g::Grid) = GridData(g,(0,0));
NodeData(g::Grid) = GridData(g,(1,1));



function Base.show(io::IO, g::Grid)
    println(io, "Grid: number of cells = ($(g.N[1]),$(g.N[2])), Δx = "*
    "$(g.Δx), xmin = ($(g.xmin[1]), $(g.xmin[2])), xmax = "*
    "($(g.N[1]*g.Δx+g.xmin[1]),$(g.N[2]*g.Δx+g.xmin[2]))") 
end


end
