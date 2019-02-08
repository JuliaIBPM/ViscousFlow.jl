
export PhysicalGrid, limits, origin, cellsize

struct PhysicalGrid{ND}
  N :: NTuple{ND,Int}
  I0 :: NTuple{ND,Int}
  Δx :: Float64
  xlim :: NTuple{ND,Tuple{Real,Real}}
end

"""
    PhysicalGrid(xlim::Tuple{Real,Real},ylim::Tuple{Real,Real},Δx::Float64)

Constructor to set up a grid connected to physical space. The region to be
discretized by the grid is defined by the limits `xlim` and `ylim`, and the
cell spacing (uniform and indentical in each direction) is specified by `Δx`.
The constructor uses this information to determine the number of
cells in each direction, expanding the given range if necessary to accommodate
an integer number. It also pads each side with a ghost cell.
It also determines the indices corresponding to the corner
of the cell to which the physical origin corresponds. Note that the corner
corresponding to the lowest limit in each direction has indices (1,1).
"""
function PhysicalGrid(xlim::Tuple{Real,Real},
                      ylim::Tuple{Real,Real},Δx::Float64)


  #= set grid spacing and the grid position of the origin
  In case the physical limits are not consistent with an integer number of dual cells, based on
  the given Δx, we adjust them outward a bit in all directions. We also seek to place the
  origin on the corner of a cell.
  =#
  xmin, xmax = xlim
  ymin, ymax = ylim
  @assert xmax >= xmin && ymax >= ymin "Maximum limits must exceed minimum limits"
  Lx = xmax-xmin
  Ly = ymax-ymin

  NX, i0, xlimnew = set_1d_grid(xmin,xmax,Δx)
  NY, j0, ylimnew = set_1d_grid(ymin,ymax,Δx)

  PhysicalGrid((NX,NY),(i0,j0),Δx,(xlimnew,ylimnew))
end

function set_1d_grid(xmin::Real,xmax::Real,Δx::Float64)
  NL, NR = floor(Int,xmin/Δx), ceil(Int,xmax/Δx)
  return NR-NL+2, 1-NL, (Δx*NL, Δx*NR)
end

"""
    size(g::PhysicalGrid,d::Int) -> Int

Return the number of cells in direction `d` in grid `g`.
"""
Base.size(g::PhysicalGrid,d::Int) = g.N[d]

"""
    size(g::PhysicalGrid) -> Tuple

Return a tuple of the number of cells in all directions in grid `g`.
"""
Base.size(g::PhysicalGrid) = g.N

"""
    length(g::PhysicalGrid,d::Int) -> Int

Return the total number of cells in grid `g`.
"""
Base.length(g::PhysicalGrid) = prod(size(g))

"""
    limits(g::PhysicalGrid,d::Int) -> Tuple

Return the minimum and maximum physical dimensions in direction `d` for grid `g`.
"""
limits(g::PhysicalGrid,d::Int) = g.xlim[d]

"""
    coordinates(w::Nodes/Edges,g::PhysicalGrid) -> Range

Return coordinate data range for type of `w`.
"""
coordinates(w,g::PhysicalGrid) = coordinates(w,dx=g.Δx,I0=g.I0)

"""
    origin(g::PhysicalGrid) -> Tuple{Int,Int}

Return a tuple of the indices of the primal node that corresponds to the
physical origin of the coordinate system used by `g`. Note that these
indices need not lie inside the range of indices occupied by the grid.
For example, if the range of physical coordinates occupied by the grid
is (1.0,3.0) x (2.0,4.0), then the origin is not inside the grid.
"""
origin(g::PhysicalGrid) = g.I0

"""
    cellsize(g::PhysicalGrid) -> Float64

Return the grid cell size of system `sys`
"""
cellsize(g::PhysicalGrid) = g.Δx
