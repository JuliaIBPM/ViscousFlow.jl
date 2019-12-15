
# These routines interpolate only within the grid, from cell centers to cell
# edges.

# Interpolations in one direction only, C -> E or E -> C

"""
    grid_interpolate!(out:ScalarGridData,in::ScalarGridData) -> ScalarGridData

Return the 1-d symmetric interpolation of scalar grid data `in` onto scalar grid
data `out`. Either `in` or `out` must be edge component data and the other must
be node data. The direction of interpolation is determined by the relationship
of `in` and `out`. For example, if `in` is dual nodes (cell by cell) and `out` is
primal x-edge components (cell by edge), then the interpolation takes place in the
y direction, since they are different types in this direction. (In the x direction,
they are of the same type (cell), so there is no interpolation in that direction.)
"""
function grid_interpolate! end


# Nodes to edge components

function grid_interpolate!(qu::XEdges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    # E x C <- C x C
    @inbounds for y in 2:NY-1, x in 1:NX-1
        qu[x,y] = (w[x,y] + w[x+1,y])/2
    end
    qu
end

function grid_interpolate!(qv::YEdges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    # C x E <- C x C
    @inbounds for y in 1:NY-1, x in 2:NX-1
        qv[x,y] = (w[x,y] + w[x,y+1])/2
    end
    qv
end

function grid_interpolate!(qu::XEdges{Primal, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    # C x E <- C x C
    @inbounds for y in 1:NY-1, x in 2:NX-1
        qu[x,y] = (w[x,y] + w[x,y+1])/2
    end
    qu
end

function grid_interpolate!(qv::YEdges{Primal, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
   # E x C <- C x C
    @inbounds for y in 2:NY-1, x in 1:NX-1
        qv[x,y] = (w[x,y] + w[x+1,y])/2
    end
    qv
end

function grid_interpolate!(qu::XEdges{Primal, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    # C x E <- E x E
    @inbounds for y in 1:NY-1, x in 2:NX-1
      qu[x,y] = (w[x-1,y] + w[x,y])/2
    end
    qu
end

function grid_interpolate!(qv::YEdges{Primal, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    # E x C <- E x E
    @inbounds for y in 2:NY-1, x in 1:NX-1
      qv[x,y] = (w[x,y-1] + w[x,y])/2
    end
    qv
end

function grid_interpolate!(qu::XEdges{Dual, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    # E x C <- E x E
    @inbounds for y in 2:NY-1, x in 1:NX-1
      qu[x,y] = (w[x,y-1] + w[x,y])/2
    end
    qu
end

function grid_interpolate!(qv::YEdges{Dual, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    # C x E <- E x E
    @inbounds for y in 1:NY-1, x in 2:NX-1
      qv[x,y] = (w[x-1,y] + w[x,y])/2
    end
    qv
end

# Edge components to nodes

function grid_interpolate!(w::Nodes{Dual, NX, NY}, qu::XEdges{Primal,NX, NY}) where {NX, NY}
    # C x C <- C x E
    @inbounds for y in 2:NY-1, x in 1:NX
        w[x,y] = (qu[x,y-1] + qu[x,y])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Dual, NX, NY}, qv::YEdges{Primal,NX, NY}) where {NX, NY}
    # C x C <- E x C
    @inbounds for y in 1:NY, x in 2:NX-1
        w[x,y] = (qv[x-1,y] + qv[x,y])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Dual, NX, NY}, qu::XEdges{Dual,NX, NY}) where {NX, NY}
  # C x C <- E x C
  @inbounds for y in 1:NY, x in 2:NX-1
      w[x,y] = (qu[x-1,y] + qu[x,y])/2
  end
    w
end

function grid_interpolate!(w::Nodes{Dual, NX, NY}, qv::YEdges{Dual,NX, NY}) where {NX, NY}
    # C x C <- C x E
    @inbounds for y in 2:NY-1, x in 1:NX
      w[x,y] = (qv[x,y-1] + qv[x,y])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Primal, NX, NY}, qu::XEdges{Dual,NX, NY}) where {NX, NY}
    # E x E <- E x C
    @inbounds for y in 1:NY-1, x in 1:NX-1
      w[x,y] = (qu[x,y] + qu[x,y+1])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Primal, NX, NY}, qv::YEdges{Dual,NX, NY}) where {NX, NY}
    # E x E <- C x E
    @inbounds for y in 1:NY-1, x in 1:NX-1
      w[x,y] = (qv[x,y] + qv[x+1,y])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Primal, NX, NY}, qu::XEdges{Primal,NX, NY}) where {NX, NY}
    # E x E <- C x E
    @inbounds for y in 1:NY-1, x in 1:NX-1
      w[x,y] = (qu[x,y] + qu[x+1,y])/2
    end
    w
end

function grid_interpolate!(w::Nodes{Primal, NX, NY}, qv::YEdges{Primal,NX, NY}) where {NX, NY}
    # E x E <- E x C
    @inbounds for y in 1:NY-1, x in 1:NX-1
      w[x,y] = (qv[x,y] + qv[x,y+1])/2
    end
    w
end
