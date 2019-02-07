"""
    cellshift!(q::Edges{Dual},w::Nodes{Dual})

Shift (by linear interpolation) the dual nodal data `w` to the edges of the dual
cells, and return the result in `q`.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Dual,w);

julia> cellshift!(q,w)
Edges{Dual,8,6} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.5  0.5  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function cellshift!(dual::Edges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 1:NX-1
        dual.u[x,y] = (w[x,y] + w[x+1,y])/2
    end

    @inbounds for y in 1:NY-1, x in 2:NX-1
        dual.v[x,y] = (w[x,y] + w[x,y+1])/2
    end
    dual
end


cellshift(nodes::Nodes{Dual,NX,NY}) where {NX,NY} = cellshift!(Edges(Dual, nodes), nodes)

"""
    cellshift!((wx::Nodes,wy::Nodes),q::Edges)

Shift (by linear interpolation) the edge data `q` (of either dual or primal
type) to the dual or primal nodes, and return the result in `wx` and `wy`. `wx`
holds the shifted `q.u` data and `wy` the shifted `q.v` data.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> wx = Nodes(Dual,(8,6)); wy = deepcopy(wx);

julia> Fields.cellshift!((wx,wy),q);

julia> wx
Nodes{Dual,8,6} data
Printing in grid orientation (lower left is (1,1))
6×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> wy
Nodes{Dual,8,6} data
Printing in grid orientation (lower left is (1,1))
6×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function cellshift!(dual::Tuple{Nodes{Dual, NX, NY},Nodes{Dual, NX, NY}}, w::Edges{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 1:NX
        dual[1][x,y] = (w.u[x,y-1] + w.u[x,y])/2
    end

    @inbounds for y in 1:NY, x in 2:NX-1
        dual[2][x,y] = (w.v[x-1,y] + w.v[x,y])/2
    end
    dual
end

function cellshift!(dual::Tuple{Nodes{Dual, NX, NY},Nodes{Dual, NX, NY}}, w::Edges{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY, x in 2:NX-1
        dual[1][x,y] = (w.u[x-1,y] + w.u[x,y])/2
    end

    @inbounds for y in 2:NY-1, x in 1:NX
        dual[2][x,y] = (w.v[x,y-1] + w.v[x,y])/2
    end
    dual
end

function cellshift!(primal::Tuple{Nodes{Primal, NX, NY},Nodes{Primal, NX, NY}}, w::Edges{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[1][x,y] = (w.u[x,y] + w.u[x+1,y])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[2][x,y] = (w.v[x,y] + w.v[x,y+1])/2
    end
    primal
end

function cellshift!(primal::Tuple{Nodes{Primal, NX, NY},Nodes{Primal, NX, NY}}, w::Edges{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[1][x,y] = (w.u[x,y] + w.u[x,y+1])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[2][x,y] = (w.v[x,y] + w.v[x+1,y])/2
    end
    primal
end

"""
    cellshift!(q::Edges{Primal},w::Nodes{Primal})

Shift (by linear interpolation) the primal nodal data `w` to the edges of the primal cells,
and return the result in `q`.
"""
function cellshift!(primal::Edges{Primal, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-1, x in 2:NX-1
        primal.u[x,y] = (w[x-1,y] + w[x,y])/2
    end

    @inbounds for y in 2:NY-1, x in 1:NX-1
        primal.v[x,y] = (w[x,y-1] + w[x,y])/2
    end
    primal
end


# I don't like this one. It is ambiguous what type of nodes are being shifted to.
nodeshift(edges::Edges{Primal,NX,NY}) where {NX,NY} = cellshift!((Nodes(Dual, (NX,NY)),Nodes(Dual, (NX,NY))),edges)
