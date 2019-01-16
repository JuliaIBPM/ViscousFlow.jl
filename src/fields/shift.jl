"""
    shift!(q::Edges{Dual},w::Nodes{Dual})

Shift (by linear interpolation) the dual nodal data `w` to the edges of the dual
cells, and return the result in `q`.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Dual,w);

julia> shift!(q,w)
Whirl.Fields.Edges{Whirl.Fields.Dual,8,6} data
u (in grid orientation):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.5  0.5  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
v (in grid orientation):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function shift!(dual::Edges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 1:NX-1
        dual.u[x,y] = (w[x,y] + w[x+1,y])/2
    end

    @inbounds for y in 1:NY-1, x in 2:NX-1
        dual.v[x,y] = (w[x,y] + w[x,y+1])/2
    end
    dual
end


shift(nodes::Nodes{Dual,NX,NY}) where {NX,NY} = shift!(Edges(Dual, nodes), nodes)

"""
    shift!((wx::Nodes,wy::Nodes),q::Edges)

Shift (by linear interpolation) the edge data `q` (of either dual or primal
type) to the dual or primal nodes, and return the result in `wx` and `wy`. `wx`
holds the shifted `q.u` data and `wy` the shifted `q.v` data.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> wx = Nodes(Dual,(8,6)); wy = deepcopy(wx);

julia> Fields.shift!((wx,wy),q);

julia> wx
Whirl.Fields.Nodes{Whirl.Fields.Dual,8,6} data
Printing in grid orientation (lower left is (1,1)):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> wy
Whirl.Fields.Nodes{Whirl.Fields.Dual,8,6} data
Printing in grid orientation (lower left is (1,1)):
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function shift!(dual::Tuple{Nodes{Dual, NX, NY},Nodes{Dual, NX, NY}}, w::Edges{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 2:NY-1, x in 1:NX
        dual[1][x,y] = (w.u[x,y-1] + w.u[x,y])/2
    end

    @inbounds for y in 1:NY, x in 2:NX-1
        dual[2][x,y] = (w.v[x-1,y] + w.v[x,y])/2
    end
    dual
end

function shift!(dual::Tuple{Nodes{Dual, NX, NY},Nodes{Dual, NX, NY}}, w::Edges{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY, x in 2:NX-1
        dual[1][x,y] = (w.u[x-1,y] + w.u[x,y])/2
    end

    @inbounds for y in 2:NY-1, x in 1:NX
        dual[2][x,y] = (w.v[x,y-1] + w.v[x,y])/2
    end
    dual
end

function shift!(primal::Tuple{Nodes{Primal, NX, NY},Nodes{Primal, NX, NY}}, w::Edges{Primal,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[1][x,y] = (w.u[x,y] + w.u[x+1,y])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[2][x,y] = (w.v[x,y] + w.v[x,y+1])/2
    end
    primal
end

function shift!(primal::Tuple{Nodes{Primal, NX, NY},Nodes{Primal, NX, NY}}, w::Edges{Dual,NX, NY}) where {NX, NY}
    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[1][x,y] = (w.u[x,y] + w.u[x,y+1])/2
    end

    @inbounds for y in 1:NY-1, x in 1:NX-1
        primal[2][x,y] = (w.v[x,y] + w.v[x+1,y])/2
    end
    primal
end


# I don't like this one. It is ambiguous what type of nodes are being shifted to.
nodeshift(edges::Edges{Primal,NX,NY}) where {NX,NY} = shift!((Nodes(Dual, (NX,NY)),Nodes(Dual, (NX,NY))),edges)
