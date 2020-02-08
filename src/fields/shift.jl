# Grid interpolation operations

# Load the 1-d interpolations on which most of below are based.
include("interpolation1d.jl")


"""
    grid_interpolate!(q::Edges{Dual},w::Nodes{Dual})

Interpolate the dual nodal data `w` to the edges of the dual
cells, and return the result in `q`.

# Example

```jldoctest
julia> w = Nodes(Dual,(8,6));

julia> w[3,4] = 1.0;

julia> q = Edges(Dual,w);

julia> grid_interpolate!(q,w)
Edges{Dual,8,6,Float64} data
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
function grid_interpolate!(q::Edges{Dual, NX, NY}, w::Nodes{Dual,NX, NY}) where {NX, NY}
    grid_interpolate!(q.u,w)
    grid_interpolate!(q.v,w)
    q
end

# This operation is not necessarily desirable, since its meaning is ambiguous
grid_interpolate(nodes::Nodes{Dual,NX,NY}) where {NX,NY} = grid_interpolate!(Edges(Dual, nodes), nodes)




"""
    grid_interpolate!((wx::Nodes,wy::Nodes),q::Edges)

Interpolate the edge data `q` (of either dual or primal
type) to the dual or primal nodes, and return the result in `wx` and `wy`. `wx`
holds the shifted `q.u` data and `wy` the shifted `q.v` data.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> wx = Nodes(Dual,(8,6)); wy = deepcopy(wx);

julia> Fields.grid_interpolate!((wx,wy),q);

julia> wx
Nodes{Dual,8,6,Float64} data
Printing in grid orientation (lower left is (1,1))
6×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.5  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> wy
Nodes{Dual,8,6,Float64} data
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
function grid_interpolate!(out::Tuple{Nodes{C, NX, NY},Nodes{C, NX, NY}}, q::Edges{D,NX, NY}) where {C<:CellType, D<:CellType, NX, NY}
    grid_interpolate!(out[1],q.u)
    grid_interpolate!(out[2],q.v)
    out
end

"""
    grid_interpolate!(q::Edges{Primal},w::Nodes{Primal})

Interpolate the primal nodal data `w` to the edges of the primal cells,
and return the result in `q`.
"""
function grid_interpolate!(q::Edges{Primal, NX, NY}, w::Nodes{Primal,NX, NY}) where {NX, NY}
    grid_interpolate!(q.u,w)
    grid_interpolate!(q.v,w)
    q
end

# (Dual/Primal) edges to (Primal/Dual) edges. These require some rethinking,
# since they lead to broadened stencils.

"""
    grid_interpolate!(v::Edges{Dual/Primal},q::Edges{Primal/Dual})

Interpolate the primal (resp. dual) edge data `q` to the
edges of the dual (resp. primal) cells, and return the result in `v`.

# Example

```jldoctest
julia> q = Edges(Primal,(8,6));

julia> q.u[3,2] = 1.0;

julia> v = Edges(Dual,(8,6));

julia> Fields.grid_interpolate!(v,q)
Edges{Dual,8,6,Float64} data
u (in grid orientation)
6×7 Array{Float64,2}:
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
 0.0  0.25  0.25  0.0  0.0  0.0  0.0
 0.0  0.25  0.25  0.0  0.0  0.0  0.0
 0.0  0.0   0.0   0.0  0.0  0.0  0.0
v (in grid orientation)
5×8 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
function grid_interpolate!(dual::Edges{Dual, NX, NY},
                primal::Edges{Primal, NX, NY}) where {NX, NY}
    uₚ = primal.u
    @inbounds for y in 2:NY-1, x in 1:NX-1
        dual.u[x,y] = (uₚ[x,y] + uₚ[x+1,y] + uₚ[x,y-1] + uₚ[x+1,y-1])/4
    end

    vₚ = primal.v
    @inbounds for y in 1:NY-1, x in 2:NX-1
        dual.v[x,y] = (vₚ[x,y] + vₚ[x-1,y] + vₚ[x,y+1] + vₚ[x-1,y+1])/4
    end
    dual
end

function grid_interpolate(primal::Edges{Primal, NX, NY}) where {NX, NY}
    grid_interpolate!(Edges(Dual, (NX, NY)), primal)
end

function grid_interpolate!(primal::Edges{Primal, NX, NY},
                dual::Edges{Dual, NX, NY}) where {NX, NY}
    uₚ = dual.u
    @inbounds for y in 1:NY-1, x in 2:NX-1
        primal.u[x,y] = (uₚ[x,y] + uₚ[x-1,y] + uₚ[x,y+1] + uₚ[x-1,y+1])/4
    end

    vₚ = dual.v
    @inbounds for y in 2:NY-1, x in 1:NX-1
        primal.v[x,y] = (vₚ[x,y] + vₚ[x+1,y] + vₚ[x,y-1] + vₚ[x+1,y-1])/4
    end
    primal
end

function grid_interpolate(dual::Edges{Dual, NX, NY}) where {NX, NY}
    grid_interpolate!(Edges(Primal, dual), dual)
end




# I don't like this one. It is ambiguous what type of nodes are being shifted to.
nodeshift(edges::Edges{Primal,NX,NY}) where {NX,NY} = grid_interpolate!((Nodes(Dual, edges),Nodes(Dual, edges)),edges)
